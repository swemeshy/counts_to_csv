use anyhow::anyhow;
use clap::Clap;
use csv;
use hdf5::types::*;
use indicatif::{ProgressBar, ProgressIterator, ProgressStyle};
use log::info;
use pyo3::prelude::*;
use pyo3::types::PyString;
use serde::Serialize;
use sprs::CsMatBase;
use std::path::PathBuf;

#[pymodule]
fn counts_to_csv(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    #[pyfn(m, "make_csv")]
    fn make_csv_py(
        adata: &PyAny,
        delimiter: &PyString,
        column_orient: &PyString,
        outfile: &PyString,
    ) -> PyResult<()> {
        use pyo3::exceptions::PyException;
        let outfile_path = PathBuf::from(outfile.to_str()?);
        let orient = {
            let x = match column_orient.to_str()? {
                "var-names" => Ok(Orient::VarNames),
                "obs-names" => Ok(Orient::ObsNames),
                _ => Err(anyhow::anyhow!(
                    "Invalid value: {}\nPossible values: obs-names, var-names",
                    column_orient
                )),
            };
            x.unwrap()
        };
        let delim = {
            let x = match delimiter.to_str()? {
                "comma" => Ok(Delimiter::Comma),
                "tab" => Ok(Delimiter::Tab),
                "colon" => Ok(Delimiter::Colon),
                "pipe" => Ok(Delimiter::Pipe),
                "semicolon" => Ok(Delimiter::Semicolon),
                _ => Err(anyhow::anyhow!(
                    "Invalid value: {}\nPossible values: comma, tab, colon, pipe, semicolon",
                    delimiter
                )),
            };
            x.unwrap()
        };
        let libopts = LibOpts {
            column_orient: orient,
            delimiter: delim,
            outfile: outfile_path,
        };
        make_csv(adata, libopts).map_err(|err| PyException::new_err(err.to_string()))
    }

    Ok(())
}

// for argument parsing
struct LibOpts {
    column_orient: Orient,
    delimiter: Delimiter,
    outfile: PathBuf,
}

/// argument enum for delimiter
#[derive(Clap)]
pub(crate) enum Delimiter {
    Comma,
    Tab,
    Colon,
    Pipe,
    Semicolon,
}

/// argument enum for column_orient
#[derive(Clap)]
pub(crate) enum Orient {
    VarNames,
    ObsNames,
}

/// trait that describes the type for the data array
pub(crate) trait ArrayDtype: H5Type + Default + Copy + Serialize {}
impl<T> ArrayDtype for T where T: H5Type + Default + Copy + Serialize {}

/// represents a row to be written to the CSV file
#[derive(Serialize)]
pub(crate) struct Row<'a, T: ArrayDtype> {
    pub(crate) name: &'a str,
    pub(crate) values: RowValIter<'a, T>,
}

/// an iterator wrapper for the indices and data of a CsVecBase object that represents one row
/// in the counts matrix
#[derive(Clone)]
pub(crate) struct RowValIter<'a, T> {
    data: std::slice::Iter<'a, T>,
    indices: std::slice::Iter<'a, usize>,
    next: Option<usize>,
    index: usize,
    stop: usize,
}

/// to create a new instance of the RowValIter iterator wrapper
impl<'a, T: ArrayDtype> RowValIter<'a, T> {
    pub(crate) fn new(values: &'a sprs::vec::CsVecBase<&[usize], &[T]>) -> RowValIter<'a, T> {
        let data_iter = values.data().iter();
        let mut idx_iter = values.indices().iter();
        let next = idx_iter.next().cloned();
        let stop = values.dim();
        RowValIter {
            data: data_iter,
            indices: idx_iter,
            next,
            index: 0,
            stop,
        }
    }
}

/// iterator implementation for RowValIter
impl<'a, T: ArrayDtype> Iterator for RowValIter<'a, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next {
            Some(ind) => {
                if self.index == ind {
                    self.next = self.indices.next().cloned();
                    self.index += 1;
                    self.data.next().cloned()
                } else {
                    self.index += 1;
                    Some(Default::default())
                }
            }
            None => {
                if self.index < self.stop {
                    self.index += 1;
                    Some(Default::default())
                } else {
                    None
                }
            }
        }
    }
}

/// serialization implementation for RowValIter
impl<'a, T: ArrayDtype> Serialize for RowValIter<'a, T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.collect_seq(self.clone())
    }
}

/// custom progress bar that shows elapsed time and percentage completion
pub(crate) fn create_progress_bar(iter_size: usize) -> ProgressBar {
    let bar = ProgressBar::new(iter_size as u64);
    bar.set_style(ProgressStyle::default_bar().template("[{elapsed}] {wide_bar} {percent}%"));
    bar
}

/// write counts matrix to csv
fn arrays_to_csv<T: ArrayDtype>(
    args: LibOpts,
    data: Vec<T>,
    indptr: Vec<usize>,
    indices: Vec<usize>,
    obs_vec: Vec<String>,
    var_vec: Vec<String>,
) -> anyhow::Result<()> {
    // get delimiter from args.delimiter
    let delimiter = match args.delimiter {
        Delimiter::Comma => b',',
        Delimiter::Tab => b'\t',
        Delimiter::Colon => b':',
        Delimiter::Pipe => b'|',
        Delimiter::Semicolon => b';',
    };

    // construct sparse matrix
    let mut counts_mtx = CsMatBase::try_new((obs_vec.len(), var_vec.len()), indptr, indices, data)?;

    // transpose matrix if needed based on column orientation specified
    let (header, first_col, row_names) = match args.column_orient {
        Orient::ObsNames => {
            counts_mtx.transpose_mut();
            counts_mtx = counts_mtx.to_csr();
            (obs_vec, "gene", var_vec)
        }
        Orient::VarNames => (var_vec, "cell", obs_vec),
    };

    // open CSV file
    info!("Writing {}", args.outfile.display());
    let mut writer = csv::WriterBuilder::new()
        .has_headers(false)
        .delimiter(delimiter)
        .from_path(args.outfile.clone())?;

    // write the column names
    writer.write_field(first_col)?;
    writer.write_record(header)?;

    // write the rows to the CSV file
    let row_iter = counts_mtx.outer_iterator();
    for (row, row_name) in row_iter
        .zip(row_names.iter())
        .progress_with(create_progress_bar(row_names.len()))
    {
        let row_val_iter = RowValIter::new(&row);
        writer.serialize(Row {
            name: row_name,
            values: row_val_iter,
        })?;
    }

    info!("Done writing {}", args.outfile.display());

    Ok(())
}

/// extract arrays from anndata and determine datatype of counts matrix
fn make_csv(adata: &PyAny, args: LibOpts) -> anyhow::Result<()> {
    let indptr = adata.getattr("X")?.getattr("indptr")?.extract()?;
    let indices = adata.getattr("X")?.getattr("indices")?.extract()?;
    let obs_vec = adata.getattr("obs_names")?.extract()?;
    let var_vec = adata.getattr("var_names")?.extract()?;
    let data = adata.getattr("X")?.getattr("data")?;

    let data_dtype: String = adata
        .getattr("X")?
        .getattr("data")?
        .getattr("dtype")?
        .getattr("name")?
        .to_string();
    match data_dtype.as_str() {
        "int8" => arrays_to_csv(args, data.extract::<Vec<i8>>()?, indptr, indices, obs_vec, var_vec),
        "int16" => arrays_to_csv(args, data.extract::<Vec<i16>>()?, indptr, indices, obs_vec, var_vec),
        "int32" => arrays_to_csv(args, data.extract::<Vec<i32>>()?, indptr, indices, obs_vec, var_vec),
        "int64" => arrays_to_csv(args, data.extract::<Vec<i64>>()?, indptr, indices, obs_vec, var_vec),
        "uint8" => arrays_to_csv(args, data.extract::<Vec<u8>>()?, indptr, indices, obs_vec, var_vec),
        "uint16" => arrays_to_csv(args, data.extract::<Vec<u16>>()?, indptr, indices, obs_vec, var_vec),
        "uint32" => arrays_to_csv(args, data.extract::<Vec<u32>>()?, indptr, indices, obs_vec, var_vec),
        "uint64" => arrays_to_csv(args, data.extract::<Vec<u64>>()?, indptr, indices, obs_vec, var_vec),
        "float32" => arrays_to_csv(args, data.extract::<Vec<f32>>()?, indptr, indices, obs_vec, var_vec),
        "float64" => arrays_to_csv(args, data.extract::<Vec<f64>>()?, indptr, indices, obs_vec, var_vec),
        _ => Err(anyhow!("Invalid data type\nPossible data types: i8, i16, i32, i64, u8, u16, u32, u64, f32, f64")),
    }
}
