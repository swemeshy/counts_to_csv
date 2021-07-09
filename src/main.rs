use anyhow::anyhow;
use clap::Clap;
use csv;
use hdf5::types::*;
use indicatif::{ProgressBar, ProgressIterator, ProgressStyle};
use log::info;
use serde::Serialize;
use simple_logger::SimpleLogger;
use sprs::CsMatBase;
use std::path::PathBuf;

#[derive(Clap)]
#[clap(about = "Write counts matrix of H5 to CSV file")]
struct Opts {
    #[clap(
        short = 'f',
        long,
        about = "file path of H5 file can be read as an AnnData,
    must have counts matrix formatted in CSR"
    )]
    h5_file: PathBuf,
    #[clap(
        short,
        long,
        arg_enum,
        about = "orient CSV with var-names as column names or obs-names as column names"
    )]
    column_orient: Orient,
    #[clap(
        short,
        long,
        arg_enum,
        about = "delimiter for CSV file",
        default_value = "comma"
    )]
    delimiter: Delimiter,
    #[clap(
        short,
        long,
        about = "file path of output CSV file",
        default_value = "out.csv"
    )]
    outfile: PathBuf,
}

/// argument enum for delimiter
#[derive(Clap)]
enum Delimiter {
    Comma,
    Tab,
    Colon,
    Pipe,
    Semicolon,
}

/// argument enum for column_orient
#[derive(Clap)]
enum Orient {
    VarNames,
    ObsNames,
}

/// trait that describes the type for the data array
trait ArrayDtype: H5Type + Default + Copy + Serialize {}
impl<T> ArrayDtype for T where T: H5Type + Default + Copy + Serialize {}

/// represents a row to be written to the CSV file
#[derive(Serialize)]
struct Row<'a, T: ArrayDtype> {
    name: &'a str,
    values: RowValIter<'a, T>,
}

/// an iterator wrapper for the indices and data of a CsVecBase object that represents one row
/// in the counts matrix
#[derive(Clone)]
struct RowValIter<'a, T> {
    data: std::slice::Iter<'a, T>,
    indices: std::slice::Iter<'a, usize>,
    next: Option<usize>,
    index: usize,
    stop: usize,
}

/// to create a new instance of the RowValIter iterator wrapper
impl<'a, T: ArrayDtype> RowValIter<'a, T> {
    fn new(values: &'a sprs::vec::CsVecBase<&[usize], &[T]>) -> RowValIter<'a, T> {
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
fn create_progress_bar(iter_size: usize) -> ProgressBar {
    let bar = ProgressBar::new(iter_size as u64);
    bar.set_style(ProgressStyle::default_bar().template("[{elapsed}] {wide_bar} {percent}%"));
    bar
}

/// write counts matrix to csv
fn data_to_csv<T: ArrayDtype>(file: hdf5::File, data: Vec<T>, args: Opts) -> anyhow::Result<()> {
    // get delimiter from args.delimiter
    let delimiter = match args.delimiter {
        Delimiter::Comma => b',',
        Delimiter::Tab => b'\t',
        Delimiter::Colon => b':',
        Delimiter::Pipe => b'|',
        Delimiter::Semicolon => b';',
    };

    // get indptr and indices arrays for creating sparse matrix
    let indptr = file.dataset("X/indptr")?.read_1d::<usize>()?.to_vec();
    let indices = file.dataset("X/indices")?.read_1d::<usize>()?.to_vec();

    // get index column name of var and obs dataframes
    let var_index_name = file
        .group("var")?
        .attr("_index")?
        .read_scalar::<VarLenUnicode>()?;
    let obs_index_name = file
        .group("obs")?
        .attr("_index")?
        .read_scalar::<VarLenUnicode>()?;

    // read var and obs index columns
    let var_vec = file
        .dataset(&format!("var/{}", var_index_name.as_str()))?
        .read_1d::<VarLenUnicode>()?
        .to_vec();
    let obs_vec = file
        .dataset(&format!("obs/{}", obs_index_name.as_str()))?
        .read_1d::<VarLenUnicode>()?
        .to_vec();

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

fn main() -> anyhow::Result<()> {
    // initiate logger
    SimpleLogger::new().init()?;

    let args = Opts::parse();

    // read file and determine counts matrix data type
    info!("Reading H5 file");
    let file = hdf5::File::open(&args.h5_file)?;
    let data = file.dataset("X/data")?;
    let data_dtype = data.dtype()?.to_descriptor()?;

    // call data_to_csv based on corresponding matrix data type
    use TypeDescriptor as TD;
    match data_dtype {
        TD::Integer(IntSize::U1) => data_to_csv(file, data.read_1d::<i8>()?.to_vec(), args),
        TD::Integer(IntSize::U2) => data_to_csv(file, data.read_1d::<i16>()?.to_vec(), args),
        TD::Integer(IntSize::U4) => data_to_csv(file, data.read_1d::<i32>()?.to_vec(), args),
        TD::Integer(IntSize::U8) => data_to_csv(file, data.read_1d::<i64>()?.to_vec(), args),
        TD::Unsigned(IntSize::U1) => data_to_csv(file, data.read_1d::<u8>()?.to_vec(), args),
        TD::Unsigned(IntSize::U2) => data_to_csv(file, data.read_1d::<u16>()?.to_vec(), args),
        TD::Unsigned(IntSize::U4) => data_to_csv(file, data.read_1d::<u32>()?.to_vec(), args),
        TD::Unsigned(IntSize::U8) => data_to_csv(file, data.read_1d::<u64>()?.to_vec(), args),
        TD::Float(FloatSize::U4) => data_to_csv(file, data.read_1d::<f32>()?.to_vec(), args),
        TD::Float(FloatSize::U8) => data_to_csv(file, data.read_1d::<f64>()?.to_vec(), args),
        _ => Err(anyhow!("Invalid data type\nPossible data types: i8, i16, i32, i64, u8, u16, u32, u64, f32, f64")),
    }
}
