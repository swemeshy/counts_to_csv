use anyhow;
use clap::Clap;
use csv;
use hdf5::types::VarLenUnicode;
use indicatif::ProgressIterator;
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
        about = "path to H5 file can be read as an AnnData,
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
    #[clap(short, long, arg_enum, about = "delimiter", default_value = "comma")]
    delimiter: Delimiter,
    #[clap(
        short,
        long,
        about = "path to output CSV file",
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
    Semicolon
}

/// argument enum for column_orient
#[derive(Clap)]
enum Orient {
    VarNames,
    ObsNames,
}

/// represents a row to be written to the CSV file
#[derive(Serialize)]
struct Row<'a> {
    name: &'a str,
    values: RowValIter<'a>,
}

/// an iterator wrapper for the indices and data of a CsVecBase object that represents one row
/// in the counts matrix
#[derive(Clone)]
struct RowValIter<'a> {
    data: std::slice::Iter<'a, f32>,
    indices: std::slice::Iter<'a, usize>,
    next: Option<usize>,
    index: usize,
    stop: usize,
}

/// to create a new instance of the RowValIter iterator wrapper
impl<'a> RowValIter<'a> {
    fn new(values: &'a sprs::vec::CsVecBase<&[usize], &[f32]>) -> RowValIter<'a> {
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
impl<'a> Iterator for RowValIter<'a> {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next {
            Some(ind) => {
                if self.index == ind {
                    self.next = self.indices.next().cloned();
                    self.index += 1;
                    self.data.next().cloned()
                } else {
                    self.index += 1;
                    Some(0.0)
                }
            }
            None => {
                if self.index < self.stop {
                    self.index += 1;
                    Some(0.0)
                } else {
                    None
                }
            }
        }
    }
}

/// serialization implementation for RowValIter
impl<'a> Serialize for RowValIter<'a> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.collect_seq(self.clone())
    }
}

fn main() -> anyhow::Result<()> {
    // initiate logger
    SimpleLogger::new().init().unwrap();

    let args = Opts::parse();

    // get delimiter from args.delimiter
    let delimiter = match args.delimiter {
        Delimiter::Comma => b',',
        Delimiter::Tab => b'\t',
        Delimiter::Colon => b':',
        Delimiter::Pipe => b'|',
        Delimiter::Semicolon => b';'
    };

    // read counts matrix of H5 file into Vectors for sparse matrix construction
    info!("Reading H5 file");
    let file = hdf5::File::open(args.h5_file)?;

    let genes_vec = file
        .dataset("var/featurekey")?
        .read_1d::<VarLenUnicode>()?
        .to_vec();
    let obs_vec = file
        .dataset("obs/barcodekey")?
        .read_1d::<VarLenUnicode>()?
        .to_vec();

    let data = file.dataset("X/data")?.read_1d::<f32>()?.to_vec();
    let indptr = file.dataset("X/indptr")?.read_1d::<usize>()?.to_vec();
    let indices = file.dataset("X/indices")?.read_1d::<usize>()?.to_vec();

    // construct sparse matrix
    let mut counts_mtx =
        CsMatBase::try_new((obs_vec.len(), genes_vec.len()), indptr, indices, data)?;

    // transpose matrix if needed based on column orientation specified
    let (header, first_col, row_names) = match args.column_orient {
        Orient::ObsNames => {
            counts_mtx.transpose_mut();
            counts_mtx = counts_mtx.to_csr();
            (obs_vec, "gene", genes_vec)
        }
        Orient::VarNames => (genes_vec, "cell", obs_vec),
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
    for (row, row_name) in row_iter.zip(row_names.iter()).progress() {
        let row_val_iter = RowValIter::new(&row);
        writer.serialize(Row {
            name: row_name,
            values: row_val_iter,
        })?;
    }

    info!("Done writing {}", args.outfile.display());

    Ok(())
}
