use anyhow;
use clap::Clap;
use csv;
use hdf5::types::VarLenUnicode;
use serde::Serialize;
use sprs::CsMatBase;
use std::path::PathBuf;

#[derive(Clap)]
#[clap(about = "Write counts matrix of H5AD to CSV file")]
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
        about = "orient CSV with var-names as columns or obs-names as columns"
    )]
    column_orient: Orient,
    #[clap(short, long, about = "delimiter", default_value = ",")]
    delimiter: char,
    #[clap(
        short,
        long,
        about = "path to output CSV file",
        default_value = "out.csv"
    )]
    outfile: PathBuf,
}

#[derive(Clap)]
enum Orient {
    VarNames,
    ObsNames,
}

#[derive(Serialize)]
struct Row<'a, T, U> where T: Iterator<Item = f32>, U: Iterator<Item = usize> {
    name: &'a str,
    values: RowValIter<'a, T, U>
}

#[derive(Clone)]
struct RowValIter<'a, T, U> where T: Iterator<Item = f32>, U: Iterator<Item = usize> {
    data: T,
    indices: U,
    index: usize,
}

impl<'a, T, U> RowValIter<'a, T, U> {
    fn new(mut values: sprs::vec::CsVecBase<&[usize], &[f32]>) -> Self {
        let data_iter = values.data().iter();
        let idx_iter = values.indices().iter();
        RowValIter {data: data_iter, indices: idx_iter, index: 0}
    }
}

impl<'a, T, U> Iterator for RowValIter<'a, T, U> where T: Iterator<Item = f32>, U: Iterator<Item = usize> {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {

    }

    fn next(&mut self) -> Option<Self::Item> {
        match self.next {
            Some((ind, val)) => {
                self.index += 1;
                if ind == self.index - 1 {
                    self.next = self.values.next();
                    Some(*val)
                } else {
                    Some(0.0)
                }
            },
            None => None,
        }
    }
}

impl Serialize for RowValIter<'_> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error> where S: serde::Serializer, {
        serializer.collect_seq(*self)
        // serializer.collect_seq(self.0.borrow_mut().by_ref())
    }
}

fn main() -> anyhow::Result<()> {
    let args = Opts::parse();

    let delimiter = args.delimiter as u8; // can't convert :(

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

    let mut counts_mtx =
        CsMatBase::try_new((obs_vec.len(), genes_vec.len()), indptr, indices, data)?;

    let (header, first_col, column_names) = match args.column_orient {
        Orient::VarNames => {
            counts_mtx.transpose_mut();
            (obs_vec, "gene", genes_vec)
        }
        Orient::ObsNames => (genes_vec, "cell", obs_vec),
    };

    // let mut row_iter = counts_mtx.outer_iterator();
    // let first_row: sprs::CsVecBase<&[usize], &[f32]> = row_iter.next().unwrap();
    // let mut first_row_iter: sprs::vec::VectorIterator<f32, usize> = first_row.iter();
    // println!("{:?}", first_row_iter.next());

    let mut writer = csv::WriterBuilder::new()
        .has_headers(false)
        .delimiter(delimiter)
        .from_path(args.outfile)?;
    let col_num = header.len();
    let row_iter = counts_mtx.outer_iterator();
    let mut i = 0;
    writer.write_field(first_col)?;
    writer.write_record(header)?;
    // use zip!!
    for row in row_iter {
        let row_val_iter = RowValIter::new(row);
        writer.serialize(Row {
            name: column_names[i].as_str(),
            values: row_val_iter,
        })?;
        // writer.write_field(column_names[i].as_str())?;
        // writer.write_record(as_bytes(&row_vec))?;
        i += 1;
        break;
    }

    Ok(())
}
