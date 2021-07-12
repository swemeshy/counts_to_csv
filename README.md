# Counts to CSV

Writes the counts matrix in an AnnData file to a CSV file. The file must be an H5 or H5AD file, and it must have the CSR formatted counts matrix in `/X/data`.

## Install as a Python package

### Installation

You must have Rust installed. Instructions to install Rust are here: [https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install)

You must have `setuptools-rust`. Install this package using `pip install setuptools-rust`.

Clone the repository by running `git clone https://github.com/swemeshy/counts_to_csv.git` in your terminal.

Then in the repository, run `pip install -e .` to install the Python package.

### Arguments

This package has only one function, `counts_to_csv`, which has the following parameters:

`adata`: the AnnData object.

`delimiter`: delimiter for the CSV file. Possible options: `comma`, `tab`, `colon`, `pipe`, `semicolon`. Default is `comma`.

`column_orient`: orient the CSV file with var-names as column names or obs-names as column names. Possible options: `var-names`, `obs-names`. Default is `var-names`.

`outfile`: file path of the output CSV file. Default is `out.csv`.

### Example Usage

```
>>> import scanpy as sc
>>> import counts_to_csv as ctc
>>> adata = sc.read('anndata.h5ad')
>>> ctc.counts_to_csv(adata, "comma", "obs-names", "anndata.csv")
```
## Compile and use Rust binary
### Installation

You must have Rust installed. Instructions to install Rust are here: [https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install)

Clone the repository by running `git clone https://github.com/swemeshy/counts_to_csv.git` in your terminal.

Then in the repository, run `cargo build --release`. The compiled binary will be located here: `target/release/counts_to_csv`

For ease of running the binary, you can either
* Add the path to the binary to your `PATH`
* Move the binary to a folder on your `PATH`
* Create a bash alias in your startup file. For example, for UNIX users, you can add this to `.bashrc`:

`alias counts_to_csv="/path/to/repo/target/release/counts_to_csv"`

Make sure to reload your startup file!
### Arguments

`-f, --h5-file <h5-file>`

file path of H5 file that must be readable as an AnnData, and must have the counts matrix in CSR format

`-c, --column-orient <column-orient>`

orient the CSV file with var-names as column names or obs-names as column names

`-d, --delimiter <delimiter>`

delimiter for the CSV file

`-o, --outfile <outfile>`

file path of the output CSV file

