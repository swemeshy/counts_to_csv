# Counts to CSV

Writes the counts matrix in an AnnData object to a CSV file. The file must be an H5 or H5AD file, and it must have the counts matrix in `/X`.

Note: Currently only writes matrices of float values.

## Installation

You must have Rust installed. Instructions to install Rust are here: [https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install)

Clone the repository by running `git clone https://github.com/swemeshy/counts_to_csv.git` in your terminal.

Then in the repository, run `cargo build --release`. The compiled binary will be located here: `target/release/counts_to_csv`

For ease of running the binary, you can add the following line to your startup file (e.g. `.bashrc`):

`alias counts_to_csv="path/to/repo/target/release/counts_to_csv"`

Make sure to reload your startup file!

## Arguments

`-c, --column-orient <column-orient>`

orient the CSV file with var-names as column names or obs-names as column names

`-d, --delimiter <delimiter>`

delimiter for the CSV file

`-f, --h5-file <h5-file>`

path to H5 file that must be readable as an AnnData, and must have the counts matrix in CSR format

`-o, --outfile <outfile>`

path and file name of the output CSV file

