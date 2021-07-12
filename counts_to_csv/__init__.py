from .counts_to_csv import *

def counts_to_csv(adata, delimiter="comma", column_orient="var-names", outfile="out.csv"):
    make_csv(adata, delimiter, column_orient, outfile)