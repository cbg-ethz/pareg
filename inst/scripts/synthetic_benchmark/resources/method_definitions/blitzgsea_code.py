import sys

import pandas as pd

import blitzgsea


def main(fname_genes, fname_terms, fname_out):
    df_genes = pd.read_csv(fname_genes, header=None, names=["gene", "p_value"])
    df_terms = pd.read_csv(fname_terms, header=None, names=["term", "gene"])

    term_dict = df_terms.groupby("term")["gene"].apply(list).to_dict()

    df_res = blitzgsea.gsea(df_genes, term_dict)
    df_res.to_csv(fname_out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
