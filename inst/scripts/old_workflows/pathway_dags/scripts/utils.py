import os
import contextlib
from abc import ABC, abstractmethod

import json
import time

import pandas as pd


@contextlib.contextmanager
def chdir(dir_: str):
    """Temporarily switch directories."""
    curdir = os.getcwd()

    os.makedirs(dir_, exist_ok=True)
    os.chdir(dir_)

    try:
        yield
    finally:
        os.chdir(curdir)


class Executor(ABC):
    def __init__(
        self,
        run_dir, result_fname, meta_fname,

        pathway_dict, reference_set,

        explicit_de_genes=None,  # genelist mode
        df_cell=None, df_expr=None, df_info=None, df_dea=None  # expression mode
    ):
        self.run_dir = run_dir
        self.result_fname = result_fname
        self.meta_fname = meta_fname

        # store input
        self.pathway_dict = pathway_dict
        self.reference_set = reference_set

        self.explicit_de_genes = explicit_de_genes

        self.df_cell = df_cell
        self.df_expr = df_expr
        self.df_info = df_info
        self.df_dea = df_dea

        # misc
        self.meta_data = {}
        self.df_result = None

    @classmethod
    def from_expression_data(
        cls,
        run_dir, fname_result, fname_meta,
        cell_fname, expr_fname, info_fname, dea_fname
    ):
        # read data
        df_cell = (pd.read_csv(cell_fname, dtype={'gene': str})
                     .set_index('node'))
        df_expr = (pd.read_csv(expr_fname, dtype={'node': str})
                     .set_index('node'))
        df_info = pd.read_csv(info_fname)
        df_dea = pd.read_csv(dea_fname, dtype={'node': str})

        # construct pathway mapping
        pathway_dict = {}
        for node in df_cell.index:
            pw = node.split('_')[0]
            pathway_dict.setdefault(pw, set()).add(node)

        # get reference geneset
        reference_set = set(df_expr.index) | set(df_cell.index)

        # contruct object
        return cls(
            run_dir, fname_result, fname_meta,
            pathway_dict, reference_set,
            df_cell=df_cell, df_expr=df_expr, df_info=df_info, df_dea=df_dea
        )

    @classmethod
    def from_gene_lists(
        cls,
        run_dir, fname_result, fname_meta,
        fname_geneset, fname_termdatabase
    ):
        # load genesets
        df_genes = pd.read_csv(fname_geneset)

        explicit_de_genes = set(df_genes.loc[
            df_genes['type'] == 'study',
            'gene'
        ].to_list())
        reference_set = set(df_genes.loc[
            df_genes['type'] == 'all',
            'gene'
        ].to_list())

        # load term database
        df_terms = pd.read_csv(fname_termdatabase)
        pathway_dict = df_terms.groupby('term')['gene'].apply(set).to_dict()

        # contruct object
        return cls(
            run_dir, fname_result, fname_meta,
            pathway_dict, reference_set,
            explicit_de_genes=explicit_de_genes
        )

    @property
    def de_genes(self):
        if self.explicit_de_genes is not None:
            return self.explicit_de_genes

        # custom computation
        pvalue_threshold = .05
        log2foldchange_threshold = 1

        sub = self.df_dea[
            (self.df_dea['pvalue'] <= pvalue_threshold) &
            (abs(self.df_dea['log2FoldChange']) >= log2foldchange_threshold)
        ]
        return set(sub['node'])

    def setup(self):
        """Can be overwritten to implement required setup actions."""
        pass

    def _finalize(self):
        # meta data
        with open(self.meta_fname, 'w') as fd:
            json.dump(self.meta_data, fd)

        # enrichment results
        if self.df_result is None:
            raise RuntimeError(
                'Result was not saved (must define `self.df_result`)')
        if set(self.df_result.columns) != {'term', 'p_value'}:
            raise RuntimeError(
                'Result must have columns: "term", "p_value" '
                f'(has: {self.df_result.columns})')

        self.df_result.to_csv(self.result_fname, index=False)

    def run(self):
        with chdir(self.run_dir):
            # build environment
            self.setup()

            # run and profile tool
            start = time.perf_counter()
            self.execute()
            self.meta_data['exec_time'] = time.perf_counter() - start

            # cleanup
            self._finalize()

    @abstractmethod
    def execute(self):
        pass
