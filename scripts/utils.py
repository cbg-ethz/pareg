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
        run_dir,
        cell_fname,
        expr_fname, info_fname, dea_fname,
        result_fname, meta_fname
    ):
        self.cell_fname = cell_fname
        self.expr_fname = expr_fname
        self.info_fname = info_fname
        self.dea_fname = dea_fname

        self.run_dir = run_dir
        self.result_fname = result_fname
        self.meta_fname = meta_fname

        # read data
        self.df_cell = (pd.read_csv(cell_fname, dtype={'gene': str})
                          .set_index('node'))
        self.df_expr = (pd.read_csv(expr_fname, dtype={'node': str})
                          .set_index('node'))
        self.df_info = pd.read_csv(info_fname)
        self.df_dea = pd.read_csv(dea_fname, dtype={'node': str})

        # compute helpful extra information
        self.pathway_dict = {}
        for node in self.df_cell.index:
            pw = node.split('_')[0]
            self.pathway_dict.setdefault(pw, set()).add(node)

        self.reference_set = \
            set(self.df_expr.columns) | set(self.df_cell.index)

        # misc
        self.meta_data = {}
        self.df_result = None

    @property
    def de_genes(self):
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
