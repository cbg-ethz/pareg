from abc import ABC, abstractmethod

import json
import time

import pandas as pd


def load_input(fname):
    return pd.read_csv(fname)


def load_terms(fname):
    return pd.read_csv(fname)


class Executor(ABC):
    def __init__(self, input_file, term_file):
        self.df_inp = load_input(input_file)
        self.df_terms = load_terms(term_file)

        self.reference_set = \
            set(self.df_inp['gene'].tolist()) | \
            set(self.df_terms['gene'].tolist())

        self.meta_data = {}
        self.df_result = None

    def setup(self):
        pass

    def finalize(self):
        # meta data
        with open('meta.json', 'w') as fd:
            json.dump(self.meta_data, fd)

        # enrichment results
        if self.df_result is None:
            raise RuntimeError(
                'Result was not saved (must define `self.df_result`)')
        if set(self.df_result.columns) != {'term', 'p_value'}:
            raise RuntimeError(
                'Result must have columns: "term", "p_value" '
                f'(has: {self.df_result.columns})')

        self.df_result.to_csv('result.csv', index=False)

    def run(self):
        # build environment
        self.setup()

        # run and profile tool
        start = time.perf_counter()
        self.execute()
        self.meta_data['exec_time'] = time.perf_counter() - start

        # cleanup
        self.finalize()

    @abstractmethod
    def execute(self):
        pass
