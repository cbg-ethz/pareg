import sys
import pandas as pd

import sh

from utils import Executor


class MyExecutor(Executor):
    def setup(self):
        self.df_inp.to_csv('genes.csv', index=False)
        self.df_terms.to_csv('terms.csv', index=False)

    def execute(self):
        sh.Rscript(
            'script.R',
            _out=sys.stdout, _err=sys.stderr)

        self.df_result = pd.read_csv('enrichment_result.csv').rename(columns={
            'name': 'term',
            'enrichment': 'p_value'
        })


if __name__ == '__main__':
    ex = MyExecutor(sys.argv[1], sys.argv[2])
    ex.run()
