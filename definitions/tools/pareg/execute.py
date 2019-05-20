import sys
import pandas as pd

import sh

from utils import Executor


class MyExecutor(Executor):
    def setup(self):
        self.df_inp.to_csv('genes.csv')
        self.df_terms.to_csv('terms.csv')

    def execute(self):
        sh.Rscript(
            'script.R',
            _out=sys.stdout, _err=sys.stderr)

        self.df_result = pd.read_csv('result.csv')


if __name__ == '__main__':
    ex = MyExecutor(sys.argv[1], sys.argv[2])
    ex.run()
