import os
import sys
import pandas as pd

import sh

from utils import Executor


class MyExecutor(Executor):
    def execute(self):
        root = os.path.dirname(os.path.realpath(__file__))
        sh.Rscript(
            os.path.join(root, 'script.R'),
            self.expr_fname, self.info_fname,
            _out=sys.stdout, _err=sys.stderr)

        self.df_result = pd.read_csv('result.csv').rename(columns={
            'pathway': 'term',
            'PValue': 'p_value'
        }).drop(columns=['NGenes', 'Direction', 'FDR'])


if __name__ == '__main__':
    ex = MyExecutor(*sys.argv[1:])
    ex.run()
