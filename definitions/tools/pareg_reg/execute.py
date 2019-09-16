import os
import sys
import pandas as pd

import sh

from utils import Executor


class MyExecutor(Executor):
    def setup(self):
        (self.df_inp.rename(columns={'p_value': 'pvalue'})
                    .to_csv('genes.csv', index=False))
        (self.df_terms.rename(columns={'term': 'name'})
                      .to_csv('terms.csv', index=False))

        assert self.df_network is not None, 'No term network loaded'
        self.df_network.to_csv('term_network.csv')

    def execute(self):
        root = os.path.dirname(os.path.realpath(__file__))
        sh.Rscript(
            os.path.join(root, 'script.R'),
            _out=sys.stdout, _err=sys.stderr)

        self.df_result = pd.read_csv('enrichment_result.csv').rename(columns={
            'name': 'term',
            'enrichment': 'p_value'
        })


if __name__ == '__main__':
    ex = MyExecutor(sys.argv[1], sys.argv[2])
    ex.set_term_network(sys.argv[3])
    ex.run()
