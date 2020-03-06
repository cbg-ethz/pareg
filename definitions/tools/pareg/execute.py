import os
import sys
import pandas as pd

import sh

from utils import Executor


class MyExecutor(Executor):
    def setup(self):
        (self.df_dea.rename(columns={'p_value': 'pvalue', 'node': 'gene'})
                    .loc[:, ['gene', 'pvalue']]
                    .to_csv('genes.csv', index=False))

        (pd.DataFrame.from_dict(self.pathway_dict, orient='index')
                     .reset_index()
                     .melt(id_vars=['index']).drop('variable', axis=1)
                     .rename(columns={'index': 'name', 'value': 'gene'})
                     .sort_values('name')
                     .to_csv('terms.csv', index=False))

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
    ex = MyExecutor(*sys.argv[1:])
    ex.run()
