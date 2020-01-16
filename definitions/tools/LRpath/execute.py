"""
Python reimplementation of the original LRPath R code
"""

import sys

import numpy as np
import pandas as pd

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import PerfectSeparationError

from utils import Executor


class MyExecutor(Executor):
    def setup(self):
        self.gene_dict = self.df_dea.set_index('node').to_dict()['pvalue']

    def execute(self):
        results = []
        for term, gset in self.pathway_dict.items():
            # print(self.genes, term, gset)

            gene_list = []
            pvalue_list = []
            for g, p in self.gene_dict.items():
                gene_list.append(g)
                pvalue_list.append(p)

            df = pd.DataFrame({
                'y': np.array([1 if g in gset else 0
                               for g in gene_list]),
                'X': -np.log10(np.array(pvalue_list))
            })
            # print(df)

            logit = sm.genmod.families.links.logit()
            fam = sm.families.Binomial(link=logit)
            model = smf.glm('y ~ X', data=df, family=fam)

            fit = model.fit()

            # print(fit.summary())
            results.append({
                'term': term,
                'log_odds': fit.params['X'],
                'p_value': fit.pvalues['X']
            })

        df_res = pd.DataFrame(results)
        self.df_result = df_res[['term', 'p_value']].copy()


if __name__ == '__main__':
    ex = MyExecutor(*sys.argv[1:])
    ex.run()
