import sys
import pandas as pd

from bioinf_common.algorithms import SetEnrichmentComputer

from utils import Executor


class MyExecutor(Executor):
    def execute(self):
        sec = SetEnrichmentComputer(
            self.pathway_dict, self.reference_set,
            alternative_hypothesis='two-sided')
        res = sec.get_terms(self.de_genes)

        self.df_result = (res[['group_name', 'p_value']]
                          .rename(columns={'group_name': 'term'})
                          .copy())


if __name__ == '__main__':
    ex = MyExecutor(*sys.argv[1:])
    ex.run()
