import sys

import numpy as np
import pandas as pd

from utils import Executor


class MyExecutor(Executor):
    def execute(self):
        term_list = list(self.pathway_dict.keys())

        self.df_result = pd.DataFrame({
            'term': term_list,
            'p_value': np.random.uniform(0, 1, size=len(term_list))
        })


if __name__ == '__main__':
    ex = MyExecutor(*sys.argv[1:])
    ex.run()
