import sys

import numpy as np
import pandas as pd

from utils import Executor


class MyExecutor(Executor):
    def execute(self):
        term_list = self.df_terms['term'].unique()

        self.df_result = pd.DataFrame({
            'term': term_list,
            'p_value': np.random.uniform(0, 1, size=term_list.size)
        })


if __name__ == '__main__':
    ex = MyExecutor(sys.argv[1], sys.argv[2])
    ex.run()
