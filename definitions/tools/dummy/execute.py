import sys

import pandas as pd

from utils import Executor


class MyExecutor(Executor):
    def execute(self):
        tmp = []
        for term in self.df_terms['name'].unique():
            tmp.append({
                'term': term,
                'pvalue': .5
            })

        self.df_result = pd.DataFrame(tmp)


if __name__ == '__main__':
    ex = MyExecutor(sys.argv[1], sys.argv[2])
    ex.run()
