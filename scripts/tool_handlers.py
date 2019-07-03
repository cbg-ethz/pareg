import collections
from abc import ABC, abstractmethod

import numpy as np


class EnrichmentScoreHandler(ABC):
    @staticmethod
    @abstractmethod
    def transform_pvalues(self, pvalues):
        pass

    @staticmethod
    @abstractmethod
    def threshold_pvalues(self, pvalues):
        pass


class PValueHandler(EnrichmentScoreHandler):
    def transform_pvalues(self, pvalues):
        tmp = np.asarray(pvalues)
        tmp[tmp == 0] = np.nan
        tmp_trans = -np.log10(tmp)

        # set "most significant" values
        #tmp_trans[np.isnan(tmp_trans)] = tmp_trans[~np.isnan(tmp_trans)].max()

        return tmp_trans

    def threshold_pvalues(self, pvalues):
        return np.asarray(pvalues) <= .05


class RegressionCoefficientHandler(EnrichmentScoreHandler):
    def transform_pvalues(self, pvalues):
        return abs(np.asarray(pvalues))

    def threshold_pvalues(self, pvalues):
        return abs(np.asarray(pvalues)) >= 1


TRANSFORMER_DICT = collections.defaultdict(PValueHandler)
TRANSFORMER_DICT['pareg'] = RegressionCoefficientHandler()
TRANSFORMER_DICT['pareg_lasso'] = RegressionCoefficientHandler()


if __name__ == '__main__':
    print(TRANSFORMER_DICT['dummy'].transform_pvalues([0.05, 0.1, 1, 0]))
    print(TRANSFORMER_DICT['dummy'].threshold_pvalues([0.05, 0.1, 1, 0]))

    print(TRANSFORMER_DICT['pareg'].transform_pvalues([0.05, 0.1, 1, 0]))
    print(TRANSFORMER_DICT['pareg'].threshold_pvalues([0.05, 0.1, 1, 0]))
