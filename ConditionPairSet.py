import numpy as np

from FoldChangeCalculator import FoldChangeCalculator
from PvalueCalculator import PvalueCalculator

moc_minus_index = 0
moc_plus_index = 1
trn_minus_index = 2
trn_plus_index = 3


# distance_func = lambda x, y: x-y
def distance_func(x, y):
    return x / y


class ConditionPairSet:
    def __init__(self, conditions):
        self.conditions = conditions

    def get_p_value_pair(self):
        calculator_pValue = PvalueCalculator(self.conditions[0], self.conditions[1])
        return calculator_pValue.calc_p_value()

    def calculate_fold_changes(self):
        calculator_fold_change = FoldChangeCalculator(self.conditions[0], self.conditions[1])
        return calculator_fold_change.calculate_fold_change()

    def delete_by_threshold(self, threshold):
        union_indeces = np.array([])

        for i in range(self.conditions.shape[0]):
            indeces = self.conditions[i].get_gene_indeces_under_threshold(threshold)
            union_indeces = np.union1d(union_indeces, indeces)

        for condition in self.conditions:
            condition.delete_and_log(union_indeces)

        return union_indeces
