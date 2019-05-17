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

    def get_p_values(self):
        calculator_pValue = PvalueCalculator(self.conditions[0], self.conditions[1])
        return calculator_pValue.calc_p_values()

    def calculate_fold_changes(self):
        calculator_fold_change = FoldChangeCalculator(self.conditions[0], self.conditions[1])
        return calculator_fold_change.calculate_fold_change()

    def delete_by_threshold(self, threshold):
        union_indeces = np.array([])

        for i in range(self.conditions.shape[0]):
            indeces = self.conditions[i].get_gene_indeces_under_threshold(threshold)
            union_indeces = np.union1d(union_indeces, indeces)

        for condition in self.conditions:
            condition.delete(union_indeces)

        return union_indeces

    def get_std_of_pair(self, normalized=True):
        std_left_condition = self.conditions[0].calc_std_for_hist(normalized)
        std_right_condition = self.conditions[1].calc_std_for_hist(normalized)
        return std_left_condition, std_right_condition

    def delete_by_std(self):
        union_indeces_condition = np.array([])
        union_indeces_conditions = np.array([])

        for i in range(self.conditions.shape[0]):
            for j, gene_repetitions in enumerate(self.conditions[i].repetitions):
                std = np.std(gene_repetitions) / np.average(gene_repetitions)
                if (std > 0.5):
                    union_indeces_condition = np.union1d(union_indeces_condition, j)
            union_indeces_conditions = np.union1d(union_indeces_conditions, union_indeces_condition)

        for condition in self.conditions:
            condition.delete(union_indeces_conditions)

        return union_indeces_conditions

    def get_gene_most_appeared(self):
        averages_condition_1 = np.average(self.conditions[0].repetitions, axis=1)
        condition_1_indeces = np.arange(len(averages_condition_1))[np.argsort(averages_condition_1)][
                              int(2 * len(averages_condition_1) / 3):]

        averages_condition_2 = np.average(self.conditions[1].repetitions, axis=1)
        condition_2_indeces = np.arange(len(averages_condition_2))[np.argsort(averages_condition_2)][
                              int(2 * len(averages_condition_2) / 3):]

        intersect_indeces = np.intersect1d(condition_1_indeces, condition_2_indeces)
        a1 = self.conditions[0].repetitions[:, 0][intersect_indeces]
        a2 = self.conditions[0].repetitions[:, 1][intersect_indeces]
        a3 = self.conditions[0].repetitions[:, 2][intersect_indeces]
        a4 = self.conditions[0].repetitions[:, 3][intersect_indeces]
        self.conditions[0].repetitions = np.array([a1, a2, a3, a4]).T

        b1 = self.conditions[1].repetitions[:, 0][intersect_indeces]
        b2 = self.conditions[1].repetitions[:, 1][intersect_indeces]
        b3 = self.conditions[1].repetitions[:, 2][intersect_indeces]
        b4 = self.conditions[1].repetitions[:, 3][intersect_indeces]
        self.conditions[1].repetitions = np.array([b1, b2, b3, b4]).T
