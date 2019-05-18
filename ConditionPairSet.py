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

        calculator_p_value = PvalueCalculator(self.conditions[0], self.conditions[1])
        self.p_values = calculator_p_value.calc_p_values()

        calculator_fold_change = FoldChangeCalculator(self.conditions[0], self.conditions[1])
        self.fold_changes = calculator_fold_change.calculate_fold_change()

    def get_p_values(self):
        return self.p_values

    def get_fold_changes(self):
        return self.fold_changes

    def delete_by_single_condition_threshold(self, threshold):
        union_indeces = np.array([])

        for i in range(self.conditions.shape[0]):
            indeces = self.conditions[i].get_gene_indeces_under_threshold(threshold)
            union_indeces = np.union1d(union_indeces, indeces)

        self.delete_by_indeces(union_indeces)

        return union_indeces

    def get_gene_indeces_by_p_value_threshold(self, max_p_value_threshold):
        p_values = self.get_p_values()
        return np.where(p_values > max_p_value_threshold)[0]

    def get_gene_indeces_by_fold_change_threshold(self, fold_change_threshold):
        fold_changes = self.get_fold_changes()
        if fold_change_threshold > 0:
            return np.where(fold_changes < fold_change_threshold)[0]
        else:
            return np.where(fold_changes > -fold_change_threshold)[0]

    def get_gene_indeces_by_pair_threshold(self, max_p_value_threshold, fold_change_threshold):
        unreliable_indeces = self.get_gene_indeces_by_p_value_threshold(max_p_value_threshold)
        insignificant_indeces = self.get_gene_indeces_by_fold_change_threshold(fold_change_threshold)

        return np.union1d(unreliable_indeces, insignificant_indeces)

    def get_std_of_pair(self, normalized=True):
        std_left_condition = self.conditions[0].calc_std_for_hist(normalized)
        std_right_condition = self.conditions[1].calc_std_for_hist(normalized)
        return std_left_condition, std_right_condition

    def delete_by_indeces(self, indeces_to_delete):
        for condition in self.conditions:
            condition.delete_by_indeces(indeces_to_delete)
        np.delete(self.p_values, indeces_to_delete, axis=0)
        np.delete(self.fold_changes, indeces_to_delete, axis=0)

    def delete_by_std(self):
        union_indeces_condition = np.array([])
        union_indeces_conditions = np.array([])

        for i in range(self.conditions.shape[0]):
            for j, gene_repetitions in enumerate(self.conditions[i].repetitions):
                std = np.std(gene_repetitions) / np.average(gene_repetitions)
                if (std > 0.5):
                    union_indeces_condition = np.union1d(union_indeces_condition, j)
            union_indeces_conditions = np.union1d(union_indeces_conditions, union_indeces_condition)

        self.delete_by_indeces(union_indeces_conditions)

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
