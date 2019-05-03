import numpy as np

moc_minus_index = 0
moc_plus_index = 1
trn_minus_index = 2
trn_plus_index = 3


# distance_func = lambda x, y: x-y
def distance_func(x, y):
    return x / y


class TestConditionSet:
    def __init__(self, conditions):
        self.conditions = conditions

    def delete_by_threshold(self, threshold):

        union_indeces = np.array([])
        for i in range(self.conditions.shape[0]):
            indeces = self.conditions[i].get_gene_indeces_under_threshold(threshold)
            union_indeces = np.union1d(union_indeces, indeces)

        for condition in self.conditions:
            condition.delete_and_log(union_indeces)

        return union_indeces

    def calculate_ratio_of_average(self):

        averages = [
            condition.get_averages_over_repetitions()
            for condition in self.conditions
        ]

        ratio_moc = distance_func(averages[moc_minus_index], averages[moc_plus_index])
        ratio_trn = distance_func(averages[trn_minus_index], averages[trn_plus_index])

        return ratio_moc, ratio_trn

    def calculate_average_of_ratio(self):

        ratio_moc = distance_func(self.get_condition_repetitions(moc_minus_index),
                                  self.get_condition_repetitions(moc_plus_index))
        ratio_trn = distance_func(self.get_condition_repetitions(trn_minus_index),
                                  self.get_condition_repetitions(trn_plus_index))

        average_moc = np.average(ratio_moc, axis=1)
        average_trn = np.average(ratio_trn, axis=1)

        return average_moc, average_trn

    def get_condition_repetitions(self, index):
        return self.conditions[index].repetitions
