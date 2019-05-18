import numpy as np
from Condition import Condition
from TestConditionSet import TestConditionSet
from ConditionPairSet import ConditionPairSet

test = ".test"
num_repetitions = 4


class ExperimentResults:
    def __init__(self):
        data_dir = "./data/"
        gene_names_file = data_dir + "gene.names" + ".npy"
        libraries_file = data_dir + "libraries" + ".npy"

        # num_genes x total_conditions
        # for each gene total_conditions frequency values, one for each condition-repetition
        libraries = np.load(libraries_file).astype(float)
        self.num_genes = libraries.shape[0]
        self.total_num_conditions = libraries.shape[1]
        self.num_conditions = (int)(self.total_num_conditions / num_repetitions)

        # num_genes x 1
        # for each gene the name of the gene
        self.gene_transcript_id = np.load(gene_names_file).astype(str)
        # self.gene_names = self.get_gene_names()
        # num_conditions x 1
        self.conditions = self.create_conditions(libraries)
        self.cpm = self.update_cpm(libraries)

    def get_gene_count(self):
        return self.gene_transcript_id.shape[0]

    def get_genes_under_threshold(self, threshold):
        union_indeces = np.array([]).astype(int)

        for i in range(self.conditions.shape[0]):
            indeces = self.conditions[i].get_gene_indeces_under_threshold(threshold)
            union_indeces = np.union1d(union_indeces, indeces)

        return union_indeces

    def delete_genes_by_index(self, indeces):
        for condition in self.conditions:
            condition.delete_by_indeces(indeces)

        self.gene_transcript_id = np.delete(self.gene_transcript_id, indeces, axis=0)
        return self.gene_transcript_id.shape[0]

    def _get_condition(self, condition_index: int) -> Condition:
        return self.conditions[condition_index]

    def create_condition_pair_set_with_delete(self, test_condition_indeces, threshold):
        condition_pair_set = self.create_condition_pair_set(test_condition_indeces)
        deleted_indeces = condition_pair_set.delete_by_single_condition_threshold(threshold)
        self.gene_transcript_id = np.delete(self.gene_transcript_id, deleted_indeces, axis=0)

        return condition_pair_set

    def create_condition_pair_set(self, test_condition_indeces):
        test_conditions = np.array([self.conditions[i] for i in test_condition_indeces])
        return ConditionPairSet(test_conditions)

    def create_test_condition_set(self, test_condition_indeces, threshold):
        test_conditions = np.array([self.conditions[i] for i in test_condition_indeces])

        test_condition_set = TestConditionSet(test_conditions)
        deleted_indeces = test_condition_set.delete_by_threshold(threshold)

        self.gene_transcript_id = np.delete(self.gene_transcript_id, deleted_indeces, axis=0)

        return test_condition_set

    def create_condition(self, condition_index, libraries):
        condition_indeces = [
            condition_index + i * self.num_conditions for i in range(num_repetitions)
        ]
        return Condition(condition_index, condition_indeces, libraries)

    def delete_by_threshold(self, threshold):
        union_indeces = np.array([])
        for i in range(self.conditions.shape[0]):
            indeces = self.conditions[i].get_gene_indeces_under_threshold(threshold)
            union_indeces = np.union1d(union_indeces, indeces)

        for condition in self.conditions:
            condition.delete(union_indeces)

        return union_indeces

    def create_conditions(self, libraries):
        return np.array([self.create_condition(i, libraries) for i in range(self.num_conditions)])

    def get_condition(self, condition_index):
        return self.conditions[condition_index]

    def get_gene_names(self):
        return self.gene_transcript_id

    def calc_averages_of_repetitions(self, condition_index):
        log_repetitions = np.log(self.conditions[condition_index].repetitions)
        sum_1 = log_repetitions[:, 0] + log_repetitions[:, 1]
        sum_2 = log_repetitions[:, 2] + log_repetitions[:, 3]
        return sum_1 / 2, sum_2 / 2

    def update_cpm(self, libraries):
        cpm = np.zeros((libraries.shape[0], libraries.shape[1]))
        for i in range(libraries.shape[1]):
            sum = np.sum(libraries[:, i])
            for j in range(libraries.shape[0]):
                cpm[j, i] = (libraries[j, i] / sum) * 1000000
        return cpm

    def delete_insignificant_genes_by_fold_change(self, cond_1, cond_2, initial_thresholds, reduction_ratio):
        max_p_value_threshold = initial_thresholds[0]
        fold_change_threshold = initial_thresholds[1]

        gene_count = self.get_gene_count()

        required_gene_count = gene_count * reduction_ratio

        condition_pair_set = self.create_condition_pair_set((cond_1, cond_2))

        while gene_count > required_gene_count:
            # find indeces of genes that are not sensitive to starvation at all

            indeces_of_genes_not_sensitive_to_starvation = condition_pair_set.get_gene_indeces_by_pair_threshold(
                max_p_value_threshold,
                fold_change_threshold)

            # delete the indeces
            gene_count = self.delete_genes_by_index(indeces_of_genes_not_sensitive_to_starvation)

            max_p_value_threshold = max_p_value_threshold * 0.9
            fold_change_threshold = fold_change_threshold / 0.9

        return max_p_value_threshold, fold_change_threshold







