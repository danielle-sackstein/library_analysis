import numpy as np
from Condition import Condition
from TestConditionSet import TestConditionSet
from ConditionPairSet import ConditionPairSet

num_repetitions = 4


class ExperimentResults:
    def __init__(self):
        data_dir = "./data/"
        gene_names_file = data_dir + "test_genes.npy"
        libraries_file = data_dir + "test_frequencies.npy"

        # num_genes x total_conditions
        # for each gene total_conditions frequency values, one for each condition-repetition
        libraries = np.load(libraries_file).astype(float)
        self.num_genes = libraries.shape[0]
        self.total_num_conditions = libraries.shape[1]
        self.num_conditions = (int)(self.total_num_conditions / num_repetitions)

        # num_genes x 1
        # for each gene the name of the gene
        self.gene_names = np.load(gene_names_file).astype(str)

        # num_conditions x 1
        self.conditions = self.create_conditions(libraries)

    def _get_condition(self, condition_index: int) -> Condition:
        return self.conditions[condition_index]

    def create_condition_pair_set(self, test_condition_indeces, threshold):
        test_conditions = np.array([self.conditions[i] for i in test_condition_indeces])

        test_condition_set = ConditionPairSet(test_conditions)
        deleted_indeces = test_condition_set.delete_by_threshold(threshold)

        self.gene_names = np.delete(self.gene_names, deleted_indeces, axis=0)

        return test_condition_set

    def create_condition(self, condition_index, libraries):
        condition_indeces = [
            condition_index + i * self.num_conditions for i in range(num_repetitions)
        ]
        return Condition(condition_index, condition_indeces, libraries)

    def create_conditions(self, libraries):
        return np.array([self.create_condition(i, libraries) for i in range(self.num_conditions)])
