import numpy as np


class Condition:
    """
    Represents all repeated libraries of one condition
    """

    def __init__(self, condition_index, condition_repetition_indeces, libraries):
        """
        condition_index : the index of the condition {0, ... num_conditions-1 }
        condition_repetition_indeces : array of num_repetitions indeces that index
        the repetitions in the libraries for this condition
        """
        self.condition_index = condition_index
        self.repetition_indeces = condition_repetition_indeces

        # num_genes x num_repetitions
        self.repetitions = libraries[:, condition_repetition_indeces]

    def get_gene_indeces_under_threshold(self, threshold):
        """
        returns the indeces of the genes (rows) for which any of the frequency values in
        any of the repetitions (columns) is lower or equal to threshold
        """
        return np.where(self.repetitions <= threshold)[0]

    def delete_and_log(self, indeces_to_delete):
        self.repetitions = np.log(np.delete(self.repetitions, indeces_to_delete, axis=0))
        x = 2

    def get_averages_over_repetitions(self):
        return np.average(self.repetitions, axis=1)
