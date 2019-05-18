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

    def delete_by_indeces(self, indeces_to_delete):
        self.repetitions = np.delete(self.repetitions, indeces_to_delete, axis=0)

    def get_averages_over_repetitions(self):
        return np.average(self.repetitions, axis=1)

    def calc_std_for_hist(self, normalized=True):
        standard_deviations = np.array([])
        for gene_repetitions in self.repetitions:
            if normalized:
                standard_deviations = np.append(standard_deviations,
                                                np.std(gene_repetitions) / np.average(gene_repetitions))
            else:
                standard_deviations = np.append(standard_deviations,
                                                np.std(gene_repetitions))
        return standard_deviations

