import numpy as np
import scipy
from scipy.stats import chi2


class PvalueCalculator:
    def __init__(self, leftCondition, rightCondition):
        self.leftCondition = leftCondition
        self.rightCondition = rightCondition

    def calc_p_values(self):
        num_genes = self.leftCondition.repetitions.shape[0]
        p_values = np.zeros(num_genes)

        for gene_row in range(num_genes):
            x = self.leftCondition.repetitions[gene_row, :]
            y = self.rightCondition.repetitions[gene_row, :]
            p_values[gene_row] = scipy.stats.chisquare(x, y, ddof=1)[1]

        return p_values


if __name__ == '__main__':
    x = np.array([1, 1, 3, 2]).astype(float)
    y = np.array([0, 1, 1, 2]).astype(float)
    a = scipy.stats.chisquare(x, y, ddof=1)

    print(a.pvalue)
