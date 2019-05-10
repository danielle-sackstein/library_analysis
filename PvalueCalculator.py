import numpy as np
import scipy
from scipy.stats import chi2
from mne.stats import bonferroni_correction, fdr_correction


# class PvalueCalculator:
#     def __init__(self, leftCondition, rightCondition):
#         self.leftCondition = leftCondition
#         self.rightCondition = rightCondition

    # def calc_p_values(self):
    #     num_genes = self.leftCondition.repetitions.shape[0]
    #     p_values = np.zeros(num_genes)
    #
    #     for gene_row in range(num_genes):
    #         x = self.leftCondition.repetitions[gene_row, :]
    #         y = self.rightCondition.repetitions[gene_row, :]
    #         p_values[gene_row] = scipy.stats.chisquare(x, y, ddof=1)[1]
    #
    #     p_values_corrected = bonferroni_correction(p_values)[1]
    #
    #     return p_values_corrected

epsilon = 0.0000001

def calc_poisson_product(values):
    vector = np.zeros((values.shape[0]))
    for i in range(values.shape[0]):
        print(i)
        mu = np.average(values[i, :])
        # pmf?
        product = np.product([scipy.stats.poisson.cdf(v, mu) for v in values[i]])
        if product == 0:
            product = epsilon
        vector[i] = product

    return vector

class PvalueCalculator:
    def __init__(self, leftCondition, rightCondition):
        self.leftCondition = leftCondition
        self.rightCondition = rightCondition

    def calc_p_values(self):
        log_likelihood_ratio = self.calc_log_likelihood_ratio()
        value_for_chisquare = -2 * log_likelihood_ratio

        # The Pvalue will  be 1 when chi2.cdf(value_for_chisquare, 1) is zero, that happens when the value_for_chisquare < 0.
        # p_value = chi2.sf(value_for_chisquare, 1)
        p_value = 1- chi2.cdf(value_for_chisquare, 1)
        return p_value

    def calc_log_likelihood_ratio(self):
        left_repetitions = self.leftCondition.repetitions
        right_repetitions = self.rightCondition.repetitions

        left_product = calc_poisson_product(left_repetitions)
        right_product = calc_poisson_product(right_repetitions)
        all_repetitions = np.concatenate((left_repetitions, right_repetitions), axis=1)
        all_product = calc_poisson_product(all_repetitions)

        # todo
        return np.log((left_product * right_product)/ all_product)


if __name__ == '__main__':
    x = np.array([1, 1, 3, 2]).astype(float)
    y = np.array([0, 1, 1, 2]).astype(float)
    a = scipy.stats.chisquare(x, y, ddof=1)

    print(a.pvalue)
