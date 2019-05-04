import numpy as np
import scipy
from scipy.stats import chi2


def calc_poisson_product(values):
    vector = np.zeros((values.shape[0]))
    for i in range(values.shape[0]):
        print(i)
        mu = np.average(values[i, :])
        vector[i] = np.product([scipy.stats.poisson.cdf(v, mu) for v in values])

    return vector


class PvalueCalculator:
    def __init__(self, leftCondition, rightCondition):
        self.leftCondition = leftCondition
        self.rightCondition = rightCondition

    def calc_p_value(self):
        log_likelihood_ratio = self.calc_log_likelihood_ratio()
        value_for_chisquare = -2 * log_likelihood_ratio

        # The Pvalue will  be 1 when chi2.cdf(value_for_chisquare, 1) is zero, that happens when the value_for_chisquare < 0.
        p_value = 1 - chi2.cdf(value_for_chisquare, 1)
        return p_value

    def calc_log_likelihood_ratio(self):
        left_repetitions = self.leftCondition.repetitions
        right_repetitions = self.rightCondition.repetitions

        left_product = calc_poisson_product(left_repetitions)
        right_product = calc_poisson_product(right_repetitions)
        all_repetitions = np.concatenate((left_repetitions, right_repetitions))
        all_product = calc_poisson_product(all_repetitions)

        return np.log(all_product / (left_product * right_product))
