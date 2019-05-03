import numpy as np
import scipy
from scipy.stats import chi2

NUM_CONDITIONS = 7
NUM_REPETITIONS = 4
EPSILON = 0.000000001

CONDITIONS_DICT = {"WellFed": 0, "moc": 1, "trained": 2, "moc-": 3, "moc+": 4, "trained-": 5,
                   "trained+": 6}


def gen_condition_pairs():
    for i in range(NUM_CONDITIONS):
        for j in range(i + 1, NUM_CONDITIONS):
            yield (i, j)


def get_repetitions(condition_index):
    for repetition_index in range(NUM_REPETITIONS):
        yield condition_index + repetition_index * NUM_CONDITIONS


def calc_poisson_product(values):
    mu = np.average(values)
    return np.product([scipy.stats.poisson.cdf(v, mu) for v in values])


# def calc_poisson_product(values):
#     mu = np.average(values)
#     product = 1
#     for v in values:
#         power = np.power(np.e, -mu) * np.power(mu, v)
#         factorial = np.math.factorial(np.round(v))
#         if power == 0:
#             power = EPSILON
#         if factorial == 0:
#             factorial = EPSILON
#
#         poisson = np.log(power) - math.log(factorial)
#         poisson = np.power(math.e, poisson)
#         # poisson = (np.power(np.e, -mu) * np.power(mu, v)) / np.math.factorial(np.round(v))
#         # poisson = scipy.stats.poisson.pmf(v, mu)
#         product = poisson * product
#         if poisson == float("inf") or product == float("inf"):
#             a = 5
#     return product


def calc_log_likelihood_ratio(left_frequencies, right_frequencies):
    left_product = calc_poisson_product(left_frequencies)
    right_product = calc_poisson_product(right_frequencies)
    all_repetitions = np.concatenate((left_frequencies, right_frequencies))
    all_product = calc_poisson_product(all_repetitions)

    return np.log(all_product / (left_product * right_product))


def calc_p_value(log_likelihood_ratio):
    value_for_chisquare = -2 * log_likelihood_ratio

    # The Pvalue will  be 1 when chi2.cdf(value_for_chisquare, 1) is zero, that happens when the value_for_chisquare < 0.
    p_value = 1 - chi2.cdf(value_for_chisquare, 1)
    return p_value


def get_repeated_frequencies(gene_row, condition_index):
    frequencies = [gene_row[i] for i in get_repetitions(condition_index)]
    return np.array(frequencies)


def gen_likelihood_ratios_for_gene(gene_row):
    for condition_pair in gen_condition_pairs():
        left_frequencies = get_repeated_frequencies(gene_row, condition_pair[0])
        right_frequencies = get_repeated_frequencies(gene_row, condition_pair[1])
        log_likelihood_ratio = calc_log_likelihood_ratio(left_frequencies, right_frequencies)
        yield condition_pair, log_likelihood_ratio, calc_p_value(log_likelihood_ratio)


def gen_likelihood_ratios_for_all_genes(input_array):
    for gene_row in input_array:
        yield np.array(list((gen_likelihood_ratios_for_gene(gene_row))))


def run(dir, input_array, genes_id):
    # TODO: do this line
    all_libraries_frequencies = input_array[:, 1:].astype(float)
    # todo delete this line
    # all_libraries_frequencies = input_array.astype(float)

    # outfile_logLikelihood_ratios = open(dir + "\\logLikelihood.ratios.txt", "w")
    outfile_logLikelihood_ratios = open(dir + "\\logLiklihood.ratios.txt", "w")
    outfile_logLikelihood_ratios.write("\tConditions tuples -> ")

    # outfile_Pvalues = open(dir + "\\Pvalues.txt", "w")
    outfile_Pvalues = open(dir + "\\Pvalues.txt", "w")
    outfile_Pvalues.write("\tConditions tuples -> ")

    for condition in gen_condition_pairs():
        outfile_logLikelihood_ratios.write(str(condition) + "\t")
        outfile_Pvalues.write(str(condition) + "\t")

    outfile_logLikelihood_ratios.write("\n")
    outfile_Pvalues.write("\n")

    outfile_logLikelihood_ratios.write("Genes Id:\n")
    outfile_Pvalues.write("Genes Id:\n")

    for j, gene_likelihood_ratio in enumerate(gen_likelihood_ratios_for_all_genes(all_libraries_frequencies)):
        outfile_logLikelihood_ratios.write(str(genes_id[j]) + "\t")
        outfile_Pvalues.write(str(genes_id[j]) + "\t")
        for library_frequency in gene_likelihood_ratio:
            outfile_logLikelihood_ratios.write(str(library_frequency[1]) + "\t")
            outfile_Pvalues.write(str(library_frequency[2]) + "\t")
            if library_frequency[2] < 0.07:
                print("The Pvalue of the gene id:{}, is: {}".format(genes_id[j], library_frequency[2]))
        outfile_logLikelihood_ratios.write("\n")
        outfile_Pvalues.write("\n")


if __name__ == "__main__":

    datadir = "./data/"
    genes_names_file = datadir + "genes.names.npy"
    frequencies_file = datadir + "libraries.frequencies.npy"

    genes_names = np.load(genes_names_file)[:, 0].astype(str)
    frequencies = np.load(frequencies_file).astype(float)

    run(dir, frequencies, genes_names)
