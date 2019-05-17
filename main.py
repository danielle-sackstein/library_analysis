import scipy
from scipy.stats import chi2

from ExperimentResults import ExperimentResults
from ShowResults import *
import pandas as pd

from geneDiffCalculator import geneDiffCalculator
import scipy as scipy
import scipy.stats
CONDITIONS = ["Well.Fed", "moc", "trained", "moc-", "moc+", "trained-", "trained+"]

cond_well_fed = 0
cond_moc = 1
cond_trained = 2
cond_moc_minus = 3
cond_moc_plus = 4
cond_trn_minus = 5
cond_trn_plus = 6
threshold = 5

repetitions_num = 4
condition_num = 7

directory = "./data/"


def save_data_to_excel(log_2_fold_changes, pValues):
    data = pd.DataFrame(columns=['fold_change', 'pValues'])
    data['fold_change'] = log_2_fold_changes
    data['pValues'] = pValues
    data.to_csv(directory + "to_volcano.csv", sep=",")


def get_condition_indeces_to_test():
    return [cond_trn_plus, cond_moc_plus, cond_trn_minus, cond_moc_minus]


def run():
    experiment_results = ExperimentResults()

    indeces = get_condition_indeces_to_test()

    test_condition_set = experiment_results.create_test_condition_set(indeces, threshold)

    ratio_moc, ratio_trn = test_condition_set.calculate_fold_changes()

    show_results(experiment_results.gene_transcript_id, ratio_moc, ratio_trn)


def create_volcano():
    experiment_results = ExperimentResults()

    condition_pair_set_indeces = (cond_moc_plus, cond_well_fed)
    # condition_pair_set_indeces = (cond_trn_plus, cond_moc_plus)
    # condition_pair_set_indeces_test = (0,1)

    condition_pair_set = experiment_results.create_condition_pair_set_with_delete(condition_pair_set_indeces, threshold)
    # condition_pair_set.get_gene_most_appeared()

    create_volcano_plot(condition_pair_set, experiment_results)

    # minus_log_pvalues = -np.log(pValues)
    # indeces = np.where(minus_log_pvalues > 1000)
    # minus_log_pvalues[indeces] = 160
    # volcano_plot(log_fold_change, minus_log_pvalues)


def find_most_expressed_gene(experiment_results, log2_fold_change, pValues):
    _logpVales = -np.log10(pValues)
    gene_pval = 2.5
    gene_log_change = 3
    indexes= []
    for i in range(len(log2_fold_change)):
        if gene_pval < _logpVales[i] and gene_log_change < log2_fold_change[i]:
            indexes.append(i)
    for i in indexes:
        gene_name = experiment_results.gene_transcript_id[i]
        create_bars_of_expression_in_conditions(gene_name)
        print(gene_name, "-logPvalue = ", _logpVales[i], "foldChange= ", log2_fold_change[i])
    # index = np.where(experiment_results.gene_transcript_id =="F45E1.6")


def create_volcano_plot(condition_pair_set, experiment_results):
    pValues = condition_pair_set.calculate_p_values()
    fold_changes = condition_pair_set.calculate_fold_changes()
    log2_fold_change = np.log2(fold_changes)
    save_data_to_excel(log2_fold_change, pValues)
    plot_volcano_from_library()

    # indeces_upregulated, indeces_downregulated = find_and_print_de_genes(experiment_results)
    # plot_upregulated_downregulated(condition_pair_set, experiment_results, indeces_upregulated, indeces_downregulated,
    #                                log2_fold_change, pValues)
    find_most_expressed_gene(experiment_results, log2_fold_change, pValues)

    # find_and_print_genes_on_y_axis(condition_pair_set, pValues, fold_changes)


def plot_upregulated_downregulated(condition_pair_set, experiment_results, indeces_upregulated, indeces_downregulated,
                                   log2_fold_change, pValues):
    indeces_upregulated = indeces_upregulated.astype(int)
    indeces_downregulated = indeces_downregulated.astype(int)

    repetition_indeces_left_cond = condition_pair_set.conditions[0].repetition_indeces
    repetition_indeces_right_cond = condition_pair_set.conditions[1].repetition_indeces

    for index in indeces_upregulated:
        gene_name = experiment_results.gene_transcript_id[index]
        # if (
        #         gene_name == "R09F10.4" or gene_name == "C53B7.4" or gene_name == "W01A11.4" or gene_name == "C29E4.7"
        #         or gene_name == "C02C2.3" or gene_name == "Y42G9A.1" or gene_name == "F17A9.4" or
        #         gene_name == "B0496.8" or gene_name == " C02C2.3" or gene_name == "F46H6.1" or gene_name == "F32A5.6"
        #         or gene_name == "C44F1.5" or gene_name == " F17A9.4" or gene_name == " F18F11.1" or gene_name == "C47B2.8"
        #         or gene_name=="W01A11.4" or gene_name == " C29E4.7"):
        repetitions_cond_0 = experiment_results.cpm[index, repetition_indeces_left_cond]
        repetitions_cond_1 = experiment_results.cpm[index, repetition_indeces_right_cond]
        save_gene_plot(repetitions_cond_0, repetitions_cond_1, gene_name)

    for index in indeces_downregulated:
        gene_name = experiment_results.gene_transcript_id[index]
        # if (
        #         gene_name == "R09F10.4" or gene_name == "C53B7.4" or gene_name == "W01A11.4" or gene_name == "C29E4.7"
        #         or gene_name == "C02C2.3" or gene_name == "Y42G9A.1" or gene_name == "F17A9.4" or
        #         gene_name == "B0496.8" or gene_name == " C02C2.3" or gene_name == "F46H6.1" or gene_name == "F32A5.6"
        #         or gene_name == "C44F1.5" or gene_name == " F17A9.4" or gene_name == " F18F11.1" or gene_name == "C47B2.8"
        #         or gene_name == "W01A11.4" or gene_name == " C29E4.7"):
        gene_name = experiment_results.gene_transcript_id[index]
        epetitions_cond_0 = experiment_results.cpm[index, repetition_indeces_left_cond]
        repetitions_cond_1 = experiment_results.cpm[index, repetition_indeces_right_cond]
        save_gene_plot(epetitions_cond_0, repetitions_cond_1, gene_name)


def print_de_genes(down_regulated_genes, up_regulated_genes):
    print("Genes that were up regulated:")
    for gene in up_regulated_genes:
        print(gene)

    print("Genes that were down regulated:")
    for gene in down_regulated_genes:
        print(gene)


def find_and_print_de_genes(experiment_results):
    gene_names = experiment_results.get_gene_names()

    log_fold_change_pvalues = pd.read_csv(directory + "to_volcano.csv", sep=",", dtype=float)
    log_fold_changes = log_fold_change_pvalues['fold_change']
    pvalues = log_fold_change_pvalues['pValues']

    down_regulated_genes = []
    up_regulated_genes = []

    indeces_upregulated = np.array([])
    indeces_downregulated = np.array([])

    for i in range(log_fold_changes.shape[0]):

        if log_fold_changes[i] <= -1.2 and -np.log(pvalues[i]) > 1:
            indeces_upregulated = np.union1d(indeces_upregulated, i)
            gene_name = gene_names[i]
            down_regulated_genes.append(gene_name)

        elif log_fold_changes[i] >= 1.2 and -np.log(pvalues[i]) > 1:
            indeces_downregulated = np.union1d(indeces_downregulated, i)
            gene_name = gene_names[i]
            up_regulated_genes.append(gene_name)

    print_de_genes(down_regulated_genes, up_regulated_genes)
    return indeces_upregulated, indeces_downregulated


def find_and_print_genes_on_y_axis(condition_pair_set, pvalues, fold_changes):
    log_fold_change = np.log2(fold_changes)
    minus_log_pvalue = -np.log10(pvalues)
    count_gents = 0
    for i in range(fold_changes.shape[0]):
        if -0.2 < log_fold_change[i] < 0.2 and minus_log_pvalue[i] > 10:
            x, y = condition_pair_set.conditions[0].repetitions[i, :], \
                   condition_pair_set.conditions[1].repetitions[i, :]
            count_gents += 1
            print(repr(x), repr(y), "pvalue= ", pvalues[i], "foldchange=", fold_changes[i])
    print("num of genes: ", count_gents)


def create_genes_expression_plot():
    experiment_results = ExperimentResults()

    condition_pair_set_indeces = (cond_trn_plus, cond_moc_plus)

    condition_pair_set = experiment_results.create_condition_pair_set_with_delete(condition_pair_set_indeces, threshold)
    log_fold_changes = np.log(condition_pair_set.calculate_fold_changes())

    # get all genes with foldChange 0
    genes_fold_zero_indeces = []
    for i, fold_change in enumerate(log_fold_changes):
        if 0.01 > fold_change > -0.05:
            genes_fold_zero_indeces.append(i)

    gene_diff_calculator = geneDiffCalculator(condition_pair_set.conditions[0], condition_pair_set.conditions[1])
    for i, gene_index in enumerate(genes_fold_zero_indeces[300:600]):
        save_gene_plot(gene_diff_calculator.create_genes_lists(gene_index)[0],
                       gene_diff_calculator.create_genes_lists(gene_index)[1], i)


def get_volcano_std():
    experiment_results = ExperimentResults()
    condition_pair_set_indeces = (cond_trn_plus, cond_moc_plus)
    condition_pair_set = experiment_results.create_condition_pair_set_with_delete(condition_pair_set_indeces, threshold)
    condition_pair_set.delete_by_std()
    condition_pair_set.delete_by_single_condition_threshold(3)
    create_volcano_plot(condition_pair_set, experiment_results)


def create_hist_of_std():
    experiment_results = ExperimentResults()
    experiment_results.delete_by_threshold(5)
    for i in range(7):
        condition = experiment_results.get_condition(i)
        std_condition = condition.calc_std_for_hist()
        plot_hist_of_std(std_condition, CONDITIONS[i])


def create_hist_of_foldchange_std():
    experiment_results = ExperimentResults()
    condition_pair_set_indeces = (cond_trn_plus, cond_moc_plus)
    condition_pair_set = experiment_results.create_condition_pair_set_with_delete(condition_pair_set_indeces, threshold)
    condition_pair_set.delete_by_single_condition_threshold(5)

    left_condition_std, right_condition_std = condition_pair_set.get_std_of_pair(normalized=False)
    left_condition_average = np.average(left_condition_std)
    right_condition_average = np.average(right_condition_std)
    average_diff = abs(right_condition_average - left_condition_average)

    noise_fold_ratio = (left_condition_std + right_condition_std) / average_diff
    sorted = np.sort(noise_fold_ratio)
    indeces_to_cutt = np.where(sorted < 0.5)

    condition_pair_set.conditions[0].repetitions = condition_pair_set.conditions[0].repetitions[indeces_to_cutt]
    condition_pair_set.conditions[1].repetitions = condition_pair_set.conditions[1].repetitions[indeces_to_cutt]

    pValues = condition_pair_set.calculate_p_values()
    fold_changes = condition_pair_set.calculate_fold_changes()

    save_data_to_excel(fold_changes, pValues)
    plot_volcano_from_library()

    # plot_hist_of_foldchange_std(left_condition_std, right_condition_std, conditions[cond_trn_plus],
    #                             conditions[cond_moc_plus])

def create_scattered_of_two_repetitions():
    experiment_results = ExperimentResults()
    max_correlation= 0

    for i in range(28):
        for j in range(i+1,28):

            first_cond = experiment_results.conditions[i].repetitions[:,j]
            second_cond = experiment_results.conditions[i].repetitions[:,j]
            correlation = scipy.stats.pearsonr(first_cond, second_cond)[0]
            if correlation == 1:
                correlation = 0
            if correlation > max_correlation:
                max_correlation = correlation


    for j in range(7):
        for i in range(j+1, 28):
            graph_1_average_1, graph_1_average_2 = calc_averages(libraries, j, j+7, i, i+7)
            graph_2_average_1, graph_2_average_2 = calc_averages(libraries, j, j+7, j+14, j+21)

            plot_scatter(graph_1_average_1, graph_1_average_2, "Number of counts of genes (in log) of the first condition",
                         graph_2_average_1, graph_2_average_2, "Number of counts of genes (in log) of the second condition",
                         "Well Fed", "Starvation")

def create_scattered_of_two_conditions():
    experiment_results = ExperimentResults()
    conditions = experiment_results.conditions
    experiment_results.delete_by_threshold(5)
    for i in range(7):
        graph_1_average_1, graph_1_average_2 = experiment_results.calc_averages_of_repetitions(i)
        graph_2_average_1, graph_2_average_2 = experiment_results.calc_averages_of_repetitions(j)
        label_x = "Num of counts (in log) in condition: {}".format(CONDITIONS[i])
        label_y = "Num of counts (in log) in condition: {}".format(CONDITIONS[j])

        plot_scatter(graph_1_average_1, graph_1_average_2, label_x, graph_2_average_1, graph_2_average_2, label_y,
                     CONDITIONS[i], CONDITIONS[j])


def gen_eyal_data_volcano():
    dir = "C:\\Users\\Danielle\\Documents\\project\\outputUnPairedAlignment\\dataBeforeIndexing\\after.TMM\\after.index\\"
    data = pd.read_csv(dir + "NegativeBinomial.csv", sep=",", header=None)
    pvalues = []
    log_fold_change = []
    for i in range(1, 15031):
        a = data[0][i].split(sep=";")[-1]
        b = data[1][i].split(sep=";")[0]
        pvalues.append(float(a + b))

        c = data[1][i].split(sep=";")[-1]
        d = data[2][i].split(sep=";")[0]
        log_fold_change.append(float(c + d))

    save_data_to_excel(np.array(log_fold_change), np.array(pvalues))
    plot_volcano_from_library()
    # minus_log_pvalues = -np.log(pvalues)
    # volcano_plot(log_fold_change, minus_log_pvalues)


def check_gene(gene_name):
    experiment_results = ExperimentResults()
    condition_pair_set_indeces = (cond_trn_minus, cond_well_fed)
    condition_pair_set = experiment_results.create_condition_pair_set_with_delete(condition_pair_set_indeces, threshold)
    gene_index = np.where(experiment_results.gene_transcript_id == gene_name)[0][0]
    if gene_index != 0:
        repetitions_cond_0 = condition_pair_set.conditions[0].repetitions[gene_index, :]
        repetitions_cond_1 = condition_pair_set.conditions[1].repetitions[gene_index, :]
        save_gene_plot(repetitions_cond_0, repetitions_cond_1, gene_name)


def create_bars_of_expression_in_conditions(gene_name):
    experiment_results = ExperimentResults()
    gene_index = np.where(experiment_results.gene_transcript_id == gene_name)[0][0]
    repetitions = np.zeros((condition_num, repetitions_num-1))
    for i in range(condition_num):
        rep_indeces = experiment_results.conditions[i].repetition_indeces
        # repetitions[i,:] = experiment_results.conditions[i].repetitions[gene_index, :]
        repetitions[i, :] = (experiment_results.cpm[gene_index, rep_indeces])[np.array([0,1,3])]
    create_eroor_bar(repetitions, gene_name)


def cascaded_filter_analysis():
    experiment_results = ExperimentResults()

    # find indeces of genes for which any of the libraries contain very low counts
    indeces_of_genes_not_valid = experiment_results.get_genes_under_threshold(threshold)

    # delete the indeces
    gene_count = experiment_results.delete_genes_by_index(indeces_of_genes_not_valid)

    initial_thresholds = (0.05, 1.2)

    final_thresholds = experiment_results.delete_genes_with_no_significant_change(
        cond_trn_minus, cond_well_fed, initial_thresholds)

    initial_thresholds = (0.05, 1.2)

    final_thresholds = experiment_results.delete_genes_with_no_significant_change(
        cond_moc_plus, cond_trn_plus, initial_thresholds)


if __name__ == '__main__':
    cascaded_filter_analysis()
    # run()
    # create_volcano()
    # create_genes_expression_plot()
    # get_volcano_std()
    # create_hist_of_std()
    # create_hist_of_foldchange_std()
    create_scattered_of_two_repetitions()
    # gen_eyal_data_volcano()
    # check_gene("T15B7.2.1")

    # create_bars_of_expression_in_conditions("F09E10.11a")
    # create_bars_of_expression_in_conditions("C05E4.9a")
