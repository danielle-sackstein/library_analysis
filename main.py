from ExperimentResults import ExperimentResults
from ShowResults import *
import pandas as pd

from geneDiffCalculator import geneDiffCalculator

cond_moc_minus = 3
cond_moc_plus = 4
cond_trn_minus = 5
cond_trn_plus = 6
threshold = 5

directory = "./data/"


def save_data_to_excel(fold_changes, pValues):
    data = pd.DataFrame(columns=['fold_change', 'pValues'])
    data['fold_change'] = fold_changes
    data['pValues'] = pValues
    data.to_csv(directory + "to_volcano.csv", sep=",")


def get_condition_indeces_to_test():
    return [cond_trn_plus, cond_moc_plus, cond_trn_minus, cond_moc_minus]


def run():
    experiment_results = ExperimentResults()

    indeces = get_condition_indeces_to_test()

    test_condition_set = experiment_results.create_test_condition_set(indeces, threshold)

    ratio_moc, ratio_trn = test_condition_set.calculate_fold_changes()

    show_results(experiment_results.gene_names, ratio_moc, ratio_trn)


def create_volcano():
    experiment_results = ExperimentResults()

    # condition_pair_set_indeces = (cond_trn_plus, cond_trn_minus)
    condition_pair_set_indeces_test = (0,1)

    condition_pair_set = experiment_results.create_condition_pair_set(condition_pair_set_indeces_test, threshold)

    create_volcano_plot(condition_pair_set, experiment_results)

    # minus_log_pvalues = -np.log(pValues)
    # indeces = np.where(minus_log_pvalues > 1000)
    # minus_log_pvalues[indeces] = 160
    # volcano_plot(log_fold_change, minus_log_pvalues)


def create_volcano_plot(condition_pair_set, experiment_results):
    pValues = condition_pair_set.get_p_values()
    fold_changes = condition_pair_set.calculate_fold_changes()
    log2_fold_change = np.log2(fold_changes)
    save_data_to_excel(log2_fold_change, pValues)
    plot_volcano_from_library()
    find_and_print_de_genes(experiment_results)


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

    for i in range(log_fold_changes.shape[0]):

        if log_fold_changes[i] <= -0.5:
            gene_name = gene_names[i]
            down_regulated_genes.append(gene_name)

        elif log_fold_changes[i] >= 0.85:
            gene_name = gene_names[i]
            up_regulated_genes.append(gene_name)

    print_de_genes(down_regulated_genes, up_regulated_genes)


def create_genes_expression_plot():
    experiment_results = ExperimentResults()

    condition_pair_set_indeces = (cond_trn_plus, cond_trn_minus)

    condition_pair_set = experiment_results.create_condition_pair_set(condition_pair_set_indeces, threshold)
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


def get_std():
    experiment_results = ExperimentResults()
    condition_pair_set_indeces = (cond_trn_plus, cond_moc_plus)
    condition_pair_set = experiment_results.create_condition_pair_set(condition_pair_set_indeces, threshold)
    condition_pair_set.delete_by_std()
    condition_pair_set.delete_by_threshold(3)
    create_volcano_plot(condition_pair_set, experiment_results)


if __name__ == '__main__':
    # run()
    create_volcano()
    # create_genes_expression_plot()
    # get_std()
