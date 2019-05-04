from ExperimentResults import ExperimentResults
from ShowResults import *

cond_moc_minus = 3
cond_moc_plus = 4
cond_trn_minus = 5
cond_trn_plus = 6
threshold = 5


def get_condition_indeces_to_test():
    return [cond_moc_plus, cond_moc_minus, cond_trn_plus, cond_trn_minus]


def run():
    experiment_results = ExperimentResults()

    indeces = get_condition_indeces_to_test()

    test_condition_set = experiment_results.create_test_condition_set(indeces, threshold)

    ratio_moc, ratio_trn = test_condition_set.calculate_fold_changes()

    show_results(experiment_results.gene_names, ratio_moc, ratio_trn)


def create_volcano():
    experiment_results = ExperimentResults()

    condition_pair_set_indeces = (cond_moc_plus, cond_moc_minus)

    condition_pair_set = experiment_results.create_condition_pair_set(condition_pair_set_indeces, threshold)

    pValues = condition_pair_set.get_p_value_pair()
    fold_changes = condition_pair_set.calculate_fold_changes()

    volcano_plot(fold_changes, -np.log(pValues))


if __name__ == '__main__':
    # run()
    create_volcano()
