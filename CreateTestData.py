import numpy as np

num_genes = 6
num_conditions = 2
num_repetitions = 4


def create_test_volcaono():
    total_num_conditions = num_conditions * num_repetitions
    test_data = np.ones((num_genes, total_num_conditions))
    test_data += 100
    left_condition_indeces = np.arange(0, total_num_conditions, 2)
    right_condition_indeces = np.arange(1, total_num_conditions, 2)

    preapare_middle_high_gene(test_data, left_condition_indeces, right_condition_indeces, 1)
    preapare_middle_low_gene(test_data, left_condition_indeces, right_condition_indeces, 451)

    preapare_left_low_gene(test_data, left_condition_indeces, right_condition_indeces, 251)
    preapare_left_high_gene(test_data, left_condition_indeces, right_condition_indeces, 251)

    preapare_right_high_gene(test_data, left_condition_indeces, right_condition_indeces, 251)
    preapare_right_low_gene(test_data, left_condition_indeces, right_condition_indeces, 251)

    return test_data


def preapare_middle_high_gene(test_data, condition_0_indeces, condition_1_indeces, variance):
    test_data[0, condition_0_indeces[:2]] = -499
    test_data[0, condition_0_indeces[2:]] = 500

    test_data[0, condition_1_indeces[:2]] = 500
    test_data[0, condition_1_indeces[2:]] = -499


def preapare_middle_low_gene(test_data, condition_0_indeces, condition_1_indeces, variance):
    test_data[1, condition_0_indeces[:2]] = 20
    test_data[1, condition_0_indeces[2:]] = 20

    test_data[1, condition_1_indeces[:2]] = 20
    test_data[1, condition_1_indeces[2:]] = 20


def preapare_left_low_gene(test_data, condition_0_indeces, condition_1_indeces, variance):
    test_data[2, condition_0_indeces[:2]] = -3500
    test_data[2, condition_0_indeces[2:]] = 3800

    test_data[2, condition_1_indeces[:2]] = 3500
    test_data[2, condition_1_indeces[2:]] = 3500


def preapare_left_high_gene(test_data, condition_0_indeces, condition_1_indeces, variance):
    test_data[3, condition_0_indeces[:2]] = -300
    test_data[3, condition_0_indeces[2:]] = 500

    test_data[3, condition_1_indeces[:2]] = 1000
    test_data[3, condition_1_indeces[2:]] = 10


def preapare_right_high_gene(test_data, condition_0_indeces, condition_1_indeces, variance):
    test_data[4, condition_0_indeces[:2]] = 10000
    test_data[4, condition_0_indeces[2:]] = -300

    test_data[4, condition_1_indeces[:2]] = 500
    test_data[4, condition_1_indeces[2:]] = -100


def preapare_right_low_gene(test_data, condition_0_indeces, condition_1_indeces, variance):
    test_data[5, condition_0_indeces[:2]] = 10000
    test_data[5, condition_0_indeces[2:]] = 10000

    test_data[5, condition_1_indeces[:2]] = 500
    test_data[5, condition_1_indeces[2:]] = 500

def preapare_right_low_gene(test_data, condition_0_indeces, condition_1_indeces, variance):
    test_data[6, condition_0_indeces[:2]] = 110
    test_data[6, condition_0_indeces[2:]] = 90

    test_data[6, condition_1_indeces[:2]] = 1010
    test_data[6, condition_1_indeces[2:]] = 910


if __name__ == '__main__':
    test_data = create_test_volcaono()


    np.save("./data/libraries.test", test_data)
    test_genes = np.array(["gene1", "gene2", "gene3", "gene4", "gene5"])
    np.save("./data/gene.names.test", test_genes)

    # random_noise = np.random.normal(loc=0.0, scale=2.0, size=test_data.shape)
    # test_data += random_noise
