import numpy as np


def setGene(test_data, gene_index, library_index, value):
    test_data[gene_index, library_index] = value


if __name__ == '__main__':
    num_genes = 3
    num_conditions = 2
    num_repetitions = 4

    total_num_conditions = num_conditions * num_repetitions
    # test_data = np.zeros((num_genes, total_num_conditions))

    test_data = np.ones((num_genes, total_num_conditions))
    test_data = test_data * 100

    # random_noise = np.random.normal(loc=0.0, scale=2.0, size=test_data.shape)
    # test_data += random_noise
    test_data[0, 0] = 1000
    test_data[0, 2] = 1000
    test_data[0, 4] = 1000
    test_data[0, 6] = 1000

    test_data[2, 0] = 10000
    test_data[2, 2] = 10000
    test_data[2, 4] = 10000
    test_data[2, 6] = 10000

    # condition_0_indeces = np.arange(0, total_num_conditions, 2)
    # condition_1_indeces = np.arange(1, total_num_conditions, 2)
    #
    # test_data[0, condition_0_indeces] = 900
    # test_data[1, condition_0_indeces] = 8
    # test_data[2, condition_0_indeces] = 10
    #
    # test_data[0, condition_1_indeces] = 8
    # test_data[1, condition_1_indeces] = 900
    # test_data[2, condition_1_indeces] = 10
    #
    # # test_data = np.transpose(np.array(test_data))
    np.save("./data/test_frequencies", test_data)
    test_genes = np.array(["gene1", "gene2", "gene3"])
    np.save("./data/test_genes", test_genes)

    # random_noise = np.random.normal(loc=0.0, scale=2.0, size=test_data.shape)
    # test_data += random_noise
