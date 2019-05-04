import numpy as np

if __name__ == '__main__':
    test_data = np.ones((100, 8))
    test_data = test_data * 100

    random_noise = np.random.normal(loc=0.0, scale=2.0, size=test_data.shape)
    test_data += random_noise
    test_data[0, 0] = 1000
    test_data[0, 2] = 1000
    test_data[0, 4] = 1000
    test_data[0, 6] = 1000

    test_data[2, 0] = 10000
    test_data[2, 2] = 10000
    test_data[2, 4] = 10000
    test_data[2, 6] = 10000

    test_data = np.transpose(np.array(test_data))
    np.save("./data/test_frequencies", test_data)
    test_genes = np.array(["gene1", "gene2", "gene3"])
    np.save("./data/test_genes", test_genes)
