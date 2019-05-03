import numpy as np
import scipy as scipy
import scipy.stats
import itertools

import matplotlib.pyplot as plt
import pandas as pd
import seabornmaster.seaborn as sns

dir = "C:\\Users\\Danielle\\Documents\\project\\outputUnPairedAlignment\\dataBeforeIndexing\\after.TMM"

NUM_CONDITIONS = 7
NUM_REPETITIONS = 4
MILLION = 1000000

LIST_OF_CONDITIONS_REPETITIONS = ['WellFed.day1','WellFed.day2','WellFed.day3','WellFed.day4',
                                  'moc.day1','moc.day2','moc.day3','moc.day4',
                                  'trained.day1','trained.day2','trained.day3','trained.day4',
                                  'moc-.day1','moc-.day2','moc-.day3','moc-.day4',
                                  'moc+.day1','moc+.day2','moc+.day3','moc+.day4',
                                  'trained-.day1','trained-.day2','trained-.day3','trained-.day4',
                                  'trained+.day1','trained+.day2','trained+.day3','trained+.day4']

LIST_OF_CONDITIONS_REPETITIONS_AFTER_DELETE = ['WellFed.day1', 'WellFed.day3','WellFed.day4',
                                  'moc.day1','moc.day2','moc.day3','moc.day4',
                                  'trained.day1','trained.day3','trained.day4',
                                  'moc-.day1','moc-.day2','moc-.day3','moc-.day4',
                                  'moc+.day1','moc+.day2','moc+.day3','moc+.day4',
                                  'trained-.day1','trained-.day2','trained-.day3','trained-.day4',
                                  'trained+.day1','trained+.day2','trained+.day3','trained+.day4']
def test():
    """
    This test checks if the heat map works properly.
    :return:
    """
    a1 = np.array((
        [1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15,
         0.10, 0.05]))
    a2 = np.array(
        [0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10,
         0.05, 0.00])
    a3 = np.array(
        [0.95, 1.00, 0.85, 0.90, 1.00, 0.80, 0.65, 0.70, 0.55, 0.60, 0.45, 0.00, 0.35, 0.40, 0.25, 0.30, 0.15, 0.20,
         0.05, 0.10])
    a4 = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    a5 = np.array(
        [0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, -0.05, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45,
         0.45, 0.45])
    a6 = np.array(
        [0.45, 0.55, 0.45, 0.55, 0.45, 0.55, 0.45, 0.55, 0.45, 0.55, -0.05, 0.55, 0.45, 0.55, 0.45, 0.55, 0.45, 0.55,
         0.45, 0.55])
    a7 = np.array(
        [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
         0.90, 0.95])
    a8 = np.array(
        [-0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,
         0.85, 0.90])
    a9 = np.array(
        [-0.05, 0.10, 0.05, 0.20, 0.15, 0.30, 0.25, 0.40, 0.35, 0.50, 0.45, 0.60, 0.55, 0.70, 0.65, 0.80, 0.75, 0.90,
         0.85, 1.00])
    all_ = np.transpose(np.array([a1, a2, a3, a4, a5, a6, a7, a8, a9]))
    to_heat_map = np.zeros((all_.shape[1], all_.shape[1]))
    for i in range(all_.shape[1]):
        for j in range(all_.shape[1]):
            to_heat_map[i, j] = scipy.stats.pearsonr(all_[:, i], all_[:, j])[0]
            print("{}\t".format(to_heat_map[i, j]), end=" ")
        print()
    plt.imshow(to_heat_map, cmap='Blues', interpolation='nearest')
    plt.show()


def get_best_correlated(corelations):
    """
    The method calculates the 4 high correlated libraries for each library.
    :param corelations:
    :return:
    """
    correlations_arranged1 = np.zeros((corelations.shape[1], 4))
    correlations_arranged2 = np.zeros((corelations.shape[1], 4))
    for i in range(corelations.shape[1]):
        max_1 = [0, 0]
        max_2 = [0, 0]
        max_3 = [0, 0]
        max_4 = [0, 0]

        for j in range(corelations.shape[1]):
            correlation = scipy.stats.pearsonr(corelations[:, i], corelations[:, j])[0]
            if correlation > max_1[1]:
                if correlation > max_2[1]:
                    if correlation > max_3[1]:
                        if correlation > max_4[1]:
                            max_4 = [j, correlation]
                        else:
                            max_3 = [j, correlation]
                    else:
                        max_2 = [j, correlation]
                else:
                    max_1 = [j, correlation]

            correlations_arranged1[i, :] = max_1[0], max_2[0], max_3[0], max_4[0]
            correlations_arranged2[i, :] = max_1[1], max_2[1], max_3[1], max_4[1]

    print(correlations_arranged1, end=" ")
    print(correlations_arranged2, end=" ")
    return correlations_arranged1


def get_repetitions(condition_index):
    """
    This method yields all the repetition indeces of a condition.
    :param condition_index:
    :return:
    """
    for repetition_index in range(NUM_REPETITIONS):
        yield condition_index + repetition_index * NUM_CONDITIONS


def correlations_arranged_by_repitions(correlations):
    """
    The method arranges the correlation matrix to a matrix which is arranged as follows:
    rows (library) \ columns - 0,7,14,21, 1,8,15,22...
    instead of - rows (library) \ columns - 0,1,2,3,4,5...
    :param correlations:
    :return:
    """
    correlations_arranged = np.zeros((correlations.shape[1], correlations.shape[0]))
    a = []

    for i in range(NUM_CONDITIONS + 1):
        a.append(list(get_repetitions(i)))

    concated = list(itertools.chain.from_iterable(a))

    for column in range(28):
        library_index_to_add = concated[column]
        correlations_arranged[column, :] = np.transpose(correlations[:, library_index_to_add])

    return np.transpose(correlations_arranged)


def normalize_frequencies(frequencies):
    """
    Normalizes the frequencies (divide each gene_count with the sum of counts in the library, and
    multiply by million.
    :param frequencies:
    :return:
    """
    updated_frequencies = np.zeros((frequencies.shape[0], frequencies.shape[1]))
    for i in range(frequencies.shape[1]):
        sum = np.sum(frequencies[:, i])
        for j in range(frequencies.shape[0]):
            updated_frequencies[j, i] = (frequencies[j, i] / sum) * MILLION
    return updated_frequencies


def run(dir):
    frequencies = np.load(dir + "\\libraries.frequencies.results.after.index.npy").astype(float)
    # frequencies = pd.read_excel(dir + "\\allLibs.xlsx").astype(float)
    # frequencies = np.array(frequencies)
    correlations_arranged = normalize_frequencies(frequencies)
    correlations_arranged = correlations_arranged_by_repitions(correlations_arranged)
    create_heat_map(correlations_arranged)
    # create_heat_map(updated_frequencies)
    # bestCorrelated = get_best_correlated(correlations_arranged)


def delete_zero_genes(first_lib, second_lib):
    first_indeces = np.where(first_lib < 2)[0]
    second_indeces = np.where(second_lib < 2)[0]

    indeces = np.union1d(first_indeces, second_indeces)

    first_lib = np.delete(first_lib, indeces)
    second_lib = np.delete(second_lib, indeces)

    return np.log(first_lib), np.log(second_lib)


def write_to_excel_doc(values, file_name):
    values_df = pd.DataFrame(values)
    writer = pd.ExcelWriter(dir + file_name, engine='xlsxwriter')
    values_df.to_excel(writer, sheet_name='Sheet1')
    writer.save()


def create_heat_map(correlations_arranged):
    """
    Creates a heat Map.
    :param correlations_arranged:
    :return:
    """
    to_heat_map = np.zeros((correlations_arranged.shape[1], correlations_arranged.shape[1]))

    for i in range(correlations_arranged.shape[1]):
        for j in range(correlations_arranged.shape[1]):

            first_lib = correlations_arranged[:, i]
            second_lib = correlations_arranged[:, j]
            delete_zero_genes(first_lib, second_lib)
            a = scipy.stats.pearsonr(first_lib, second_lib)[0]

            to_heat_map[i, j] = a
            print(to_heat_map[i, j])
    
    write_to_excel_doc(to_heat_map, "\\correlations1.xlsx")

    # vactor = np.array(
    #     [0.9, 0.9, 0.9, 0.9, 0.72, 0.71, 0.74, 0.73, 0.72, 0.8, 0.65, 0.7, 0.73, 0.73, 0.72, 0.71, 0.69, 0.71, 0.75, 0.72, 0.71,
    #      0.72, 0.71, 0.73, 0.73, 0.71, 0.74, 0.71])
    #
    # to_heat_map[1, :] = vactor
    # to_heat_map[:, 1] = vactor

    # for i in range(to_heat_map.shape[1]):
    #     for j in range(to_heat_map.shape[1]):
    #         if (to_heat_map[i, j] < 0.6):
    #             to_heat_map = np.delete(to_heat_map, (j), axis=0)
    #             to_heat_map = np.delete(to_heat_map, (j), axis=1)

    # I deleted rows and columns 1,8 because the correlation was low
    # to_heat_map = np.delete(to_heat_map, (1), axis=0)
    # to_heat_map = np.delete(to_heat_map, (1), axis=1)
    # to_heat_map = np.delete(to_heat_map, (8), axis=0)
    # to_heat_map = np.delete(to_heat_map, (8), axis=1)

    # plt.imshow(to_heat_map, cmap='Blues', interpolation='nearest')
    to_heat_map_df = pd.DataFrame(to_heat_map)
    # to_heat_map_df.columns = LIST_OF_CONDITIONS_REPETITIONS_AFTER_DELETE
    cmap = sns.cubehelix_palette(dark=0, light=1, as_cmap=True)
    sns.heatmap(to_heat_map_df, vmin=0.3, vmax=1, cmap=cmap)
    # plt.imshow(to_heat_map, cmap='Blues', interpolation='nearest')
    plt.show()



if __name__ == "__main__":
    # test()
    run("C:\\Users\\Danielle\\Documents\\project\\outputUnPairedAlignment\\dataBeforeIndexing\\after.TMM")
    # run("C:\\Users\\Danielle\\Documents\\SharedUbuntu\\outputPairedAlignment\\data\\after.TMM")
