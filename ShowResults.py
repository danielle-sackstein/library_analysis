import numpy as np
import matplotlib.pyplot as plt

CONDITIONS = ["Well Fed", "mock", "trained", "mock(-)", "mock(+)", "trained(-)", "trained(+)"]

def get_outlier_indeces(x, y):
    x_threshold = 1.5
    y_threshold = 0.8

    compare = x / y
    indeces_x_outliers = np.argwhere(compare > x_threshold)
    indeces_y_outliers = np.argwhere(compare < y_threshold)
    return indeces_x_outliers, indeces_y_outliers


def label_genes(ax1, gene_names, indeces_x_outliers, x, y):
    for i in indeces_x_outliers:
        ax1.annotate(gene_names[i], (x[i], y[i]), size=5, bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5))


def plot_scatter(graph_1_lib_1, graph_1_lib_2, label_x, graph_2_lib_1, graph_2_lib_2, labal_y, left_condition,
                 right_condition):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.scatter(graph_1_lib_1, graph_1_lib_2, s=10, c='red',
                label="{} X {} ".format(left_condition, left_condition))
    ax1.scatter(graph_2_lib_1, graph_2_lib_2, s=10, c='black',
                label="{} X {} ".format(right_condition, right_condition))
    ax1.legend(loc='lower right')

    plt.title('Difference in genes count of two {} and {}'.format(left_condition, right_condition))

    plt.xlabel(label_x)
    plt.ylabel(labal_y)
    plt.savefig('scattered_diff/Scatter.Plot.of{}.and.{}.png'.format(left_condition, right_condition))

    plt.show()
    plt.clf()


def show_results(gene_names, ratio_plus, ratio_minus):
    title_moc = "log(moc plus) / log(moc minus)"
    title_trn = "log(trained plus) / log(trained minus)"

    title_plus = "log(trained plus) / log(moc plus)"
    title_minus = "log(trained minus) / log(moc minus)"

    # indeces_x_outliers, indeces_y_outliers = get_outlier_indeces(ratio_plus, ratio_minus)

    plot_scatter(
        "ratio of averages in conditions",
        ratio_plus, title_minus,
        ratio_minus, title_plus,
        gene_names
    )

    print_results(gene_names, indeces_x_outliers, indeces_y_outliers)


def print_results(gene_names, indeces_x_outliers, indeces_y_outliers):
    print("Genes that are expressed in moc more than in trained")
    for i in indeces_x_outliers:
        print(gene_names[i])
    print("Genes that are expressed in trained more than in moc")
    for i in indeces_y_outliers:
        print(gene_names[i])


def volcano_plot(fold_changes, pValues):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.scatter(fold_changes, pValues, s=10, c='red')

    # ax1.legend(loc='lower right')
    # ax1.set_xlim(-10,10)
    # ax1.set_ylim(-1,1)

    plt.title("volcano Plot")
    plt.xlabel("log(fold_changes)")
    plt.ylabel("-log(pvalue)")

    plt.savefig('Scatter Plot')
    plt.show()
    plt.clf()


def save_gene_plot(x, y, i):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # ax1.plot([1, 2, 3, 4], x, label="Trained")
    # ax1.plot([1, 2, 3, 4], y, label="Mock Trained")
    ax1.plot([1, 2, 3, 4], x, label="Starved")
    ax1.plot([1, 2, 3, 4], y, label="Well Fed")
    ax1.legend(loc='best')
    plt.xlabel("days")
    plt.ylabel("CPM")

    # plt.show()
    plt.savefig("starved.unstarved.genes/gene.plot.{}.png".format(i))


def plot_hist_of_std(condition_std, cond_name):
    range = np.arange(condition_std.shape[0])
    sorted = np.sort(condition_std)
    plt.bar(range, sorted, align='center', alpha=0.5)
    plt.xlabel("gene num")
    plt.ylabel("(std \ averages) gene ")
    plt.title("(std \ averages) vs gene of condition {}".format(cond_name))
    plt.savefig("histogram.of.{}.png".format(cond_name))
    plt.clf()
    # plt.show()


def plot_hist_of_foldchange_std(left_condition_std, right_condition_std, left_cond_name, right_cond_name):
    left_condition_average = np.average(left_condition_std)
    right_condition_average = np.average(right_condition_std)
    average_diff = abs(right_condition_average - left_condition_average)

    std_divided_with_avr = (left_condition_std + right_condition_std) / average_diff
    sorted = np.sort(std_divided_with_avr)
    sorted = sorted[np.where(sorted < 2)]

    print(sorted.shape)
    range = np.arange(sorted.shape[0])

    plt.bar(range, sorted, align='center', alpha=0.5)
    plt.xlabel("gene num")
    plt.ylabel("(std1 + st2) \ (|average1 - average2|)")
    plt.title("(std1 + st2) \ (|average1 - average2|)")
    plt.savefig("new_histogram.png")
    plt.clf()
    # plt.show()


from matplotlib.axes import Axes

from pylab import *

def create_eroor_bar(repetitions, gene_name):
    means = np.mean(repetitions, axis=1)
    std = np.std(repetitions, axis=1)
    x = np.arange(1,8)
    colors = ["salmon", "indianred","indianred", "darkred", "darkred","darkred","darkred"]
    legend = ["Well Fed", "Recovery","Recovery", "Starved", "Starved","Starved","Starved"]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for i in range(7):
        if i == 0:
            ax1.bar(x[i], means[i],  color=colors[i], yerr=np.array(std)[i], label=legend[i])
        if i == 1:
            ax1.bar(x[i], means[i],  color=colors[i], yerr=np.array(std)[i], label=legend[i])
        if i == 3:
            ax1.bar(x[i], means[i],  color=colors[i], yerr=np.array(std)[i], label=legend[i])
        else:
            ax1.bar(x[i], means[i],  color=colors[i], yerr=np.array(std)[i])

    ax1.legend(loc='upper left')
    plt.title("Cpm of each condition of significant in starvation gene {}".format(gene_name))
    plt.xlabel("Condition")
    plt.ylabel("Cpm")
    # plt.show()
    plt.savefig("./interesting.genes/gene {} expressions.png".format(gene_name))
    print(repetitions)


from bioinfokit.bioinfokit import visuz


def plot_volcano_from_library():
    visuz.volcano(table="./data/to_volcano.csv", lfc="fold_change", pv="pValues")
    visuz.involcano(table="./data/to_volcano.csv", lfc="fold_change", pv="pValues")


