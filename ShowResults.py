import numpy as np
import matplotlib.pyplot as plt


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


def plot_scatter(title, x, label_x, y, label_y, gene_names, indeces_x_outliers, indeces_y_outliers) :

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.scatter(x, y, s=10, c='red', label="Gene expression")
    ax1.scatter(x[indeces_x_outliers], y[indeces_x_outliers], s=10, c='blue', label="x outliers")
    ax1.scatter(x[indeces_y_outliers], y[indeces_y_outliers], s=10, c='yellow', label="y outliers")

    label_genes(ax1, gene_names, indeces_x_outliers, x, y)
    label_genes(ax1, gene_names, indeces_y_outliers, x, y)

    ax1.legend(loc='lower right')

    plt.title('Difference in genes count of two conditions ' + title)
    plt.xlabel(label_x)
    plt.ylabel(label_y)

    plt.savefig('Scatter Plot')
    plt.show()
    plt.clf()


def show_results(gene_names, ratio_moc, ratio_trn):

    title_moc = "log(moc plus) / log(moc minus)"
    title_trn = "log(trained plus) / log(trained minus)"

    indeces_x_outliers, indeces_y_outliers = get_outlier_indeces(ratio_moc, ratio_trn)

    plot_scatter(
        "ratio of averages",
        ratio_moc, title_moc,
        ratio_trn, title_trn,
        gene_names,
        indeces_x_outliers, indeces_y_outliers)

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

    ax1.scatter(fold_changes, pValues, s=10, c='red', label="Gene expression")

    ax1.legend(loc='lower right')

    plt.title('Difference in genes count of two conditions ' + title)
    plt.xlabel("fold_changes")
    plt.ylabel("pvalues")

    plt.savefig('Scatter Plot')
    plt.show()
    plt.clf()
