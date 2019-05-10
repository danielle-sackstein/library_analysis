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


def plot_scatter(title, x, label_x, y, label_y, gene_names) :

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.scatter(x, y, s=10, c='red', label="gene expression in each ratio of conditions")
    # ax1.scatter(x[indeces_x_outliers], y[indeces_x_outliers], s=10, c='blue', label="x outliers")
    # ax1.scatter(x[indeces_y_outliers], y[indeces_y_outliers], s=10, c='yellow', label="y outliers")

    # label_genes(ax1, gene_names, indeces_x_outliers, x, y)
    # label_genes(ax1, gene_names, indeces_y_outliers, x, y)

    ax1.legend(loc='best')

    plt.title(title)
    plt.xlabel(label_x)
    plt.ylabel(label_y)

    plt.savefig('Scatter Plot')
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

def save_gene_plot(x, y,i):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot([1,2,3,4],x, label="trained plus")
    ax1.plot([1,2,3,4],y, label="trained minus")
    ax1.legend(loc='best')
    plt.xlabel("days")
    plt.ylabel("expression")

    # plt.show()
    plt.savefig("gene.plots/gene.plot.{}.png".format(i))


from bioinfokit.bioinfokit import visuz

def plot_volcano_from_library():
    visuz.volcano(table="./data/to_volcano.csv", lfc="fold_change", pv="pValues")
    visuz.involcano(table="./data/to_volcano.csv", lfc="fold_change", pv="pValues")
