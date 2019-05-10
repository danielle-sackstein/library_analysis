class geneDiffCalculator:
    def __init__(self, first_condition, second_condition):
        self.first_condition = first_condition
        self.second_condition = second_condition

    def create_genes_lists(self, gene_index):
        first_cond_gene_repetitions = self.first_condition.repetitions[gene_index]
        second_cond_gene_repetitions = self.second_condition.repetitions[gene_index]
        return first_cond_gene_repetitions, second_cond_gene_repetitions
