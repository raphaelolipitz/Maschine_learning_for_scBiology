import numpy as np
import pandas as pd
from scipy import io, stats
import seaborn as sns
from matplotlib import pyplot as plt

# Convenience functions to handle data import
def import_expression_data(mtx_gene_exp_fil, barcodes_fil, features_fil):

    gene_expression = io.mmread(mtx_gene_exp_fil)
    barcodes = pd.read_csv(barcodes_fil, header=None)[0]
    features = pd.read_csv(features_fil, header=None)[0]

    gene_expression = pd.DataFrame(gene_expression.toarray(), index=barcodes, columns=features)

    return gene_expression

# Insert your code into the function shells below
def task_1_a(relative_path_to_data = '../exercise_1_data/'):

    gene_expression = import_expression_data(relative_path_to_data+'expression_data_1.mtx', 
                                             relative_path_to_data+'expression_data_1_barcodes.tsv',
                                             relative_path_to_data+'expression_data_1_genes.tsv')
    cell_annotations = pd.read_csv(relative_path_to_data+'expression_data_1_metadata.tsv', sep='\t', index_col=0)

    organ_identity = {'A': '', 'B': '', 'C': '', 'D': '', 'E': '', 'F': ''}
    cell_type_identity = {'T1': '', 'T2': '', 'T3': '', 'T4': '', 'T5': ''}

    '''Insert your code here and update the dictionaries with your identification. 
    Take care to use the same annotations are used in the paper.'''

    organ_identity = pd.Series(organ_identity).to_frame()
    cell_type_identity = pd.Series(cell_type_identity).to_frame()

    return organ_identity, cell_type_identity

def task_1_b(relative_path_to_data = '../exercise_1_data/'):

    gene_expression = import_expression_data(relative_path_to_data+'expression_data_2.mtx', 
                                             relative_path_to_data+'expression_data_2_barcodes.tsv',
                                             relative_path_to_data+'expression_data_2_genes.tsv')
    
    cell_annotations = pd.DataFrame(index=gene_expression.index.values, columns=['phenotype'])
    
    '''Insert code here and update the cell_annotations dataframe with the phenotypes.
       Take care to use the same annotations are used in the paper.'''

    return cell_annotations

def task_2_a(relative_path_to_data = '../exercise_1_data/'):

    cell_sequences = np.load(relative_path_to_data+'tex_sampling.npy')
    cell_annotations = pd.read_csv(relative_path_to_data+'expression_data_3_metadata.txt', sep='\t', index_col=0)

    upper_trajectory_louvain_ordering = []
    lower_trajectory_louvain_ordering = []

    '''Insert code here and append each list with the sequence of louvain clusters.'''

    return upper_trajectory_louvain_ordering, lower_trajectory_louvain_ordering

def task_2_b(relative_path_to_data = '../exercise_1_data/'):

    gene_expression = import_expression_data(relative_path_to_data+'expression_data_3.mtx', 
                                             relative_path_to_data+'expression_data_3_barcodes.tsv',
                                             relative_path_to_data+'expression_data_3_genes.tsv')
    
    cell_annotations = pd.read_csv(relative_path_to_data+'expression_data_3_metadata.txt', sep='\t', index_col=0)

    gene_correlations = pd.DataFrame(index=gene_expression.columns.values, columns=['correlation', 'p-values'])

    '''Insert your code here.'''

    return gene_correlations

# Run the code after updating all the functions
def main(relative_path_to_data = '../exercise_1_data/', output_path=''):

    #Task 1A output
    t1_organ_identity, t1_cell_type_identity = task_1_a(relative_path_to_data)
    t1_organ_identity.to_csv(output_path+'t1_organ_identity.txt', header=False)
    t1_cell_type_identity.to_csv(output_path+'t1_celltype_identity.txt', header=False)

    t1b_cell_annotations = task_1_b(relative_path_to_data)
    t1b_cell_annotations.to_csv(output_path+'t1b_cell_annotations.txt', sep='\t')

    t2a_upper_trajectory_louvain_ordering, t2a_lower_trajectory_louvain_ordering = task_2_a(relative_path_to_data)
    pd.DataFrame(t2a_upper_trajectory_louvain_ordering).to_csv(output_path+'t2a_upper_louvain_sequence.txt', header=False)
    pd.DataFrame(t2a_lower_trajectory_louvain_ordering).to_csv(output_path+'t2a_lower_louvain_sequence.txt', header=False)

    t2b_gene_correlations = task_2_b(relative_path_to_data)
    t2b_gene_correlations.sort_values('correlation').iloc[:10].to_csv(output_path+'t2b_gene_correlations.txt', sep='\t')

    sns.scatterplot(t2b_gene_correlations.correlation, -np.log(t2b_gene_correlations['p-values']))
    plt.savefig('t2b_volcano_plot.png')

    return t1_cell_type_identity, t1_cell_type_identity, t1b_cell_annotations, t2a_upper_trajectory_louvain_ordering, t2a_lower_trajectory_louvain_ordering, t2b_gene_correlations

if __name__ == "__main__":

    # Update paths if required
    t1_cell_type_identity, t1_cell_type_identity, t1b_cell_annotations, t2a_upper_trajectory_louvain_ordering, 
    t2a_lower_trajectory_louvain_ordering, t2b_gene_correlations = main('../exercise_1_data/', './')





    




