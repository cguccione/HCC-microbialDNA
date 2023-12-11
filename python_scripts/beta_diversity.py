'''
Helper functions used to calculate basic
beta diversity
'''

import pandas as pd
import subprocess
import qiime2 as q2
from qiime2 import Artifact, Metadata, Visualization
from qiime2.plugins.feature_table.methods import rarefy
from qiime2.plugins import diversity
from qiime2.plugins.diversity_lib.methods import weighted_unifrac, unweighted_unifrac
from qiime2.plugins.diversity.visualizers import alpha_rarefaction
import biom
from skbio.stats.distance import permanova

#WOL2 Taxonomy
wol2_taxonomy = q2.Artifact.import_data('Phylogeny[Rooted]', '/Users/cguccion/Dropbox/Storage/HelpfulLabDocs/taxonomy_trees/WOL2/tree.nwk')

def beta_decoide(ft, meta, metric, permutations, decoide_min_feature_count):
    
    '''Calculates beta diversity for unweighted 
    & weighted unifraq

    Parameters
    ---------
    ft: q2 FeatureTable[Frequency] object
        Qiime2 taxonomy object
    
    meta: pd df
        metadata
        
    metric: str
        The column in sample_meta used to
        calculate beta diversity across 
        samples
    
    Returns
    -------
    beta_result_d: q2 Visulization object
        Qiime2 visulaization object of 
        decoide
    
    Notes
    -----
    '''
    
    #Convert meta to q2 object
    sample_meta = q2.Metadata(meta.set_index('sample_name'))
    
    #'''
    #Save ft for command line processing
    ft.save('deicode_processing.qza')
    
    #Run deicode from command line using qza created above
    #!qiime deicode rpca --i-table deicode_processing.qza --p-min-sample-count 500 --p-min-feature-count $decoide_min_feature_count --o-biplot deicode_biplot.qza --o-distance-matrix deicode_distance_test.qza
    
    command = [
    'qiime',
    'deicode',
    'rpca',
    '--i-table',
    'deicode_processing.qza',
    '--p-min-sample-count',
    '500',
    '--p-min-feature-count',
    str(decoide_min_feature_count),
    '--o-biplot',
    'deicode_biplot.qza',
    '--o-distance-matrix',
    'deicode_distance_test.qza']

    # Run the command
    subprocess.run(command)
    
    #Import biplot back into python
    rpca_biplot = Artifact.load('deicode_biplot.qza')
    
    #Import biplot back into python
    rpca_distance_matrix = Artifact.load('deicode_distance_test.qza')
    
    #Calculate permanova 
    beta_result_d = diversity.actions.beta_group_significance(distance_matrix=rpca_distance_matrix,
                                                              metadata=sample_meta.get_column(metric),
                                                              method = 'permanova', pairwise = True,
                                                              permutations=permutations)
                                                   
    return(beta_result_d)  

def table_prep(table, rarefaction=None):

    '''Convert data into qiime2 objects
    and rarefy if needed

    Parameters
    ---------
    table: pandas df
        Taxonomy table of samples and 
        corresponding microbial abudances
    
    rarefaction(optional): int
        Level to rarefy table to
    
    Returns
    -------
    ft: q2 FeatureTable[Frequency] object
        Qiime2 taxonomy object
    
    Notes
    -----
    '''
    
    #Convert table to q2 object
    ft = Artifact.import_data("FeatureTable[Frequency]", table.T)
    
    #(Optional) Rarefaction
    if rarefaction != None:
        ft = rarefy(table=ft, sampling_depth = rarefaction)
        ft = ft.rarefied_table
    
    return(ft)

def beta_unifraq(ft, sample_meta, metric, permutations):
    
    '''Calculates beta diversity for unweighted 
    & weighted unifraq

    Parameters
    ---------
    ft: q2 FeatureTable[Frequency] object
        Qiime2 taxonomy object
    
    sample_meta: q2 Metadata Object
        Qiime2 metadata object
        
    metric: str
        The column in sample_meta used to
        calculate beta diversity across 
        samples
    
    Returns
    -------
    beta_result_uw: q2 Visulization object
        Qiime2 visulaization object of 
        unweighted unifraq
    
    beta_result_w: q2 Visulization object
        Qiime2 visulaization object of 
        weighted unifraq
    
    Notes
    -----
    '''
    
    #Create Un-Weighted Unifraq distance matrix
    unweighted_unifrac_distance_matrix, = unweighted_unifrac(table=ft, phylogeny = wol2_taxonomy)

    #Calculate Beta Diversity: un-weighted unifraq
    beta_result_uw = diversity.actions.beta_group_significance(distance_matrix=unweighted_unifrac_distance_matrix,
                                                              metadata=sample_meta.get_column(metric),
                                                              method = 'permanova', pairwise = True,
                                                              permutations = permutations)
    
    #Create weighted Unifraq distance matrix
    weighted_unifrac_distance_matrix, = weighted_unifrac(table=ft, phylogeny = wol2_taxonomy)

    #Calculate Beta Diversity: uweighted unifraq
    beta_result_w = diversity.actions.beta_group_significance(distance_matrix=weighted_unifrac_distance_matrix,
                                                              metadata=sample_meta.get_column(metric),
                                                              method = 'permanova', pairwise = True,
                                                             permutations = permutations)
    
    return(beta_result_uw, beta_result_w)

def all_beta(table_genome_rs210, table_genome_wol2, meta, rarefaction, metric, permutations=999, numRares=10, decoide_min_feature_count=10):
    
    '''Convert data into qiime2 objects
    and rarefy if needed

    Parameters
    ---------
    table: pandas df
        Taxonomy table of samples and 
        corresponding microbial abudances
        
    meta: pandas df
        Metadata for samples in table
    
    rarefaction(optional): int
        Level to rarefy table to
    
    Returns
    -------
    ft: q2 FeatureTable[Frequency] object
        Qiime2 taxonomy object
    
    sample_meta: q2 Metadata Object
        Qiime2 metadata object
    
    Notes
    -----
    '''
    
    #Convert meta to q2 object
    sample_meta = q2.Metadata(meta.set_index('sample_name'))
    
    ##Create all Q2 Tables needed for Calculations##
    
    #RS210 Non-Rare
    ft_genome_rs210 = table_prep(table_genome_rs210)
    
    #WOL2 Non-Rare
    ft_genome_wol2 = table_prep(table_genome_wol2)
    
    ##Call all combinations of Beta Diversity##
    
    rs210_rpca_genome = beta_decoide(ft_genome_rs210, meta, metric, permutations, decoide_min_feature_count)
    
    wol2_rpca_genome = beta_decoide(ft_genome_wol2, meta, metric, permutations, decoide_min_feature_count)
    
    #Rarifed samples should be run 10x times and averaged
    
    #Create df to store averages across runs
    df_wol2_uwUni_genome = pd.DataFrame(columns = ['p-value', 'pseudo-F', 'Sample Size'])
    df_wol2_wUni_genome = pd.DataFrame(columns = ['p-value', 'pseudo-F', 'Sample Size'])
    df_wol2_rpca_rare_genome = pd.DataFrame(columns = ['p-value', 'pseudo-F', 'Sample Size'])
    df_rs210_rpca_rare_genome = pd.DataFrame(columns = ['p-value', 'pseudo-F', 'Sample Size'])
        
    for i in range(0, numRares):
        df_wol2_uwUni_genome, df_wol2_wUni_genome, df_wol2_rpca_rare_genome, df_rs210_rpca_rare_genome = rare_multiplier(df_wol2_uwUni_genome,
                                                                                                                         df_wol2_wUni_genome,
                                                                                                                         df_wol2_rpca_rare_genome,
                                                                                                                         df_rs210_rpca_rare_genome, 
                                                                                                                         table_genome_rs210, 
                                                                                                                         table_genome_wol2,
                                                                                                                         meta, rarefaction, 
                                                                                                                         metric, decoide_min_feature_count,
                                                                                                                         permutations=999)

    print('Wol2 unWeighted Unifrac Genome')
    display(df_wol2_uwUni_genome)
    print(df_wol2_uwUni_genome.mean())
    print()
    print('Wol2 Weighted Unifrac Genome')
    display(df_wol2_wUni_genome)
    print(df_wol2_wUni_genome.mean())
    print()
    print('Wol2 RPCA Rare Genome')
    display(df_wol2_rpca_rare_genome)
    print(df_wol2_rpca_rare_genome.mean())
    print()
    print('RS210 RPCA Rare Genome')
    display(df_rs210_rpca_rare_genome)
    print(df_rs210_rpca_rare_genome.mean())
    
    '''
    #RS210 Rare - adding in to rarify across RPCA
    ft_genome_rs210_rare = table_prep(table_genome_rs210, rarefaction=rarefaction)
    
    #WOL2 Rare
    ft_genome_wol2_rare = table_prep(table_genome_wol2, rarefaction=rarefaction)
    
    wol2_uwUni_genome, wol2_wUni_genome = beta_unifraq(ft_genome_wol2_rare, sample_meta, metric, permutations)
    
    #Rarefied RPCA
    wol2_rpca_rare_genome = beta_decoide(ft_genome_wol2_rare, sample_meta, metric, permutations)
    rs210_rpca_rare_genome = beta_decoide(ft_genome_rs210_rare, sample_meta, metric, permutations)
    '''
    
    #return(rs210_rpca_genome, wol2_uwUni_genome, wol2_wUni_genome, wol2_rpca_genome, wol2_rpca_rare_genome, rs210_rpca_rare_genome)
    return(rs210_rpca_genome, wol2_rpca_genome, df_wol2_uwUni_genome, df_wol2_wUni_genome, df_wol2_rpca_rare_genome, df_rs210_rpca_rare_genome)
    
def rare_helper(vis_output, df):
    
    #Save output
    vis_output.visualization.export_data("temp_vis_output")
    
    #Import vis output
    perm = pd.read_csv('temp_vis_output/permanova-pairwise.csv')
    
    #Extract values
    p_value = perm['p-value'][0]
    pseudoF = perm['pseudo-F'][0]
    size = perm['Sample size'][0]
    
    #Add new column to pandas df
    df.loc[len(df)] = [p_value, pseudoF, size]
    
    return(df)

def rare_multiplier(df_wol2_uwUni_genome, df_wol2_wUni_genome, df_wol2_rpca_rare_genome,
                    df_rs210_rpca_rare_genome, table_genome_rs210, table_genome_wol2,
                    meta, rarefaction, metric, decoide_min_feature_count, permutations=999):
    
    #Convert meta to q2 object
    sample_meta = q2.Metadata(meta.set_index('sample_name'))
    
    #Rare Tables
    ft_genome_rs210_rare = table_prep(table_genome_rs210, rarefaction=rarefaction)
    ft_genome_wol2_rare = table_prep(table_genome_wol2, rarefaction=rarefaction)
    
    #Unifrac Metrics
    wol2_uwUni_genome, wol2_wUni_genome = beta_unifraq(ft_genome_wol2_rare, sample_meta, metric, permutations)
    
    df_wol2_uwUni_genome = rare_helper(wol2_uwUni_genome, df_wol2_uwUni_genome)
    df_wol2_wUni_genome = rare_helper(wol2_wUni_genome, df_wol2_wUni_genome)
    
    #Rarefied RPCA
    wol2_rpca_rare_genome = beta_decoide(ft_genome_wol2_rare, meta, metric, permutations, decoide_min_feature_count)
    rs210_rpca_rare_genome = beta_decoide(ft_genome_rs210_rare, meta, metric, permutations, decoide_min_feature_count)
    
    df_wol2_rpca_rare_genome = rare_helper(wol2_rpca_rare_genome, df_wol2_rpca_rare_genome)
    df_rs210_rpca_rare_genome = rare_helper(rs210_rpca_rare_genome, df_rs210_rpca_rare_genome)
    
    return(df_wol2_uwUni_genome, df_wol2_wUni_genome, df_wol2_rpca_rare_genome, df_rs210_rpca_rare_genome)

def rare_multiplier_species(df_wol2_rpca_rare_genome, df_rs210_rpca_rare_genome, table_genome_rs210, table_genome_wol2,
                    meta, rarefaction, metric, decoide_min_feature_count, permutations=999):
    
    #Convert meta to q2 object
    sample_meta = q2.Metadata(meta.set_index('sample_name'))
    
    #Rare Tables
    ft_genome_rs210_rare = table_prep(table_genome_rs210, rarefaction=rarefaction)
    ft_genome_wol2_rare = table_prep(table_genome_wol2, rarefaction=rarefaction)
    
    
    #Rarefied RPCA
    wol2_rpca_rare_genome = beta_decoide(ft_genome_wol2_rare, meta, metric, permutations, decoide_min_feature_count)
    rs210_rpca_rare_genome = beta_decoide(ft_genome_rs210_rare, meta, metric, permutations, decoide_min_feature_count)
    
    df_wol2_rpca_rare_genome = rare_helper(wol2_rpca_rare_genome, df_wol2_rpca_rare_genome)
    df_rs210_rpca_rare_genome = rare_helper(rs210_rpca_rare_genome, df_rs210_rpca_rare_genome)
    
    return(df_wol2_rpca_rare_genome, df_rs210_rpca_rare_genome)

def alpha_rare_curve(table, max_depth, meta):
    
    q2_meta = Metadata(meta) 
    
    alph_vis = alpha_rarefaction(table = table, max_depth = max_depth, metadata = q2_meta)
    
    return(alph_vis)


def all_beta_species(table_genome_rs210, table_genome_wol2, meta, rarefaction, metric, permutations=999, numRares=10, decoide_min_feature_count=10):
    
    '''Convert data into qiime2 objects
    and rarefy if needed

    Parameters
    ---------
    table: pandas df
        Taxonomy table of samples and 
        corresponding microbial abudances
        
    meta: pandas df
        Metadata for samples in table
    
    rarefaction(optional): int
        Level to rarefy table to
    
    Returns
    -------
    ft: q2 FeatureTable[Frequency] object
        Qiime2 taxonomy object
    
    sample_meta: q2 Metadata Object
        Qiime2 metadata object
    
    Notes
    -----
    '''
    
    #Convert meta to q2 object
    sample_meta = q2.Metadata(meta.set_index('sample_name'))
    
    ##Create all Q2 Tables needed for Calculations##
    
    #RS210 Non-Rare
    ft_genome_rs210 = table_prep(table_genome_rs210)
    
    #WOL2 Non-Rare
    ft_genome_wol2 = table_prep(table_genome_wol2)
    
    ##Call all combinations of Beta Diversity##
    
    rs210_rpca_genome = beta_decoide(ft_genome_rs210, meta, metric, permutations, decoide_min_feature_count)
    
    wol2_rpca_genome = beta_decoide(ft_genome_wol2, meta, metric, permutations, decoide_min_feature_count)
    
    #Rarifed samples should be run 10x times and averaged
    
    #Create df to store averages across runs
    df_wol2_uwUni_genome = pd.DataFrame(columns = ['p-value', 'pseudo-F', 'Sample Size'])
    df_wol2_wUni_genome = pd.DataFrame(columns = ['p-value', 'pseudo-F', 'Sample Size'])
    df_wol2_rpca_rare_genome = pd.DataFrame(columns = ['p-value', 'pseudo-F', 'Sample Size'])
    df_rs210_rpca_rare_genome = pd.DataFrame(columns = ['p-value', 'pseudo-F', 'Sample Size'])
        
    for i in range(0, numRares):
        df_wol2_rpca_rare_genome, df_rs210_rpca_rare_genome = rare_multiplier_species(
                                                                                                                    df_wol2_rpca_rare_genome,
                                                                                                                         df_rs210_rpca_rare_genome, 
                                                                                                                         table_genome_rs210, 
                                                                                                                         table_genome_wol2,
                                                                                                                         meta, rarefaction, 
                                                                                                                         metric, decoide_min_feature_count,
                                                                                                                         permutations=999)

    print('Wol2 RPCA Rare Genome')
    display(df_wol2_rpca_rare_genome)
    print(df_wol2_rpca_rare_genome.mean())
    print()
    print('RS210 RPCA Rare Genome')
    display(df_rs210_rpca_rare_genome)
    print(df_rs210_rpca_rare_genome.mean())
    
    return(rs210_rpca_genome, wol2_rpca_genome, df_wol2_rpca_rare_genome, df_rs210_rpca_rare_genome)