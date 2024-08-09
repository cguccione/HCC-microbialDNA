'''
Helper functions used to calculate basic
alpha and beta diversity
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
from qiime2.plugins.diversity.pipelines import alpha
from scipy.stats import mannwhitneyu

def alpha_rare_curve(table, max_depth, meta):
    
    q2_meta = Metadata(meta) 
    
    alph_vis = alpha_rarefaction(table = table, max_depth = max_depth, metadata = q2_meta)
    
    return(alph_vis)

def beta_decoide(ft, meta, metric, permutations, return_biplot=False):
    
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
        
    return_biplot: bool (False)
        Returns the biplot instead 
        of the permaonva ouputs 
    
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
    
    #Save ft for command line processing
    ft.save('deicode_processing.qza')

    command = [
    'qiime',
    'deicode',
    'rpca',
    '--i-table',
    'deicode_processing.qza',
    '--p-min-sample-count',
    '10',
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
    
    if return_biplot == True:
        return(rpca_biplot)
                                                   
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


def rare_multiplier_fig(df_rs210_rpca_rare_genome, table_genome_rs210, meta, rarefaction, metric, decoide_min_feature_count, permutations=999):
    
    #Convert meta to q2 object
    sample_meta = q2.Metadata(meta.set_index('sample_name'))
    
    #Rare Tables
    ft_genome_rs210_rare = table_prep(table_genome_rs210, rarefaction=rarefaction)
    
    #Rarefied RPCA
    rs210_rpca_rare_genome = beta_decoide(ft_genome_rs210_rare, meta, metric, permutations, decoide_min_feature_count)
    df_rs210_rpca_rare_genome = rare_helper(rs210_rpca_rare_genome, df_rs210_rpca_rare_genome)
    
    return(df_rs210_rpca_rare_genome)

def beta_figures(table_genome_rs210,  meta, rarefaction, metric, permutations=999, numRares=10, decoide_min_feature_count=10, return_biplot=False):
    
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
        
    return_biplot: bool (False)
        Returns the distance matrix instead 
        of the permaonva ouputs in RPCA
    
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
    
    #Rarifed samples should be run 10x times and averaged
    
    #Create df to store averages across runs
    df_rs210_rpca_rare_genome = pd.DataFrame(columns = ['p-value', 'pseudo-F', 'Sample Size'])
        
    for i in range(0, numRares):
        df_rs210_rpca_rare_genome = rare_multiplier_fig(df_rs210_rpca_rare_genome, table_genome_rs210, meta, 
                                                        rarefaction, metric, decoide_min_feature_count, permutations=999)

    print('RS210 RPCA Rare Genome')
    display(df_rs210_rpca_rare_genome)
    print(df_rs210_rpca_rare_genome.mean())
    
    return(df_rs210_rpca_rare_genome)
    

def alpha_rare_multiplier_fig(df_rs210_alpha_rare_genome, table_genome_rs210, meta, rarefaction, metric):
    
    #Convert meta to q2 object
    sample_meta = q2.Metadata(meta.set_index('sample_name'))
    
    #Metadata
    meta = meta.set_index('sample_name', drop = True)
    
    #Rare Tables
    ft_genome_rs210_rare = table_prep(table_genome_rs210, rarefaction=rarefaction)
    
    #Alpha diverstiy, Shannon
    alpha_result_s = alpha(table=ft_genome_rs210_rare, metric='shannon')
    alpha_diversity_s = alpha_result_s.alpha_diversity
    alpha_series_s = alpha_diversity_s.view(pd.Series)
    
    df = pd.merge(meta, alpha_series_s, right_index = True, left_index = True)
    
    metric_list = list(set(df[metric]))
    
    #Shannon p-values
    list_0 = list(df.loc[df[metric] ==  metric_list[0]]['shannon_entropy'])
    list_1 = list(df.loc[df[metric] ==  metric_list[1]]['shannon_entropy'])
    
    #Shannon p-values
    _, pnorm = mannwhitneyu(list_0, list_1)
    
    #Add to df
    df_rs210_alpha_rare_genome.loc[len(df_rs210_alpha_rare_genome)] = [pnorm]
    
    return(df_rs210_alpha_rare_genome)

def alpha_figures(table_genome_rs210, meta, rarefaction, metric, numRares=10):
    
    '''Same as beta_figures but with alpha'''
    
    #Convert meta to q2 object
    sample_meta = q2.Metadata(meta.set_index('sample_name'))
    
    #Create df to store averages across runs
    df_rs210_alpha_rare_genome = pd.DataFrame(columns = ['p-value'])
        
    for i in range(0, numRares):
        df_rs210_alpha_rare_genome = alpha_rare_multiplier_fig(df_rs210_alpha_rare_genome, table_genome_rs210, meta, rarefaction, metric)

    print('RS210 Alpha Rare Genome')
    display(df_rs210_alpha_rare_genome)
    print(df_rs210_alpha_rare_genome.mean())
    
    return(df_rs210_alpha_rare_genome)