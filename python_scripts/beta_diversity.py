'''
Helper functions used to calculate basic
beta diversity
'''

import pandas as pd

import qiime2 as q2
from qiime2 import Artifact, Metadata
from qiime2.plugins.feature_table.methods import rarefy
from qiime2.plugins import diversity
from qiime2.plugins.diversity_lib.methods import weighted_unifrac, unweighted_unifrac

#WOL2 Taxonomy
wol2_taxonomy = q2.Artifact.import_data('Phylogeny[Rooted]', '/Users/cguccion/Dropbox/Storage/HelpfulLabDocs/taxonomy_trees/WOL2/tree.nwk')

'''
Helper functions used to calculate basic
beta diversity
'''

import pandas as pd

import qiime2 as q2
from qiime2 import Artifact, Metadata
from qiime2.plugins.feature_table.methods import rarefy
from qiime2.plugins import diversity
from qiime2.plugins.diversity_lib.methods import weighted_unifrac, unweighted_unifrac

#WOL2 Taxonomy
wol2_taxonomy = q2.Artifact.import_data('Phylogeny[Rooted]', '/Users/cguccion/Dropbox/Storage/HelpfulLabDocs/taxonomy_trees/WOL2/tree.nwk')

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

def beta_decoide(ft, sample_meta, metric, permutations):
    
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
    beta_result_d: q2 Visulization object
        Qiime2 visulaization object of 
        decoide
    
    Notes
    -----
    '''
    
    #Save ft for command line processing
    ft.save('deicode_processing.qza')
    
    #Run deicode from command line using qza created above
    !qiime deicode rpca --i-table deicode_processing.qza --o-biplot deicode_biplot.qza --o-distance-matrix deicode_distance_test.qza
    
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

def all_beta(table_genome_rs210, table_genome_wol2, meta, rarefaction, metric, permutations=999):
    
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
    
    #WOL2 Rare
    ft_genome_wol2_rare = table_prep(table_genome_wol2, rarefaction=rarefaction)
    
    
    ##Call all combinations of Beta Diversity##
    
    rs210_rpca_genome = beta_decoide(ft_genome_rs210, sample_meta, metric, permutations)
    
    wol2_uwUni_genome, wol2_wUni_genome = beta_unifraq(ft_genome_wol2_rare, sample_meta, metric, permutations)
    
    wol2_rpca_genome = beta_decoide(ft_genome_wol2, sample_meta, metric, permutations)
    
    return(rs210_rpca_genome, wol2_uwUni_genome, wol2_wUni_genome, wol2_rpca_genome)