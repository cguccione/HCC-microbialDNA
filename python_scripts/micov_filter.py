'''
Functions to run micov filter
'''

#Basic imports
import pandas as pd
import biom
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

def choose_taxa(database_name):
    '''Pull taxonomy table for corresponding database'''
    
    if database_name == 'RS210':
        taxa = pd.read_csv('/Users/cguccion/Dropbox/Storage/HelpfulLabDocs/taxonomy_trees/RS210/RS210.txt', sep='\t', names=['GOTU', 'taxa'])
    elif database_name == 'WOL2':
        taxa = pd.read_csv('/Users/cguccion/Dropbox/Storage/HelpfulLabDocs/taxonomy_trees/WOL2/lineages.txt', sep='\t', names=['GOTU', 'taxa'])
    elif database_name == 'WOL':
        taxa = pd.read_csv('/Users/cguccion/Dropbox/Storage/HelpfulLabDocs/taxonomy_trees/WOL/ranks.tsv', sep='\t')
        taxa['taxa'] = taxa['kingdom'] + ' ' + taxa['phylum'] + ' ' + taxa['order'] + ' ' + taxa['family'] + ' ' + taxa['genus'] + ' ' + taxa['species']
        taxa=taxa.rename(columns={'genome':'GOTU'})
    else:
        print('Taxa not available')
        
    return(taxa)
        
def results_import(run_name, biom_name):
    '''Import all needed df from run'''
    
    biom_table = biom.load_table(biom_name)
    
    #Genome Length
    gen_len = pd.read_csv(run_name + '_subset/' + run_name + '_subset.coverage', sep='\t')

    #Covered Positions
    pos= pd.read_csv(run_name + '_subset/' + run_name + '_subset.covered_positions', sep='\t')

    #Biom Table
    biom_df = biom_table.to_dataframe(dense=True)

    #Coverages
    cov= pd.read_csv(run_name + '_subset/coverage_percentage.txt', sep='\t', names=['GOTU', 'coverage_percent'])
    
    return(gen_len, pos, biom_df, cov)

def og_micov_plot(df, genome_length, genome_id, taxa_name):
    '''Create classic Micov plot'''
    
    # Plotting
    plt.figure(figsize=(25, 8))

    # Get unique sample ids
    samples = df['sample_id'].unique()
    
    # Plotting with axes flipped
    for i, sample in enumerate(samples):
        sample_df = df[df['sample_id'] == sample]
        plt.scatter([i + 1] * len(sample_df), sample_df['start'], label=sample, s=5)

    plt.xticks(range(1, len(samples) + 1), samples)
    plt.xlim(0, len(samples) + 1)
    plt.ylim(0, genome_length+100000)
    plt.ylabel('Genome Position')
    plt.xlabel('Sample')
    plt.title(str(genome_id) + ': ' + str(taxa_name) + ', Length: ' + str(genome_length))
    plt.xticks(rotation=90)
    plt.show()

def micov_main(genome_id, gen_len, pos, biom_df, cov, taxa, pass_, fail_, bin_size, plot):
    '''Creates dfs needed to run calculations, then produces thresholds for results '''
    
    #Status on if genome_id passes or fails
    current=None

    ''' Pull Basic Info on Genome ID '''
    #Genome Name
    taxa_name, =  taxa[taxa['GOTU'] == genome_id]['taxa']

    #Print coverage
    coverage=cov[cov['GOTU'] == genome_id].iloc[0]['coverage_percent']
    
    #Genome length
    genome_length = gen_len[gen_len['genome_id'] == genome_id].iloc[0]['length']
    
    #Greater than 10% coverage -- keep regardless (Zebra paper)
    if coverage > 10:
        pass_.append(genome_id)
        current='Pass -- Coverage'
        return(pass_, fail_, None, genome_length, taxa_name)

    #Df with just hits against this genome_id
    df = pos[pos['genome_id']== genome_id]
    df = df.copy() #avoids python warnings
    
    ''' Created bin_df with all bins and if it is connected to one below it
        + Updates df which contains all hits/reads to have info about which bin(s) it hit'''
    #Create empty bin df
    bins = pd.IntervalIndex.from_breaks(range(0, genome_length + bin_size + bin_size, bin_size))
    bin_df = pd.DataFrame({
            'bin_start': [interval.left for interval in bins],
            'bin_stop': [interval.right for interval in bins],
            'sample_hits': [set() for _ in range(len(bins))],
            'read_hits': [0 for _ in range(len(bins))]})
    
    df_hit_index_list = []
    df_hit_start_stop_list = []
    
    for i, row in df.iterrows():

        #Storage variables
        hit_index_list=[]
        hit_start_stop_list=[]

        #Find the first bin in bins_df
        mask = (bin_df['bin_start'] <= row['start']) & (bin_df['bin_stop'] > row['start'])
        hit = bin_df[mask]
        
        hit_start, = hit['bin_start']
        hit_stop, = hit['bin_stop']
        hit_index, = hit.index

        #Add first bin to storage lists
        hit_index_list.append(hit_index)
        hit_start_stop_list.append([hit_start, hit_stop])

        #Update the bin list to have this sample_id
        bin_df.loc[hit_index, 'sample_hits'].add(row['sample_id'])
        #Update the bin list to have this read count
        bin_df.loc[hit_index, 'read_hits'] += 1
        
        #Check to see if this should be in other bins
        while hit_stop < row['stop']:

            hit_index+=1
            hit_start = bin_df.loc[hit_index, 'bin_start']
            hit_stop = bin_df.loc[hit_index, 'bin_stop']

            #Add current bin to storage lists
            hit_index_list.append(hit_index)
            hit_start_stop_list.append([hit_start, hit_stop])

            #Update the bin list to have this sample_id
            bin_df.loc[hit_index, 'sample_hits'].add(row['sample_id'])
            #Update the bin list to have this read count [NOTE: this is double counting reads then]
            bin_df.loc[hit_index, 'read_hits'] += 1

        df_hit_index_list.append(hit_index_list)
        df_hit_start_stop_list.append(hit_start_stop_list)

    df['hit_index_list'] = df_hit_index_list
    df['hit_start_stop_list'] = df_hit_start_stop_list

    #Remove empty rows in bin_df
    bin_df = bin_df[bin_df['sample_hits'].apply(lambda x: len(x) > 0)]

    #Add column with sample numbers
    bin_df['num_samples'] = bin_df['sample_hits'].apply(len)

    #Sort df by bin_start instead of num_samples
    bin_df = bin_df.sort_values(by='bin_start')
    bin_df = bin_df.reset_index(drop=True)

    # Check if bin_start in the current row equals bin_stop in the prior row
    bin_df['prior_bin_stop'] = bin_df['bin_stop'].shift(1)
    bin_df['is_connected'] = (bin_df['bin_start'] == bin_df['prior_bin_stop'])

    '''Uses collected info to create thresholds for pass/fail '''
    print(genome_id, taxa_name)
    print('Coverage Percentage: ' + str(coverage) + '%')

    #Less than 5 hits across BIOM table -- then remove GOTU 
    if biom_df.loc[genome_id].sum() < 5:
        fail_.append(genome_id)
        current='Fail -- less than 5 BIOM hits'
        
    #Less than 10 total reads
    elif len(df) < 30:
        fail_.append(genome_id)
        current='Fail -- less than 30 reads'
    
    if current == None:
        
        reads_top_three, one_bin_only_reads, group_bins = micov_calc(genome_id, bin_df, df, biom_df)
        
        #Less than three total bins -- remove (no way this graph looks good)
        if len(group_bins) < 4:
            fail_.append(genome_id)
            current='Fail -- less than 3 bins'

            if plot==True:
                og_micov_plot(df, genome_length, genome_id, taxa_name)
            
        #Over 75% of reads are part of a sample which only hit a single bin
        elif one_bin_only_reads > 75:
            fail_.append(genome_id)
            current='Fail -- over 75% of reads are part of a sample which only hit a single bin'

            if plot==True:
                og_micov_plot(df, genome_length, genome_id, taxa_name)

        #Over 75% of reads only hit the top 3 bins with the most samples
        elif reads_top_three > 75:
            fail_.append(genome_id)
            current='Fail - too many reads in top 3 bins'

            if plot==True:
                og_micov_plot(df, genome_length, genome_id, taxa_name)

        else:
            pass_.append(genome_id)
            current='Pass'

            if plot==True:
                og_micov_plot(df, genome_length, genome_id, taxa_name)

    print('Status:', current)
    print('----------------------------')
    return(pass_, fail_, df, genome_length, taxa_name)
    
def micov_calc(genome_id, bin_df, df, biom_df):
    '''Calculates/Reports the distrubtion of reads across groups '''
    
    ''' Determine the number of reads in each group '''
    main_list=[] #List of lists (of blocks)
    main_samples=[] #list of sets (of samples)
    main_reads=[] #List of # of reads
    current=[]
    current_samples=set()
    current_reads=0
    
    for index, row in bin_df.iterrows():
        
        if index > int(len(bin_df)-2):
            next_row_connect=False
        else:
            next_row_connect = bin_df.iloc[index+1]['is_connected']
        
        #Solo group
        if row['is_connected'] == False and next_row_connect == False:
            #Add row to main list
            main_list.append([[row['bin_start'], row['bin_stop']]])
            main_samples.append(row['sample_hits'])
            main_reads.append(row['read_hits'])
        
        #At the end of a group
        elif row['is_connected'] == True and next_row_connect == False:
            #Add row to current list
            current.append([row['bin_start'], row['bin_stop']])
            current_samples = current_samples.union(row['sample_hits'])
            current_reads+= row['read_hits']
            
            #Add group to main list
            main_list.append(current)
            main_samples.append(current_samples)
            main_reads.append(current_reads)
            current=[]
            current_samples=set()
            current_reads=0
            
        #The the middle of a group
        else:
            #Add row to current list
            current.append([row['bin_start'], row['bin_stop']])
            current_samples = current_samples.union(row['sample_hits'])
            current_reads+= row['read_hits']
        
    main_samples_counts=[]
    for i in main_samples:
        main_samples_counts.append(len(i))
    
    #Some reads fall into more than 1 bin so they are double counted -- this accounts for that
    adj_read_counts=bin_df['read_hits'].sum()
    print('Adj, read counts', adj_read_counts)
        
    '''Printing various facts now that they are calculated'''
    #print('Mains list', main_list)
    #print('Main samples',main_samples)
    print('Number of samples per group', main_samples_counts)
    print('Number of bins', len(bin_df))
    print('total # samples', len(set(df['sample_id'])))
    print('Number of total reads', len(df))
    print('Number of hits in BIOM tables', biom_df.loc[genome_id].sum())
    print()
    
    '''Calculate the percentage of samples in the bin with the most samples'''
    
    #Create df with groups of bins and samples in each 
    group_bins = pd.DataFrame({
    'bin_group': main_list,
    'samples': main_samples,
    'sample_counts': main_samples_counts,
    'read_counts': main_reads
    })
    group_bins = group_bins.sort_values(by='sample_counts', ascending=False)
    group_bins = group_bins.reset_index()
    
    '''
    #The number of samples in the bin group with the largest amount of samples
    top_one = sum(sorted(main_samples_counts, reverse=True)[:1])
    top_one_ = group_bins.loc[0, 'sample_counts']
    big_bin = (top_one/ len(set(df['sample_id']))) *100
    print('Percentage of samples in the bin group with the most samples', big_bin, '%')
    '''
    
    #The number of reads in the top 3 bin groups (sorted by # of samples)
    reads_top_three = ((group_bins.loc[:2, 'read_counts']/adj_read_counts).sum()) *100
    print('Percentage of reads in the top 3 bin groups with the largest samples', reads_top_three, '%')
    
    #The number of samples in the top 3 bin groups with largest amount of samples
    if len(group_bins) >= 4:
        #Percentage of samples in the three largest groups with the most samples
        top_three = len(group_bins.loc[0, 'samples'] | group_bins.loc[1, 'samples'] | group_bins.loc[0, 'samples'])
        three_bin = (top_three/ len(set(df['sample_id']))) *100
        
        #Percentage of samples which only exist in the three largest groups with the most samples
        top_three_rows = group_bins.iloc[:3]
        remain_rows = group_bins.iloc[3:]
        
        top_three_rows_set = set.union(*top_three_rows['samples'])
        remain_rows_set = set.union(*remain_rows['samples'])
        
        only_top_three_set = set()
        for i in top_three_rows_set:
            if i not in remain_rows_set:
                only_top_three_set.add(i)
        
        only_top_three = (len(only_top_three_set)/ len(set(df['sample_id'])))*100
        
    else:
        three_bin = None
        only_top_three = None
        
    #print('Percentage of samples in the three largest groups with the most samples', three_bin, '%')
    #print('Percentage of samples which only exist in the three largest groups with the most samples', only_top_three, '%')
    
    '''Calculate the percentage of samples which only exist in one bin group'''
    #Count the number of times each sample occurs in a bin group
    sample_per_group = defaultdict(int)
    
    # Iterate over each set in the list
    for current_group in main_samples:
        for sample in current_group:
            sample_per_group[sample] += 1    
    
    #List of samples which only hit one bin
    one_bin_samples_list = [key for key, value in sample_per_group.items() if value == 1]
    #Amount of reads which are in those samples
    one_bin_only_reads = (len(df[df['sample_id'].isin(one_bin_samples_list)])/ len(df)) *100
    print('Percentage of reads from samples who only fall into one bin group', one_bin_only_reads, '%')
            
    # Percentage of the number of samples that only hit one bin
    one_bin_sample = (sum(1 for value in sample_per_group.values() if value == 1) / len(sample_per_group)) * 100
    # Percentage of the number of samples that hit 3 bins or less
    three_bin_sample = (sum(1 for value in sample_per_group.values() if value <= 3) / len(sample_per_group)) * 100
    
    #print('Percentage of samples which only exist in one bin group', one_bin_sample, '%')
    #print('Percentage of samples which only exist in three bin groups or less', three_bin_sample, '%')
    
    return(reads_top_three, one_bin_only_reads, group_bins)    
    
def loop_micov(gen_len, pos, biom_df, cov, taxa, plot=False, bin_size=10000): #10000
    '''Loop through all items in biom_df and report results on passing/failing gotus'''
    pass_ = []
    fail_ = []

    for genome_id in list(biom_df.index):
        if genome_id in list(gen_len['genome_id']):
            #'''
            try:
                pass_, fail_, df, genome_length, taxa_name = micov_main(genome_id, gen_len, pos, biom_df, cov, taxa, pass_, fail_, plot=plot, bin_size=bin_size)
            except Exception as e:
                continue
            #'''
            #pass_, fail_, df, genome_length, taxa_name = micov_main(genome_id, gen_len, pos, biom_df, cov, taxa, pass_, fail_, plot=plot, bin_size=bin_size)
            
            print('Length of pass', len(pass_))
            print('Length of fail', len(fail_))
    
    return(pass_, fail_)