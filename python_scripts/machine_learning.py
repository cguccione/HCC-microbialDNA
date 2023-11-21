'''
Helper functions used for machine learning
and feature selection
'''

import pandas as pd
import numpy as np

from statistics import mean
import warnings
import graphviz
import pydot

from sklearn.feature_selection import chi2, VarianceThreshold, SelectKBest, f_classif, mutual_info_classif, RFE, SelectFromModel
from sklearn.model_selection import train_test_split, StratifiedKFold, LeaveOneOut, ShuffleSplit, cross_val_score, cross_val_predict
from sklearn.linear_model import RidgeCV, LassoCV, Ridge, Lasso, LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import RocCurveDisplay, auc
from sklearn import metrics, tree
import matplotlib.pyplot as plt

#WOL2 Taxonomy
wol2 = pd.read_csv('/Users/cguccion/Dropbox/Storage/HelpfulLabDocs/taxonomy_trees/WOL2/lineages.txt', sep='\t', header=None, index_col=0)
wol2_dict = wol2[1].to_dict() 

#ML Helper Functions

def prepML_convertTable2NonDuplicates(table, meta, metric):
    
    """
    The goal of this function is to convert samples 
    into correct format for ML.
    
    Since we have non-indepedant samples,
    we need to group them to insure we don't
    seperate samples from same person in test/train.
    The easiest way to do this is by removing all
    _2 (duplicates) from our dataset. We can put 
    this non-duplicates df into our desired ML
    splitting techquie. 
    
    Parameters
    ----------
    table : pandas DataFrame
        Taxonomy table with microbial counts listed as each row, and
        sample names as each coulumn
    meta : pandas DataFrame
        Metadata related to the samples names in table
    metric : string
        Column in meta file that the data should be split 
        into for ML. Ex. If the column name was 'Prognosis',
        then ML might try and determine differences between
        'Non-Progressors' and 'Progressors'. 
    
    Returns
    -------
    X : pandas DataFrame
        ML ready table
        Correctly sorted to match y. 
    y_nd : list
        ML read list with all duplicates removed
        Correctly sorted to match X_nd.
    X : pandas DatFrame
        ML ready table with everything. Taxonomy table with
        sample names listed as each row, and microbial counts 
        listed as each coulumn. Not sorted.
    """
    
    #Convert table into correct format for ML, by putting
    #samples as the rows and micorbes as columns
    X = table.T
            
    #Sort meta in ascending order to make sure it matches X
    meta = meta.sort_values('sample_name', ascending=True)
        
    #Create an X based on meta
    X_sort = X[X.index.isin(list(meta['sample_name']))]
    #Sort X in ascending order to make sure it matches
    #y and metadata
    X_sort = X_sort.sort_index(ascending=True)
        
    #Calculate corresponding y in same order
    y_sort = meta[metric].tolist()
    
    return(X_sort, y_sort, X)

def testTrainSplit_convertDF2NonDuplicates(X_nd, y_nd, X, meta, metric, False_outcome):
    
    """
    The goal of this function is to handle the 
    duplicates in our data whle spliting the data 
    into test and train sets.
    
    Since we have non-indepedant samples,
    we need to group them to insure we don't
    seperate samples from same person in test/train.
    The easiest way to do this is by removing all
    .2 (duplicates) from our dataset and putting this
    into cross_val.split(), then once we are given the
    approriate splittings, add the non-duplicates of our
    df back in.
    
    Parameters
    ----------
    X_nd : Pandas Dataframe
        Taxonomy table with sample naems listed in each column,
        and micorbial counts listed in each row. Duplicates removed.
    y_nd : list
        Corresponding outcome data for X_nd
    X : Pandas DataFrame
        Same as X_nd but with all samples.
    meta : pandas DataFrame
        Metadata related to the samples names in table
    metric : string
        Column in meta file that the data should be split 
        into for ML. Ex. If the column name was 'Prognosis',
        then ML might try and determine differences between
        'Non-Progressors' and 'Progressors'. 
    False_outcome : string
        String which keeps track of which should be labled as 
        'False' vs. 'True' when creating y list in ML
    
    Returns
    -------
    X_train : Pandas DataFrame
        ML ready array with all micorbial counts for
        each sample. 
    y_train : Pandas DataFrame
        ML ready array with all microbial counts for
        each sample. 
    X_test : Pandas DataFrame
        ML ready array with all micorbial counts for
        each sample. 
    y_test : Pandas DataFrame
        ML ready array with all microbial counts for
        each sample. 
    """
    
    #Create test/train split based on non-duplicated (nd) meta
    #Stratify insures that test/train groups have a balanced number 
    #of each metric (ex. Background/Tumor)
    X_train_nd, X_test_nd, y_train_nd, y_test_nd = train_test_split(X_nd, y_nd, test_size=0.30, stratify = y_nd)

    #Create X_train and y_train df/lists based on the non-duplicate
    #datasets, but with duplicates. 
        
    #Create a list with the duplicates of things in X_train
    X_train_list = [i + '.2' for i in list(X_train_nd.index)]
    #Add the non-duplicates to the list
    X_train_list = X_train_list + list(X_train_nd.index)
        
    #Build df based upon this list with both duplicates 
    #and non-duplicates
    X_train = X[(X.index.isin(X_train_list))]
        
    #Build y_train list based off of X_train
    y_train = []
    for i in list(X_train.index):
        y_temp, = meta[meta['sample_name'] == i][metric]
        y_train.append(y_temp)
    y_train = [False if i == False_outcome else True for i in y_train]
        
    #Create X_test and y_test df/lists based on the non-duplicate
    #datasets, but with duplicates. 

    #Create a list with the duplicates of things in X_est
    X_test_list = [i + '.2' for i in list(X_test_nd.index)]
    #Add the non-duplicates to the list
    X_test_list = X_test_list + list(X_test_nd.index)

    #Build df based upon this list with both duplicates 
    #and non-duplicates
    X_test = X[(X.index.isin(X_test_list))]
        
    #Build y_train list based off of X_test
    y_test = []
    for i in list(X_test.index):
        y_temp, = meta[meta['sample_name'] == i][metric]
        y_test.append(y_temp)
    y_test = [False if i == False_outcome else True for i in y_test]
    
    return(X_train, y_train, X_test, y_test)

def crossValSplit_convertArrays2NonDuplicates(current_X_nd_array, X, meta, metric):
    
    """
    The goal of this function is to take the arrays created
    from the cross_val.split() function in scikit learn and
    add our duplicated samples.
    
    Since we have non-indepedant samples,
    we need to group them to insure we don't
    seperate samples from same person in test/train.
    The easiest way to do this is by removing all
    .2 (duplicates) from our dataset and putting this
    into cross_val.split(), then once we are given the
    approriate splittings, add the non-duplicates of our
    df back in.
    
    Parameters
    ----------
    current_X_nd_array : numpy array
        Array of arrays with each subarray composed of:
        ['sample_name', microbe_1_count, microbe_2_count, ...].
        This is the X_nd_array[train] or X_nd_array[test] output
        from cross_val.split().
    X : pandas DatFrame
        ML ready table with everything. Taxonomy table with
        sample names listed as each row, and microbial counts 
        listed as each coulumn. Not sorted.
    meta : pandas DataFrame
        Metadata related to the samples names in table
    metric : string
        Column in meta file that the data should be split 
        into for ML. Ex. If the column name was 'Prognosis',
        then ML might try and determine differences between
        'Non-Progressors' and 'Progressors'. 
    
    Returns
    -------
    X_temp : Numpy array
        ML ready array with all micorbial counts for
        each sample. This will be used as X_train, or
        X_test depedning on inputs. 
    y_temp : Numpy array
        ML ready array with all microbial counts for
        each sample. This will be used as y_train, or
        y_test depending on inputs. 

    """
    
    #Reset and calculate X_temp, y_temp for current split
    X_temp = np.empty([0, (len(current_X_nd_array[0])-1)])
    y_temp_list = []
    
    for i in range(0, len(current_X_nd_array)):
        sample_name = current_X_nd_array[i][0] #Find correlated sample
        nd_array = current_X_nd_array[i][1:] #Remove sample name from nd array
        #Find duplicated array based on sample name
        dup_array = X.loc[[sample_name + '.2']].to_numpy() 
 
        #Add both non-duplicated and duplicated array to X_train
        X_temp = np.append(X_temp, [nd_array], axis = 0)
        X_temp = np.append(X_temp, dup_array, axis = 0)
            
        #Add both non-duplicated and duplicated outcomes to y_train
        y_current, = meta[meta['sample_name'] == sample_name][metric]
        y_temp_list.append(y_current)
        y_current, = meta[meta['sample_name'] == sample_name + '.2'][metric]
        y_temp_list.append(y_current)
    
    #Convert y_train_list to numpy array 
    y_temp = np.array(y_temp_list)
    
    return(X_temp, y_temp)

def plot_ROC(ax, tprs, mean_fpr, aucs, fn):
    
    '''
    Plotting from: https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
    '''
    
    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )
    
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(
        mean_fpr,
        tprs_lower,
        tprs_upper,
        color="grey",
        alpha=0.2,
        label=r"$\pm$ 1 std. dev.",
    )
    
    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title="Receiver operating characteristic example",
    )
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    #Saving plot as svg
    #plt.rcParams['svg.fonttype'] = 'none'
    #plt.savefig('ML_ROC_curves/' + str(fn) + '_gradientBoosting_ML_cv.svg', format = 'svg')
    
    
    plt.show()

def call_all_ML(df2run, df2run_fn, meta, metric, cross_validation_splits=3, gb_max_depth = 3, gb_n_estimators = 100, loops = 100, gb_min_samples_leaf = 5):
    
    """
    Call all ML models at once in a uniform matter.
    
    Parameters
    ----------
    df2run : list
        A list of dataframes that need to be run with ML.
    df2run_fn : list
        A list of the names of all dataframes that need to 
        be run with ML.
    meta : pandas DataFrame
        Metadata related to the samples names in table
    metric : string
        Column in meta file that the data should be split 
        into for ML. Ex. If the column name was 'Prognosis',
        then ML might try and determine differences between
        'Non-Progressors' and 'Progressors'. 
    --Optional Parameters--
    cross_validation_splits : int
        'n_splits' argument in StratifiedKFold.
        The total number of folds desired for the dataset.
        Deafult=3
    gb_max_depth : int
        'max_depth' argument in GradientBoostingClassifier
        The maximum depth limits the numebr of notdes in 
        the tree. Default=3
    gb_n_estimators : int
        'n_estimators' argument in GradientBoostingClassifier
        The number of boosting stages to preform. Default = 100
    loops: int
        The number of times we run the classic Test/TrainML model. 
        We take the average across all loops and report that score 
        at the end. Default = 100
        
    Returns
    -------
    """
    
    for i in range(0, len(df2run)):
        
        print('----------------------------------------------------')
        print('Dataframe Type:', df2run_fn[i])
        print()
        
        ##Avearge Outputs##
                   
        #Classic Test/Train
        gb_TT = gradientBoosting_ML_testTrain(df2run[i], meta, metric, 
                                             gb_max_depth=gb_max_depth, gb_n_estimators=gb_n_estimators,
                                             loops = loops, gb_min_samples_leaf=gb_min_samples_leaf)
        print('Gradient Boosting Test/Train Average:', gb_TT)
        print()
        
        #Leave one out Cross Validation
        gb_loo = gradientBoosting_ML_loo(df2run[i], meta, metric,
                                         gb_max_depth=gb_max_depth, gb_n_estimators=gb_n_estimators,
                                        gb_min_samples_leaf=gb_min_samples_leaf)
        print('Gradient Boosting Leave One Out Average:', gb_loo)
        print()
        
        ##Graph Outputs##
        
        print('Stratifed K-fold Cross Validation')
        #Stratifed K-fold Cross Validation
        gradientBoosting_ML_cv(df2run[i], meta, metric, df2run_fn[i],
                               gb_max_depth=gb_max_depth, gb_n_estimators=gb_n_estimators,
                               cross_validation_splits=cross_validation_splits,
                               gb_min_samples_leaf=gb_min_samples_leaf)
        print()
        
        print('Shuffle Split Cross Validation')
        #Shuffle Split Cross Validation
        gradientBoosting_ML_ss(df2run[i], meta, metric, df2run_fn[i],
                               gb_max_depth=gb_max_depth, gb_n_estimators=gb_n_estimators,
                               cross_validation_splits=cross_validation_splits,
                               gb_min_samples_leaf=gb_min_samples_leaf)
        
        print()
        
    return()        
    
    
#ML

def gradientBoosting_ML_testTrain(table, meta, metric, gb_max_depth=2, gb_n_estimators=100, loops = 100, gb_min_samples_leaf=5):
    
    """
    Uses Gradient Boosting ML to determine differences across listed 
    metric using the table and metadata. 
    
    Prior to ML, basic split into test and train datasets using the 
    scikit learn built in train_test_split
    
    Parameters
    ----------
    table : pandas DataFrame
        Taxonomy table with microbial counts listed as each row, and
        sample names as each coulumn
    meta : pandas DataFrame
        Metadata related to the samples names in table
    metric : string
        Column in meta file that the data should be split 
        into for ML. Ex. If the column name was 'Prognosis',
        then ML might try and determine differences between
        'Non-Progressors' and 'Progressors'. 
    --Optional Parameters ---
    gb_max_depth : int
        'max_depth' argument in GradientBoostingClassifier
        The maximum depth limits the numebr of notdes in 
        the tree. Default=3
    gb_n_estimators : int
        'n_estimators' argument in GradientBoostingClassifier
        The number of boosting stages to preform. Default = 100
    loops: int
        The number of times we run the ML model. We take the 
        average across all loops and report that score at the
        end. Default = 100
    
    Returns
    -------
    float
        Average perfomance of ML classifier across all loops.
    
    Note
    ----
    ML based on the following: https://queirozf.com/entries/visualizing-machine-learning-models-examples-with-scikit-learn-and-matplotlib
    """
    
    #Prep ML and convert table to non-duplicates
    X_nd, y_nd, X = prepML_convertTable2NonDuplicates(table, meta, metric)
    
    #Need to determine what is True/False in ML model and keep
    #this consistant throughout the rest of the code
    possible_outcomes = list(set(y_nd))
    True_outcome = possible_outcomes[0]
    False_outcome = possible_outcomes[1]
    
    #convert y_nd list to True/False
    y_nd = [False if i == False_outcome else True for i in y_nd]
    
    auc_score_list = [] #List ot average at the end of all loops
    for i in range(0, loops):
        
        #Create gradient boosting classifier model
        gbc = GradientBoostingClassifier(max_depth = gb_max_depth, n_estimators = gb_n_estimators, min_samples_leaf = gb_min_samples_leaf)
        
        ##Split data into test / train
        X_train, y_train, X_test, y_test = testTrainSplit_convertDF2NonDuplicates(X_nd, y_nd, X, meta, metric, False_outcome)
        
        #Fit the model based upon training data
        gbc.fit(X_train, y_train)
    
        #Extract predictions from model
        y_pred = gbc.predict_proba(X_test)[:,1]

        #Find false and true positive rate
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
    
        #Caclulate AUC score
        auc_score = metrics.auc(fpr, tpr)
        auc_score_list.append(auc_score)

    
    return(mean(auc_score_list))

def gradientBoosting_ML_cv(table, meta, metric, fn, gb_max_depth=2, gb_n_estimators=100, cross_validation_splits=5, gb_min_samples_leaf = 5):
    
    """
    Uses Gradient Boosting ML to determine differences across listed 
    metric using the table and metadata. 
    
    Uses cross validation to iterate multiple ML runs. Uses stratified
    K-folds cross-validator to provide folds which preserve the 
    percentage of samples for each class.
    
    Parameters
    ----------
    table : pandas DataFrame
        Taxonomy table with microbial counts listed as each row, and
        sample names as each coulumn
    meta : pandas DataFrame
        Metadata related to the samples names in table
    metric : string
        Column in meta file that the data should be split 
        into for ML. Ex. If the column name was 'Prognosis',
        then ML might try and determine differences between
        'Non-Progressors' and 'Progressors'. 
    fn : string
        Description of the file which will be used to name
        the output png/pdf
    --Optional Parameters ---
    gb_max_depth : int
        'max_depth' argument in GradientBoostingClassifier
        The maximum depth limits the numebr of notdes in 
        the tree. Default=3
    gb_n_estimators : int
        'n_estimators' argument in GradientBoostingClassifier
        The number of boosting stages to preform. Default = 100
    cross_validation_splits : int
        'n_splits' argument in StratifiedKFold.
        The total number of folds desired for the dataset.
        Deafult=5
    
    Returns
    -------
    graph
        Graph showing ROC curve with perfomance of ML classifier 
        across cross_validation_splits
    
    Note
    ----
    ML based on the following: https://queirozf.com/entries/visualizing-machine-learning-models-examples-with-scikit-learn-and-matplotlib
    Plotting based on the following: https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
    """
    
    #Create gradient boosting classifier model
    gbc = GradientBoostingClassifier(max_depth = gb_max_depth, n_estimators = gb_n_estimators, min_samples_leaf = gb_min_samples_leaf)
    
    #Convert to numpy to be compatable in ML
    X = table.T
    X_array = X.to_numpy()
    y = meta[metric].tolist()
    y_array = np.array(y)
    
    #Calculate cross validation size to ensure we don't make them too small
    print('Cross validation testing size =', len(y)/cross_validation_splits)
    #Provides train/test indices to split data in train/test sets
    cross_val = StratifiedKFold(n_splits=cross_validation_splits)
    
    #Prep for plotting
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    fig, ax = plt.subplots()
    
    '''Loop through test & train arrays caclulated by cross_val, and fix
    them by adding in duplicates. Then, preforming actual ML learning on each group.''' 
    for fold, (train, test) in enumerate(cross_val.split(X_array, y_array)):
        
        #Run ML on current split
        gbc.fit(X_array[train], y_array[train])
        
        #Calculate visulization for ROC curves
        viz = RocCurveDisplay.from_estimator(
            gbc,
            X_array[test],
            y_array[test],
            name="ROC fold {}".format(fold),
            alpha=0.3,
            lw=1,
            ax=ax,
        )
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)
    
    #Create ROC curve
    plot_ROC(ax, tprs, mean_fpr, aucs, fn)
    
    return()

def gradientBoosting_ML_loo(table, meta, metric, gb_max_depth=2, gb_n_estimators=100, gb_min_samples_leaf=5):
    
    """
    Uses Gradient Boosting ML to determine differences across listed 
    metric using the table and metadata. 
    
    Uses leave one out cross validation to iterate multiple ML runs. 
    Leave one out cross validatoin mean each sample is used once as a test set
    while the remaining samples form the training set. In our case, we still
    keep the pairs together so will test every pair. 
    
    Parameters
    ----------
    table : pandas DataFrame
        Taxonomy table with microbial counts listed as each row, and
        sample names as each coulumn
    meta : pandas DataFrame
        Metadata related to the samples names in table
    metric : string
        Column in meta file that the data should be split 
        into for ML. Ex. If the column name was 'Prognosis',
        then ML might try and determine differences between
        'Non-Progressors' and 'Progressors'. 
    --Optional Parameters ---
    gb_max_depth : int
        'max_depth' argument in GradientBoostingClassifier
        The maximum depth limits the numebr of notdes in 
        the tree. Default=3
    gb_n_estimators : int
        'n_estimators' argument in GradientBoostingClassifier
        The number of boosting stages to preform. Default = 100
    
    Returns
    -------
    float
        Average of the mean accuracy of the test set for every subgroup
    
    Note
    ----
    ML based on the following: https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.LeaveOneOut.html
    """

    #Convert to numpy to be compatable in ML
    X = table.T
    X_array = X.to_numpy()
    y = meta[metric].tolist()
    y_array = np.array(y)
    
    #Provides train/test indices to split data in train/test sets
    loo = LeaveOneOut()
    
    #Create gradient boosting classifier model
    gbc = GradientBoostingClassifier(max_depth = gb_max_depth, n_estimators = gb_n_estimators, min_samples_leaf= gb_min_samples_leaf)    
    
    #Loop through test & train arrays caclulated by LeaveOneOut, and fix
    #them by adding in duplicates. Then, preforming actual ML learning on each group.
    score_list = []
    for fold, (train, test) in enumerate(loo.split(X_array, y_array)):
        
        #Run ML on current split
        gbc.fit(X_array[train], y_array[train])
        
        #Add resulting mean accurarcy of test data
        score_list.append(gbc.score(X_array[test], y_array[test]))
        
    
    return(mean(score_list))
    

def gradientBoosting_ML_ss(table, meta, metric, fn, gb_max_depth=2, gb_n_estimators=100, cross_validation_splits=5, test_size = 0.3, gb_min_samples_leaf = 6):
    
    """
    Uses Gradient Boosting ML to determine differences across listed 
    metric using the table and metadata. 
    
    Uses cross validation to iterate multiple ML runs. Uses random 
    permutations to create test and train datasets. Does not 
    guarantee that all folds are different, but will likely happen.
    
    Parameters
    ----------
    table : pandas DataFrame
        Taxonomy table with microbial counts listed as each row, and
        sample names as each coulumn
    meta : pandas DataFrame
        Metadata related to the samples names in table
    metric : string
        Column in meta file that the data should be split 
        into for ML. Ex. If the column name was 'Prognosis',
        then ML might try and determine differences between
        'Non-Progressors' and 'Progressors'. 
    fn : string
        Description of the file which will be used to name
        the output png/pdf
    --Optional Parameters ---
    gb_max_depth : int
        'max_depth' argument in GradientBoostingClassifier
        The maximum depth limits the numebr of notdes in 
        the tree. Default=3
    gb_n_estimators : int
        'n_estimators' argument in GradientBoostingClassifier
        The number of boosting stages to preform. Default = 100
    cross_validation_splits : int
        'n_splits' argument in StratifiedKFold.
        The total number of folds desired for the dataset.
        Deafult=5
    test_size : float, <1
        This number is the percentage of our dataset that is
        used for testing. Default = 0.3
    
    Returns
    -------
    graph
        Graph showing ROC curve with perfomance of ML classifier 
        across cross_validation_splits
    
    Note
    ----
    ML based on the following: https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.ShuffleSplit.html
    Plotting based on the following: https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
    """
    
    #Create gradient boosting classifier model
    gbc = GradientBoostingClassifier(max_depth = gb_max_depth, n_estimators = gb_n_estimators, min_samples_leaf = gb_min_samples_leaf)
    
    #Convert to numpy to be compatable in ML
    X = table.T
    X_array = X.to_numpy()
    y = meta[metric].tolist()
    y_array = np.array(y)
    
    #Calculate cross validation size to ensure we don't make them too small
    print('Cross validation testing size =', (len(y) * (1- test_size))/cross_validation_splits)
    #Provides train/test indices to split data in train/test sets
    shuffle = ShuffleSplit(n_splits=cross_validation_splits, test_size=test_size)
    
    #Prep for plotting
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    fig, ax = plt.subplots()
    
    '''Loop through test & train arrays caclulated by cross_val, and fix
    them by adding in duplicates. Then, preforming actual ML learning on each group.''' 
    for fold, (train, test) in enumerate(shuffle.split(X_array, y_array)):
        
        #Run ML on current split
        gbc.fit(X_array[train], y_array[train])
        
        #Calculate visulization for ROC curves
        viz = RocCurveDisplay.from_estimator(
            gbc,
            X_array[test],
            y_array[test],
            name="ROC fold {}".format(fold),
            alpha=0.3,
            lw=1,
            ax=ax,
        )
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)
    
    #Create ROC curve
    plot_ROC(ax, tprs, mean_fpr, aucs, fn)
    
    return()

#ML Feature Selection Models : Helper Functions

def add_taxonomy_col(df, database_name):
    
    '''
    Converts gOTUs into Full Taxonmoy 
    '''
    
    if database_name == 'wol2':
        print('Converting WOLr2 Taxonomy')
        
        # Convert gOTUS into readable taxonomy
        df['gOTU'] = df['taxonomy']
        df['taxonomy'] = df['taxonomy'].map(wol2_dict)
        
    else:
        print('Database not available')
        
    return(df)

def chi2_topHits(X_train, y_train, X, k=10, print_topK = False, fn = False, taxa='None'):
    
    '''
    Calculates Univariate Chi2 testing
    
    This score can be used to select the n_features features with the
    highest values for the test chi-squared statistic from X, which must
    contain only non-negative features such as booleans or frequencies
    (e.g., term counts in document classification), relative to the classes.
    
    Recall that the chi-square test measures dependence between stochastic
    variables, so using this function "weeds out" the features that are the
    most likely to be independent of class and therefore irrelevant for
    classification.
    
    Parameters
    ----------
    X_train : Pandas DataFrame
        Df with the micorbial counts for all samples in 
        testing dataset.
    y_train : Pandas DataFrame
        Df with all the results of samples in testing data.
    X : Pandas df
        Taxonomy table with microbial counts listed as each row, and
        sample names as each coulumn
    --Optional Parameters ---
    k : int
        The number of top features we which to output.
        Default: 10
    print_topK : bool
        If the top k features should be printed to the screen.
        Defaults: False
    fn : string
        Description of the file which will be used to name
        the output tsv. If labled 'False' nothing will be saved.
        Default: False
    taxa : string
        Convert taxonomy from gOTU into acutal string to improve
        readablity. Default: None
    
    Returns
    -------
    
    Note
    ----
    Sources: https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.chi2.html#sklearn.feature_selection.chi2
    '''
    print('\nChi2 Univariate Feature Selection Significant Microbes\n')
    
    #Create chi2 model to select top k features
    sel_chi2 = SelectKBest(chi2, k=k)
    
    #Run chi2 model on trainning dataset
    X_train_chi2 = sel_chi2.fit_transform(X_train, y_train)
    
    #Loop through df and find which microbes are in top k of chi2
    if print_topK == True:
        for i in range(0, len(X.columns)):
            if sel_chi2.get_support()[i] == True:
                print(X.columns[i])

            
    #Create output of chi2 df with all p-values and chi2 values
    if fn != False:
        chi_stat, p_value = chi2(X_train, y_train)
        chi_df = pd.DataFrame(columns = ['taxonomy','Chi2', 'p-value'])
        for i in range(0, len(X.columns)):
            if str(chi_stat[i]) != 'nan':
                chi_df = chi_df.append({'taxonomy': X.columns[i], 'Chi2': chi_stat[i], 'p-value':p_value[i]}, ignore_index = True)
        
        #Convert gOTU into taxonomy
        if taxa != 'None':
            chi_df = add_taxonomy_col(chi_df, taxa)
        
        chi_df = chi_df.sort_values(by ='p-value', ascending=True)
        chi_df = chi_df.reset_index(drop = True)
        chi_df.to_csv('outputs/feature_selection/chi2/' + str(fn) + '_Univariate_Chi2.tsv', sep = '\t', index = False)

    return()

def call_all_feature_selection(df2run, df2run_fn, meta, metric, taxa='None'):
    
    """
    Call all feature selections at once in a uniform matter.
    
    Parameters
    ----------
    df2run : list
        A list of dataframes that need to be run with ML.
    df2run_fn : list
        A list of the names of all dataframes that need to 
        be run with ML.
    meta : pandas DataFrame
        Metadata related to the samples names in table
    metric : string
        Column in meta file that the data should be split 
        into for ML. Ex. If the column name was 'Prognosis',
        then ML might try and determine differences between
        'Non-Progressors' and 'Progressors'. 
    --Optional Parameters--
    taxa : string
        Convert taxonomy from gOTU into acutal string to improve
        readablity. Default: None
        
    Returns
    -------
    """
    
    for i in range(0, len(df2run)):
        
        print('----------------------------------------------------')
        print('Dataframe Type:', df2run_fn[i])
        print()
                   
        #Univariate Chi2
        univariate_feature_selection(df2run[i], meta, metric, 
                                    fn = df2run_fn[i], taxa=taxa)
        print()
        
        #Random Forest Classifier
        random_forest_feature_selection(df2run[i], meta, metric, 
                                        fn = df2run_fn[i], taxa=taxa)
        print()
        
        #Decsion Tree
        #dt = decisionTreeClassifier(df2run[i], meta, metric, 
        #                            fn = df2run_fn[i])

        #graphviz.Source(dt)
        
        
    return()        
    
#ML Feature Selection Models
def univariate_feature_selection(table, meta, metric, k=10, print_topK=False, fn = False, taxa='None'):
    
    """
    Extracts top features using univariate techniques. 
    Univariate feature selection -- Selects the best features based on statstical test
    
    Compares each feature (one at a time) to target variable, see if there is a stat 
    significant relationship. 
    
    Univariate feature selection examines each feature individually to determine the 
    strength of the relationship of the feature with the response variable. These 
    methods are simple to run and understand and are in general particularly good for 
    gaining a better understanding of data (but not necessarily for optimizing the 
    feature set for better generalization). There are lot of different options for 
    univariate selection.
    
    Really good article on which test is best: 
    https://towardsdatascience.com/mistakes-in-applying-univariate-feature-selection-methods-34c43ce8b93d
    
    More info: https://scikit-learn.org/stable/modules/feature_selection.html#univariate-feature-selection
    https://blog.datadive.net/selecting-good-features-part-i-univariate-selection/
    
    Parameters
    ----------
    table : pandas DataFrame
        Taxonomy table with microbial counts listed as each row, and
        sample names as each coulumn
    meta : pandas DataFrame
        Metadata related to the samples names in table
    metric : string
        Column in meta file that the data should be split 
        into for ML. Ex. If the column name was 'Prognosis',
        then ML might try and determine differences between
        'Non-Progressors' and 'Progressors'. 
    --Optional Parameters ---
    k : int
        The number of top features we which to output.
        Default: 10
    print_topK : bool
        If the top k features should be printed to the screen.
        Defaults: False
    fn : string
        Description of the file which will be used to name
        the output tsv. If labled 'False' nothing will be saved.
        Default: False
    taxa : string
        Convert taxonomy from gOTU into acutal string to improve
        readablity. Default: None
    
    Returns
    -------
    
    Note
    ----
    """
    
    #Prep ML and convert table to non-duplicates
    X_nd, y_nd, X = prepML_convertTable2NonDuplicates(table, meta, metric)
    
    #Need to determine what is True/False in ML model and keep
    #this consistant throughout the rest of the code
    possible_outcomes = list(set(y_nd))
    True_outcome = possible_outcomes[0]
    False_outcome = possible_outcomes[1]
    
    #convert y_nd list to True/False
    y_nd = [False if i == False_outcome else True for i in y_nd]
    
    #Split data into test / train
    X_train, y_train, X_test, y_test = testTrainSplit_convertDF2NonDuplicates(X_nd, y_nd, X, meta, metric, False_outcome)
    
    #Variance Threshold -- removes all features whose variance are zero
    sel_variance_threshold = VarianceThreshold()
    X_train_remove_variance = sel_variance_threshold.fit_transform(X_train)
    
    chi2_topHits(X_train, y_train, X, k=k, print_topK = print_topK, fn = fn, taxa=taxa)
    
    return()

def random_forest_feature_selection(table, meta, metric, fn=False, taxa='None'):
    
    """
    Feature selection from Random Forest
    
    
    Parameters
    ----------
    table : pandas DataFrame
        Taxonomy table with microbial counts listed as each row, and
        sample names as each coulumn
    meta : pandas DataFrame
        Metadata related to the samples names in table
    metric : string
        Column in meta file that the data should be split 
        into for ML. Ex. If the column name was 'Prognosis',
        then ML might try and determine differences between
        'Non-Progressors' and 'Progressors'. 
    --Optional Parameters ---
    fn : string
        Description of the file which will be used to name
        the output tsv. If labled 'False' nothing will be saved.
        Default: False
    taxa : string
        Convert taxonomy from gOTU into acutal string to improve
        readablity. Default: None
    *Can add in RandomForestClassifer Parameters*
    
    Returns
    -------
    
    Note
    ----
    Guide: https://towardsdatascience.com/feature-selection-using-python-for-classification-problem-b5f00a1c7028
    """
    print('\nSelect from Model, Random Forest Feature Selection Significant Microbes\n')
    
    #Prep ML and convert table to non-duplicates
    X_nd, y_nd, X = prepML_convertTable2NonDuplicates(table, meta, metric)
    
    #Need to determine what is True/False in ML model and keep
    #this consistant throughout the rest of the code
    possible_outcomes = list(set(y_nd))
    True_outcome = possible_outcomes[0]
    False_outcome = possible_outcomes[1]
    
    #convert y_nd list to True/False
    y_nd = [False if i == False_outcome else True for i in y_nd]
    
    #Split data into test / train
    X_train, y_train, X_test, y_test = testTrainSplit_convertDF2NonDuplicates(X_nd, y_nd, X, meta, metric, False_outcome)
    
    #Create Random Forest Tree for feature selection, fit data based upon this
    model_tree = RandomForestClassifier(n_estimators=50, min_samples_split =5, max_depth = 10)
    model_tree.fit(X_train, y_train)
    
    #Create tree from the model we fit above
    # Features whose importance is greater or equal to the threshold are kept while the others are discarded.
    sel_model_tree = SelectFromModel(estimator=model_tree, prefit=True, threshold='mean')  
    
    #Use our new tree on our X_train dataset
    X_train_sfm_tree = sel_model_tree.transform(X_train)
    
    #Loop through random forest tree, extract each bacteria and its tree importance, than sort
    sel_model_tree_df = pd.DataFrame(columns = ['Tree_Importance', 'taxonomy'])
    for i in range(0, len(X.columns)):
        if sel_model_tree.get_support()[i] ==True:
            sel_model_tree_df = sel_model_tree_df.append({'Tree_Importance' : model_tree.feature_importances_[i], 'taxonomy' : X.columns[i]}, ignore_index=True)
    sel_model_tree_df = sel_model_tree_df.sort_values(by=['Tree_Importance'], ascending = False)
    
    #Convert gOTU into taxonomy
    if taxa != 'None':
        sel_model_tree_df = add_taxonomy_col(sel_model_tree_df, taxa)
    
    #Save output if requested
    if fn != False:
        sel_model_tree_df.to_csv('outputs/feature_selection/random_forest/' + str(fn) + '_topFeaturesRandomForest.tsv' , sep='\t', index = False)
    
    return()

def decisionTreeClassifier(table, meta, metric, fn, min_samples_split=5):
    
    #Excellent resource on how to read decsion trees: https://towardsdatascience.com/scikit-learn-decision-trees-explained-803f3812290d 
    
    #Create list with metric to perfrom
    #ML on and convert list into True/False
    #instead of tumor/Background, 
    #and create list of all microbes
    X = table.T
    y = meta[metric].tolist()
    #y = [False if i == y[0] else True for i in y]
    
    #Build and run ML model
    dtr = tree.DecisionTreeClassifier(min_samples_split = min_samples_split)
    dtr_f = dtr.fit(X, y) #Full dataset
    
    #Creating tree file in dot
    dot_data = tree.export_graphviz(dtr_f, out_file='decision_trees/' + fn + '.dot', filled=True, feature_names = X.columns, class_names = pd.Series(y).drop_duplicates().tolist())
    
    #Display the tree in jupyter
    with open('decision_trees/' + fn + '.dot') as f:
        dot_graph = f.read()
        
    #Save the tree as a pdf
    (graph,) = pydot.graph_from_dot_file('decision_trees/' + fn + '.dot')
    graph.write_pdf('decision_trees/' + fn + '.pdf')
    
    #Save the tree as png
    graph.write_png('decision_trees/' + fn + '.png')
    
    ### Use Decsion Tree as actual ML 
    #Nodes in Decsion trees are determined by the feature that produces a split with purest subsets.

    #Run the Cross Validation with ENTIRE dataset... cross validation breaks it up
    #**NOTE.. this is done after fitting, compared to with test/train it's split before
    #** You can change the scoring  varaible, but I didn't notice much, scoring= 'r2'
    scores = cross_val_score(dtr_f, X, y, cv= 10) #CV stands for the number of folds with default of 5

    #Compare predictions made from  cross validation with actual results
    #**Note that more data can be checked here because we can use the whole dataset to check
    #**Compared to before where we could just use the test dataset to check 
    predictions = cross_val_predict(dtr_f, X, y, cv= 10)

    #Print the results
    print ('Score:' + str(scores)) #Returns cross validation scores
    print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
    #**Also, gives the mean score and the 95% confidence interval of score 
    #**because there are 10 scores instead of 1
    
    return(dot_graph)