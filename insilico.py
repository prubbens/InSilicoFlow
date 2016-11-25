# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 17:37:37 2016

@author: prubbens
"""
##############################################################################
###Import packages############################################################
##############################################################################


''' Imported packages for python ''' 
import numpy as np
import pandas as pd
import scipy as sc
import pylab as plt
import itertools
import warnings
import time
from os import listdir

''' Imported packages from scikit-learn '''
from sklearn.cross_validation import train_test_split 
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier

''' Define plotting style, ignore future warnings, start stopwatch '''
plt.style.use('ggplot') #Plotting style for matplotlib
warnings.simplefilter(action = "ignore", category = FutureWarning) #Do not display futurewarnings
start_time = time.time() #Start stopwatch to determine runtime

##############################################################################
###Read-in metadata###########################################################
##############################################################################

''' Put species names in an excel sheet to later on annotate in silico communities '''
path = 'DataInsilicoComm/'
datalist_singlespecies = sorted(listdir(path)) 
list_species = pd.read_excel('listspecies.xlsx', index_col = 0, header = 0 )
list_species = list_species.index.tolist()

##############################################################################
####Functions#################################################################
##############################################################################

'''Return number of combinations '''
'Input: Number of populations (N), community complexity (k)'
'Output: Number of possible combinations'
def get_number_of_combinations(N, k): 
    return sc.misc.comb(N, k)
    
''' Filter out features you don't want to use for your classifier '''
'Input: pandas dataframe'
'Output: list of features'
def get_features(df): 
    features = list(df.columns)
    if(len(features) == 15): 
        features.remove('species')
    features.remove('Time')
    features.remove('Width')
    return features
    
''' If there are technical replicates, pool them, and subsample a number of cells '''    
'Input: index of bacterial population (idx), number of cells (n_cell), list of filenames (datalist), number of replicates (n_rep)' 
'Output: dataframe (df)'
def get_subsample_ax_cul_pool(idx, datalist, n_cell, n_rep, species_id): 
    df = pd.DataFrame()    
    for i in np.arange(0, n_rep): 
        df_subsample = pd.read_csv(path + datalist[int(n_rep*idx) + i], index_col = 0)  
        df =  pd.concat([df, df_subsample], axis = 0, ignore_index = True)  
    df_sampled = df.sample(n_cell, random_state = 567, replace = False)
    df_sampled['species'] = species_id
    return df_sampled
    
''' Split in silico community into a training and validation/test set '''
'Input: dataframe'    
def get_train_test(df): 
    features = get_features(df)
    x_train, x_test, y_train, y_test = train_test_split(df[features], df['species'], test_size = 0.3, random_state = 26)
    return x_train, x_test, y_train, y_test
    
''' Calculate ROC '''    
def return_roc_auc_score(y_true, y_scores):
    return metrics.roc_auc_score(y_true, y_scores, average = None)
    
''' Calculate accuracy '''
def return_accuracy_score(y_true, y_pred): 
    return metrics.accuracy_score(y_true, y_pred)

''' Train Linear Discriminant analysis on training set and evaluate on test set '''
'x_train: dataframe training set'
'y_train: labels of training set'
'x_test: dataframe test set'
def perform_lda(x_train, x_test, y_train): 
    lda = LinearDiscriminantAnalysis(solver = 'lsqr')
    lda.fit(x_train, y_train)
    return lda.predict(x_test), lda.predict_proba(x_test)[:,1]
    
''' Train Random Forest classifier on training set and evaluate on test set '''
'x_train: dataframe training set'
'y_train: labels of training set'
'x_test: dataframe test set'    
def perform_RF(x_train, x_test, y_train): 
    rf = RandomForestClassifier(n_estimators = 200, criterion = 'gini')
    rf.fit(x_train, y_train)
    return rf.predict(x_test), rf.predict_proba(x_test)[:,1]
    
''' Check performance on a held-out test set for all possible combinations when S=2 '''
'Input: list of datafile names (filenames), list of species names (list_species), number of cells to sample (n_sample), number of replicates (n_rep)'
'Output: dataframe (df) containing peformances (accuracy and AUC) along with community names'  
def perform_RF_2species_allcomb(filenames, list_species, n_sample, n_rep): 
    N = 20
    k = 2
    noc = int(get_number_of_combinations(N, k))
    df = pd.DataFrame()
    scores = np.zeros(noc)
    auc = np.zeros(noc)
    listspecies_i = list()
    listspecies_j = list()
    teller = 0
    for i in np.arange(0,N): 
        for j in np.arange(i+1,N): 
            print(teller)
            df_ax_cul1 = get_subsample_ax_cul_pool(i, filenames, n_sample, n_rep, 0)            
            df_ax_cul2 = get_subsample_ax_cul_pool(j, filenames, n_sample, n_rep, 1)
            df_comm = pd.concat([df_ax_cul1, df_ax_cul2], axis = 0, ignore_index = True)
            x_train, x_test, y_train, y_test = get_train_test(df_comm)
            y_pred, y_pred_proba = perform_RF(x_train, x_test, y_train)
            auc[teller] = return_roc_auc_score(y_test, y_pred_proba)
            scores[teller] = return_accuracy_score(y_test, y_pred)
            listspecies_i.append(list_species[i])
            listspecies_j.append(list_species[j])
            teller += 1
    df['species_i'] = listspecies_i
    df['species_j'] = listspecies_j
    df['accuracy'] = scores
    df['AUC'] = auc
    return df
    
''' Check performance on a held-out test set for 150 random combinations for S to choose '''
'Input: list of datafile names (filenames), species richness (S), number of cells to sample (n_sample), number of replicates (n_rep)'
'Output: dataframe (df) containing peformances (accuracy) along with the indices of every bacterial population constituting the in silico community'     
def get_perf_art_mock_nspecies(filenames, S, n_sample, n_rep): 
    np.random.seed(32)
    n_comm = 150
    N = 20
    perf = np.zeros(n_comm)
    poss_comb = list(itertools.combinations(np.arange(0, N), S))
    n_poss_comb = len(poss_comb)
    idx_poss_comb = np.arange(0, n_poss_comb)
    subset_idx_poss_comb = np.random.choice(idx_poss_comb, n_comm, replace = False)
    teller = 0
    for idx in subset_idx_poss_comb: 
        df = pd.DataFrame()
        comp_mock = poss_comb[idx]
        for ax in comp_mock:
            df_ax_cul = get_subsample_ax_cul_pool(ax, filenames, n_sample, n_rep, ax)
            df = pd.concat([df, df_ax_cul], axis = 0, ignore_index = True)
        x_train, x_test, y_train, y_test = get_train_test(df)
        y_pred, y_pred_proba = perform_RF(x_train, x_test, y_train)
        perf[teller] = return_accuracy_score(y_test, y_pred)
        teller += 1
    mocks = [poss_comb[i] for i in subset_idx_poss_comb]
    df_meta = pd.DataFrame()
    df_meta['composition'] = mocks
    df_meta['performance'] = perf
    #list_ax_cul = np.random.choice(20, k, replace = False)
    #get_ax_cul(idx, n_subsample, datalist)
    return df_meta
    
df = perform_RF_2species_allcomb(datalist_singlespecies, list_species, 5000, 2)
#df = get_perf_art_mock_nspecies(datalist_singlespecies, 5, 5000, 2)    

df.to_csv('resultsRF_repA_repIF.csv')

print("--- %s seconds ---" % (time.time() - start_time)) 
