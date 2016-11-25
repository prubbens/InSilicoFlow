# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 11:13:05 2016

@author: prubbens
"""

##############################################################################
###Import packages############################################################
##############################################################################

''' Imported packages for python ''' 
import numpy as np
import pandas as pd
import pylab as plt
import time
import warnings
import time
from os import listdir

''' Imported packages from scikit-learn '''
from sklearn.cross_validation import train_test_split 
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier

''' Define plotting style, ignore future warnings, start stopwatch '''
plt.style.use('ggplot')
warnings.simplefilter(action = "ignore", category = FutureWarning) #Do not display futurewarnings
start_time = time.time() #Start stopwatch to determine runtime

##############################################################################
###Read-in metadata and define global variables###############################
##############################################################################


'''
Read-in metadata: 
This data is our expected outcome concerning the in vitro communities in the abundance gradient 
'''
targetabundances = pd.read_excel('targetabundances.xlsx', index_col = 'File')

''' Global variables'''
Nax = 5000
Ncomm = 10000
invitro_combination = 2 #0: low, 1: medium, 2: high

##############################################################################
###Functions##################################################################
##############################################################################


''' Return path of microbial community of interest ''' 
''' Input: '''
'comb == 0: Pseudomonas putida -- Pseudomonas fluorescens (initial low performance)'
'comb == 1: Agrobacter rhizogenes -- Janthinobacterium sp. B3 (initial medium performance)'
'comb == 2: Shewanella oneidensis -- Micrococcus luteus (initial high performance)'
def get_path(comb): 
    if comb == 0: 
        path_insilico = 'PC_filtered_1_11/' 
        path_invitro = '1_11_filtered/'
        nrep_insilico = 3
        print(path_invitro)
    elif comb == 1: 
        path_insilico = 'PC_filtered_rerun_16_25/'
        path_invitro = '16_25_rerun_filtered/'
        nrep_insilico = 4
    else: 
        path_insilico = 'PC_filtered_3_17/'
        path_invitro = '3_17_filtered/'
        nrep_insilico = 3
    nrep_invitro = 3
    return path_insilico, nrep_insilico, path_invitro, nrep_invitro

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
    
''' Return in silico community containing two bacterial populations '''
'datalist: list of filenames containing individual bacterial populations'
'n_subsample: number of cells to sample per bacterial population'
'path: path to directory'
'nrep: number of replicates'
def get_insilico_comm(datalist, n_subsample, path, nrep): 
    df0 = pd.DataFrame()
    df1 = pd.DataFrame()
    for i in np.arange(0, nrep):        
        df_singlespecies = pd.read_csv(path + datalist[i], index_col = 0)
        df_singlespecies['species'] = 0
        df0 = pd.concat([df0,df_singlespecies], axis = 0, ignore_index = True)
    df0 = df0.sample(n_subsample, random_state = 903, replace = False)
    for j in np.arange(nrep, 2*nrep): 
        df_singlespecies = pd.read_csv(path + datalist[j], index_col = 0)
        df_singlespecies['species'] = 1
        df1 = pd.concat([df1,df_singlespecies], axis = 0, ignore_index = True)  
    df1 = df1.sample(n_subsample, random_state = 1503, replace = False)
    df = pd.concat([df0, df1], axis = 0, ignore_index = True)
    return df

''' Return in silico community containing two bacterial populations in varying abundances '''
'datalist: list of filenames containing individual bacterial populations'
'n_sample: number of cells to sample per bacterial population'
'abun: relative abundance (between 0 and 1) of first bacterial population'
'path: path to directory'
'nrep: number of replicates'
def get_insilico_comm_abun(datalist, n_sample, abun, path, nrep): 
    df0 = pd.DataFrame()
    df1 = pd.DataFrame()
    for i in np.arange(0, nrep):        
        df_singlespecies = pd.read_csv(path + datalist[i], index_col = 0)
        df_singlespecies['species'] = 0
        df0 = pd.concat([df0,df_singlespecies], axis = 0, ignore_index = True)
    df0 = df0.sample(int(n_sample*abun), random_state = 27, replace = False)
    for j in np.arange(nrep, 2*nrep): 
        df_singlespecies = pd.read_csv(path + datalist[j], index_col = 0)
        df_singlespecies['species'] = 1
        df1 = pd.concat([df1,df_singlespecies], axis = 0, ignore_index = True)  
    df1 = df1.sample(int(n_sample*(1.-abun)), random_state = 633, replace = False)
    df = pd.concat([df0, df1], axis = 0, ignore_index = True)
    return df
    
''' Sample cells of synthetic bacterial community '''
'idx: 0 or 1 (first or second bacterial population)'
'datalist: list of filenames containing synthetic communities'
'n_sample: number of cells to sample per invitro community'
'path: path to directory'
'nrep: number of replicates'
def get_invitro_comm(idx, datalist, n_sample, path, nrep): 
    df = pd.DataFrame()
    for i in np.arange(int(idx*nrep), int((idx+1)*nrep)): 
        df_rep = pd.read_csv(path + datalist[i], index_col = 0)
        df = pd.concat([df, df_rep], axis = 0, ignore_index = True)
    if(df.shape[0] > n_sample): 
        df = df.sample(n_sample, random_state = 5495, replace = False)
    return df
    
''' Calculate relative abundance p and alpha diversity parameters D1 and D2 ''' 
'cluster: array of (predicted) cell labels'
'n_clust: number of different clusters'
def calc_D1_D2(cluster, n_clust):
    cluster = pd.DataFrame(cluster)
    cluster.columns = ['clust']
    pp = 0
    p_d1 = 0.     
    p_d2 = 0.
    min_clust = np.amin(cluster.clust)
    max_clust = np.amax(cluster.clust)
    for i in np.arange(min_clust, max_clust + 1): 
        p = np.float64(cluster.loc[cluster['clust'] == i].shape[0]/cluster.shape[0])
        p_d1 += p*np.log(p)
        p_d2 += p**2.
        if pp == 0.: 
            pp = p
    return pp, np.exp(-1.*p_d1), 1./p_d2 

''' Split in silico community into a training and validation/test set '''    
'Input: dataframe'
def get_train_test(df): 
    features = list(df.columns)
    art_x_train, art_x_test, art_y_train, art_y_test = train_test_split(df[features[0:12]], df['species'], test_size = 0.3, random_state = 588)
    return art_x_train, art_x_test, art_y_train, art_y_test
    
''' Train Linear Discriminant analysis on training set and evaluate on test set '''
'x_train: dataframe training set'
'y_train: labels of training set'
'x_test: dataframe test set'
def perform_lda(x_train, x_test, y_train): 
    lda = LinearDiscriminantAnalysis(solver = 'lsqr')
    lda.fit(x_train, y_train)
    return lda.predict(x_test)
    
''' Perform LDA with dimensionality reduction '''
def perform_lda_dimred(x_train, x_test, y_train): 
    lda = LinearDiscriminantAnalysis(solver = 'svd')
    lda.fit(x_train, y_train)
    return lda.transform(x_test)
    
''' Plot feature importances from Random Forest classifier '''
'features: list of features' 
'feature importances: RF.feature_importances_'
def plot_feature_importances(features, feature_importances): 
    df = pd.DataFrame(feature_importances, index = features)
    df.sort(axis = 1, ascending = False, inplace = True)
    df.columns = ['feature_importances']
    pos = np.arange(0, len(features)) + 0.5
    plt.figure(figsize=(20,12))
    plt.barh(pos, df.feature_importances, color = 'darkorange', align = 'center')    
    plt.yticks(pos, df.index)
    plt.xlabel('Importance')
    plt.title('Feature Importances')
    plt.axis([0,0.25,0,12])
    #plt.show()
    plt.savefig('RF_featureimportances_2Species_3.png')
    
''' Train using Linear Discriminant Analysis on training set and return classifier '''
'x_train: dataframe training set'
'y_train: labels of training set'
def return_lda_class(x_train, y_train):
    lda = LinearDiscriminantAnalysis(solver = 'svd')
    lda.fit(x_train, y_train)
    return lda

''' Train Random Forest classifier on training set and return classifier '''
'x_train: dataframe training set'
'y_train: labels of training set'
def return_RF_class(x_train, y_train):
    rf = RandomForestClassifier(n_estimators = 200, criterion = 'gini', random_state = 6)
    rf.fit(x_train, y_train)
    return rf
    
''' Train Linear Discriminant analysis on training set and evaluate (labels) on test set '''
'x_train: dataframe training set'
'y_train: labels of training set'
'x_test: dataframe test set'
def perform_RF(x_train, x_test, y_train): 
    rf = RandomForestClassifier(n_estimators = 200, criterion = 'gini', random_state=3)
    rf.fit(x_train, y_train)
    #plot_feature_importances(features, rf.feature_importances_)
    return rf.predict(x_test)

''' Train Linear Discriminant analysis on training set and evaluate (probabilities) on test set '''
'x_train: dataframe training set'
'y_train: labels of training set'
'x_test: dataframe test set'
def perform_RF_scores(x_train, x_test, y_train): 
    rf = RandomForestClassifier(n_estimators = 200, criterion='gini', random_state=3)
    rf.fit(x_train, y_train)
    #std_fi = np.std([tree.feature_importances_ for tree in rf.estimators_], axis=0)
    #plot_feature_importances(features, rf.feature_importances_, std_fi)
    return rf.predict_proba(x_test)[:,1]

''' Create in silico abundance gradient and analyze it using classifier trained on initial in silico community '''
'datalist: list of filenames containing individual bacterial populations'
'Nax: number of cells to sample per population in in silico community'
'Nabun: total number of cells to sample for in silico communities making up an abundance gradient'
'path_insilico: path to directory containing in silico communities'
'nrep_insilico: number of replictates'
def perform_insilico_analysis(datalist_insilico, Nax, Nabun, path_insilico, nrep_insilico): 
    df_insilico = get_insilico_comm(datalist_insilico, Nax, path_insilico, nrep_insilico)
    features = get_features(df_insilico)
    features = get_features(df_insilico)
    clf = return_RF_class(df_insilico[features], df_insilico.species)
    percentages = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
    noc = len(percentages)
    p = np.zeros(noc)
    D1 = np.zeros(noc)
    D2 = np.zeros(noc)
    dummy = 0
    for pct in percentages: 
        df_insilico_abun = get_insilico_comm_abun(datalist_insilico, Nabun, pct, path_insilico, nrep_insilico)
        pred = clf.predict(df_insilico_abun[features])
        p[dummy], D1[dummy], D2[dummy] = calc_D1_D2(pred, 2)
        dfresult = pd.DataFrame(percentages, columns = ['Theoretical Abundances'])
        dummy += 1
    dfresult['D1'] = D1
    dfresult['D2'] = D2
    dfresult['p0'] = p
    dfresult.sort(columns='Theoretical Abundances', inplace=True)
    dfresult.to_csv('dfinsilico_abun.csv')
    return dfresult

''' Create in silico abundance gradient and analyze it using classifier trained on initial in silico community '''
'datalist_insilico: list of filenames containing individual bacterial populations'
'datalist_invitro: list of filenames containing synthetic bacterial communities in varying abundances'
'Nax: number of cells to sample per population in in silico community'
'Ninvitro: total number of cells to sample for in vitro communities making up an abundance gradient'
'path_insilico: path to directory containing in silico communities'
'path_invitro: path to directory containing in silico communities'
'nrep_insilico: number of replictates of individual bacterial populations'
'nrep_invitro: number of replictates of in vitro communities'  
'noc: number of communities making up an abundance gradient (13 in the paper)'
'targetabundances: metadata containing outcome (in vitro created) abundances'
def perform_invitro_analysis(datalist_insilico, datalist_invitro, Nax, Ninvitro, path_insilico, path_invitro, nrep_insilico, nrep_invitro, noc, targetabundances): 
    df_insilico = get_insilico_comm(datalist_insilico, Nax, path_insilico, nrep_insilico)
    features = get_features(df_insilico)
    clf = return_RF_class(df_insilico[features], df_insilico.species)
    theor_abundances = [10,1,20,30,40,50,5,60,70,80,90,95,99]
    p = np.zeros(noc)
    D1 = np.zeros(noc)
    D2 = np.zeros(noc)
    N = np.zeros(noc)
    for i in np.arange(0,noc): 
        df_invitro = get_invitro_comm(i, datalist_invitro, Ninvitro, path_invitro, nrep_invitro)
        pred = clf.predict(df_invitro[features])
        p[i], D1[i], D2[i] = calc_D1_D2(pred, 2)
        N[i] = df_invitro.shape[0]
    dfresult = pd.DataFrame(theor_abundances, columns = ['Theoretical Abundances'])
    dfresult['p0'] = p
    dfresult['D1'] = D1
    dfresult['D2'] = D2
    dfresult['N'] = N
    dfresult.sort(columns='Theoretical Abundances', inplace=True)
    dfresult['Target abundances'] = targetabundances
    dfresult.to_csv('dfinvitro.csv')
    return dfresult

##############################################################################
###Call functions#############################################################
##############################################################################


path_insilico, nrep_insilico, path_invitro, nrep_invitro = get_path(invitro_combination)
datalist_insilico = sorted(listdir(path_insilico))
datalist_invitro = sorted(listdir(path_invitro))
noc = int(len(datalist_invitro)/nrep_invitro)
df_result_insilico = perform_insilico_analysis(datalist_insilico, Nax, Ncomm, path_insilico, nrep_insilico)
df_result_invitro = perform_invitro_analysis(datalist_insilico, datalist_invitro, Nax, Ncomm, path_insilico, path_invitro, nrep_insilico, nrep_invitro, noc, targetabundances.iloc[:,2*invitro_combination].values)