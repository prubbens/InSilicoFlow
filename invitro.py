# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 11:13:05 2016

@author: prubbens
"""

#import sys
import numpy as np
import pandas as pd
import pylab as plt
from sklearn.cross_validation import train_test_split 
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.grid_search import GridSearchCV
#from sklearn import preprocessing
plt.style.use('ggplot')
from os import listdir

''' Read-in metadata '''
targetabundances = pd.read_excel('targetabundances.xlsx', index_col = 'File')

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

def get_features(df): 
    features = list(df.columns)
    if(len(features) == 15): 
        features.remove('species')
    features.remove('Time')
    features.remove('Width')
    return features    
    
''' Return Artificially Created Mock '''
''' 
Get in-silico mock from axenic cultures
'''    
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
    
def get_insilico_comm_abun(datalist, n_subsample, abun, path, nrep): 
    df0 = pd.DataFrame()
    df1 = pd.DataFrame()
    for i in np.arange(0, nrep):        
        df_singlespecies = pd.read_csv(path + datalist[i], index_col = 0)
        df_singlespecies['species'] = 0
        df0 = pd.concat([df0,df_singlespecies], axis = 0, ignore_index = True)
    df0 = df0.sample(int(n_subsample*abun), random_state = 27, replace = False)
    for j in np.arange(nrep, 2*nrep): 
        df_singlespecies = pd.read_csv(path + datalist[j], index_col = 0)
        df_singlespecies['species'] = 1
        df1 = pd.concat([df1,df_singlespecies], axis = 0, ignore_index = True)  
    df1 = df1.sample(int(n_subsample*(1.-abun)), random_state = 633, replace = False)
    df = pd.concat([df0, df1], axis = 0, ignore_index = True)
    return df
    
''' Return Biologically Created Mock '''
    
''' Pool replicates before sampling synthetic community '''
def get_invitro_comm(idx, datalist, n_sample, path, nrep): 
    df = pd.DataFrame()
    for i in np.arange(int(idx*nrep), int((idx+1)*nrep)): 
        df_rep = pd.read_csv(path + datalist[i], index_col = 0)
        df = pd.concat([df, df_rep], axis = 0, ignore_index = True)
    if(df.shape[0] > n_sample): 
        df = df.sample(n_sample, random_state = 5495, replace = False)
    return df
    
''' Calculate relative abundance p and diversity parameters D1 and D2 '''
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

''' Return binned data '''
''' Split Artificial Mock into a training and validation/test set '''    
def get_train_test(art_mock): 
    features = list(art_mock.columns)
    art_x_train, art_x_test, art_y_train, art_y_test = train_test_split(art_mock[features[0:12]], art_mock['species'], test_size = 0.3, random_state = 588)
    return art_x_train, art_x_test, art_y_train, art_y_test
    
''' Perform LDA '''
def perform_lda(x_train, x_test, y_train): 
    lda = LinearDiscriminantAnalysis(solver = 'lsqr')
    lda.fit(x_train, y_train)
    return lda.predict(x_test)
    
''' Perform LDA with dimensionality reduction '''
def perform_lda_dimred(x_train, x_test, y_train): 
    lda = LinearDiscriminantAnalysis(solver = 'svd')
    lda.fit(x_train, y_train)
    return lda.transform(x_test)
    
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
    
def tune_RF(x_train, y_train): 
    features = list(x_train.columns)
    param_grid = [{'max_features': np.arange(1,len(features))}, {'criterion': ['gini', 'entropy']}]
    clf_rf = GridSearchCV(RandomForestClassifier(n_estimators = 200), param_grid, cv = 10)
    clf_rf.fit(x_train, y_train)
    grid_scores = clf_rf.grid_scores_
    best_params = clf_rf.best_params_    
    return grid_scores, best_params

def return_lda_class(x_train, y_train):
    lda = LinearDiscriminantAnalysis(solver = 'svd')
    lda.fit(x_train, y_train)
    return lda
    
def return_RF_class(x_train, y_train):
    rf = RandomForestClassifier(n_estimators = 200, criterion = 'gini', random_state = 6)
    rf.fit(x_train, y_train)
    return rf
    
def perform_RF(x_train, x_test, y_train): 
    #features = list(x_train.columns)
    rf = RandomForestClassifier(n_estimators = 200, criterion = 'gini', random_state=3)
    rf.fit(x_train, y_train)
    #plot_feature_importances(features, rf.feature_importances_)
    return rf.predict(x_test), rf.predict_proba(x_test)[:,1]
    
def perform_RF_scores(x_train, x_test, y_train): 
    #features = list(x_train.columns)
    rf = RandomForestClassifier(n_estimators = 200, criterion='gini', random_state=3)
    rf.fit(x_train, y_train)
    #std_fi = np.std([tree.feature_importances_ for tree in rf.estimators_], axis=0)
    #plot_feature_importances(features, rf.feature_importances_, std_fi)
    return rf.predict_proba(x_test)[:,1]

def perform_insilico_analysis(datalist_insilico, Nax, path_insilico, nrep_insilico): 
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
        df_insilico_abun = get_insilico_comm_abun(datalist_insilico, Nax, pct, path_insilico, nrep_insilico)
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
    dfresult['D1'] = D1
    dfresult['D2'] = D2
    dfresult['N'] = N
    dfresult['p0'] = p
    dfresult.sort(columns='Theoretical Abundances', inplace=True)
    dfresult['Target abundances'] = targetabundances
    dfresult.to_csv('dfinvitro.csv')
    return dfresult


Nax = 5000
Ncomm = 10000
invitro_combination = 2 #0: low, 1: medium, 2: high
path_insilico, nrep_insilico, path_invitro, nrep_invitro = get_path(invitro_combination)
datalist_insilico = sorted(listdir(path_insilico))
datalist_invitro = sorted(listdir(path_invitro))
noc = int(len(datalist_invitro)/nrep_invitro)
df_result_insilico = perform_insilico_analysis(datalist_insilico, Ncomm, path_insilico, nrep_insilico)
df_result_invitro = perform_invitro_analysis(datalist_insilico, datalist_invitro, Nax, Ncomm, path_insilico, path_invitro, nrep_insilico, nrep_invitro, noc, targetabundances.iloc[:,2*invitro_combination].values)