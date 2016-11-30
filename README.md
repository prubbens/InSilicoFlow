# Overview

This repository contains the code used for the paper "*Flow cytometric single-cell identification of populations in
synthetic bacterial communities*". 

## Dependencies
InSilicoFlow depends on the following packages: 
- [Python 3](https://www.python.org/) (version 3.4.5)
- [numpy](http://www.numpy.org/) (version 1.11.1)
- [scikit-learn](http://scikit-learn.org/stable/) (version 0.18)
- [pandas](http://pandas.pydata.org/) (version 0.18.1)

## insilico.py
Script used to perform the first part of the analysis explained in the paper. 

## invitro.py
Script used to retrieve the composition of a synthetic bacterial community. 

## Data availability
Our data will be made available in .fcs-format on the [flowRepository](https://flowrepository.org/), 
using the identifiers FR-FCM-XXXX and FR-FCM-XXXX. It has been preprocessed following the robust digital gating strategy by [Prest et al. (2013)](http://www.sciencedirect.com/science/article/pii/S0043135413008361). 

## fcstocsv.py
The output of this preprocessing still is .fcs; therefore the data first needs to be transformed to a .csv-format; fcstocsv.py takes care of that. It makes use of [fcsparser](https://github.com/eyurtsev/fcsparser).  
