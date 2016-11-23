# Overview

This repository contains all the code used for the paper "*Flow cytometric single-cell identification of populations in
synthetic bacterial communities*". 

## Dependencies
InSilicoFlow depends on the following packages: 
- [Python 3](https://www.python.org/) (version )
- [numpy](http://www.numpy.org/) (version )
- [scikit-learn](http://scikit-learn.org/stable/) (version )
- [pandas](http://pandas.pydata.org/) (version )


## fcstocsv.py
The analysis starts from the output of a flow cytometer, after which it has been preprocessed following the robust digital gating strategy by [Prest etl al. (2013)](http://www.sciencedirect.com/science/article/pii/S0043135413008361). As the output of this preprocessing still is .fcs, the data first needs to be transformed to a .csv-format; fcstocsv.py takes care of that. 

## insilico.py
Script used to perform the first part of the analysis explained in the paper. 

## invitro.py
Script used to retrieve the composition of a synthetic bacterial community. 
