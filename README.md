# Overview

This repository contains all the code used for the paper "*Flow cytometric single-cell identification of populations in
synthetic bacterial communities*". 

## fcstocsv.py
The analysis starts from the output of a flow cytometer, after which it has been preprocessed following the robust digital gating strategy by [Prest etl al. (2013)](http://www.sciencedirect.com/science/article/pii/S0043135413008361). As the output of this preprocessing still is .fcs, the data first needs to be transformed to a .csv-format; fcstocsv.py takes care of that. 
