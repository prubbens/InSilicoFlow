# Overview

This repository contains the code and data that is used to automatically identify the taxonomy of bacterial single-cells based on flow cytometry. If you find the research and/or data useful, please consider citing: 

Rubbens, P., Props, R., Boon, N., Waegeman, W. (2017). [Flow Cytometric Single-Cell Identification of Populations in Synthetic Bacterial Communities](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169754).  *PLOS ONE* **12**(1):e0169754.
## Dependencies
InSilicoFlow depends on the following packages: 
- [Python 3](https://www.python.org/) (version 3.4.5)
- [numpy](http://www.numpy.org/) (version 1.11.1)
- [scikit-learn](http://scikit-learn.org/stable/) (version 0.18)
- [pandas](http://pandas.pydata.org/) (version 0.18.1)

## insilico.py
Script used to perform the first part of the analysis explained in the paper. It mainly performs the following steps: 
- Choose two or m bacterial populations of interest. 
- Aggregate data coming from replicates. 
- Aggregate data coming from the populations of interest (creating the so-called *in silico community*). 
- Use 70% of the community to train a classifier (e.g., LDA or Random Forests). 
- Use the other 30% to evaluate the performance of a classifier (expressed in for example the accuracy, AUC, ...). 

## invitro.py
Script used to retrieve the composition of a synthetic bacterial community. It mainly performs the following steps: 
- Create in silico community, representing the synthetic community of interest. 
- Train classifier on full community. 
- Evaluate performance by retrieving community composition for various in silico/in vitro communities. 

## Data availability
Our data can be found in .csv-format on this repo. Additionally, our data is made available in .fcs-format on the [flowRepository](https://flowrepository.org/), 
using the identifiers FR-FCM-ZZSH (axenic cultures) and FR-FCM-ZZSG. It has been preprocessed following the robust digital gating strategy by [Prest et al. (2013)](http://www.sciencedirect.com/science/article/pii/S0043135413008361). 

## fcstocsv.py
If one wants to start from the raw .fcs-files, it can be transformed to .csv-format using fcstocsv.py; therefore the data first needs to be transformed to a .csv-format; fcstocsv.py takes care of that. It makes use of [fcsparser](https://github.com/eyurtsev/fcsparser).  
