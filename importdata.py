# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import fcsparser
import pandas as pd
import os 

'''Import data files'''
for file in os.listdir("Growth Curve/filtered/FCS/Cells"):
    if file.endswith(".fcs"):
        dummystring = str('Growth Curve/filtered/FCS/Cells/' + file)
        meta, data = fcsparser.parse(dummystring, meta_data_only=False, reformat_meta=False)
        file = file[:-4]
        pd_data = pd.DataFrame(data)
        dummystring2 = str("Growth Curve/filtered/CSV/Cells/" + file + ".csv")
        pd_data.to_csv(dummystring2)
        print(file)