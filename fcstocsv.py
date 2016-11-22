# -*- coding: utf-8 -*-

import fcsparser
import pandas as pd
import os 

'''Import data files'''
for file in os.listdir("Directory"):
    if file.endswith(".fcs"):
        dummystring = str('Directory/' + file)
        meta, data = fcsparser.parse(dummystring, meta_data_only=False, reformat_meta=False)
        file = file[:-4]
        pd_data = pd.DataFrame(data)
        dummystring2 = str("Directory/" + file + ".csv")
        pd_data.to_csv(dummystring2)
        print(file)
