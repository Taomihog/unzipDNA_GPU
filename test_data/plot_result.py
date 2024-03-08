# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 20:49:16 2023

@author: taomihog
"""
import seaborn as sns
import matplotlib.pyplot as plt
# Set the color palette
sns.set_palette("husl")  # You can use different palette names like "husl", "deep", "colorblind", etc.

import pandas as pd
import numpy as np
import time
import os

# pyResult = pd.read_csv("standard.txt",delimiter='\t')

fig,ax = plt.subplots(figsize=(16, 4))

folder_path = os.getcwd()
prefix = "out"

files_with_prefix = [f for f in os.listdir(folder_path) if f.startswith(prefix)]

# ax.plot(pyResult['# zMean'],pyResult['FMean'], "k", linewidth = 2, label = "python code")
# total_extension_nm, DNA_extension_avg_nm, force_avg_pN, force_sd_pN, average_bp_unzipped_avg, bp_unzipped_sd
for file in files_with_prefix:
    print(file) 
    cppResult = pd.read_csv(file)
    # print(cppResult.columns)
    ax.plot(cppResult['total_extension_nm'],cppResult['force_avg_pN'], '--', 
            label = "GPU result, gene = "+file[len(prefix):-len(".csv")])


plt.ylim(9.8,18.2)
# plt.ylim(0,30)

plt.xlim(650,5500)
plt.xlabel("extension (nm)")
plt.ylabel("force (pN")
plt.legend(loc='best')
plt.savefig("result.svg", dpi = 600)
plt.show()

