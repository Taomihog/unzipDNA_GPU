# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 20:49:16 2023

@author: taomihog
"""
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_palette("husl")

import pandas as pd
import numpy as np
import os

def linear_interpolate(x, y, x0):
    y0 = np.interp(x0, x, y, left = -1.0, right = -1.0)
    return y0

plot = True
fmax = 20
fmin = 10
min_bp_trace_to_save = 3000
max_bp_trace_to_save = 5000

# pyResult = pd.read_csv("standard.txt",delimiter='\t'

fig,ax = plt.subplots(figsize=(16, 50))

folder_path = os.getcwd() + "/parsed_unzip_curves"
prefix = "out"
y_shift = 0
n_traces = 0

files_with_prefix = [f for f in os.listdir(folder_path)]

for file in files_with_prefix:
    cppResult = pd.read_csv(folder_path + '/' + file)
    x = cppResult['total_extension_nm'].values
    y = cppResult['force_avg_pN'].values
    
    # remove the high force region
    il = 0
    ih = len(x) - 1
    cnt = 0
    while(True):
        cnt += 1
        imid = (int)(il + (ih - il)//2)
        if y[imid] < fmax:
            il = imid
        else:
            ih = imid
        if ih - il == 1:
            break;
    
    x = x[:ih]
    y = y[:ih]
    while(ih - 2 >= 0 and y[ih - 1] - y[ih - 2] > 0):
        ih -= 1
    if ih + 20 < len(x):
        x = x[:ih + 20]
        y = y[:ih + 20]
    
    # remove the low force region
    il = 0
    ih = len(x) - 1
    cnt = 0
    while(True):
        cnt += 1
        imid = (int)(il + (ih - il)//2)
        if y[imid] < fmin:
            il = imid
        else:
            ih = imid
        if ih - il == 1:
            break;
    
    x = x[il:]
    y = y[il:]
    il = 0
    while(il + 1 < len(y) and y[il] - y[il + 1] < 0):
        il += 1
    if il >= 20:
        x = x[il - 20:]
        y = y[il - 20:]
    
    # interpolate the integer extension start at zero
    xmax = x[-1] - x[0]
    x0 = [x for x in range(xmax)]
    y0 = linear_interpolate(x, y, x0)
    if xmax < min_bp_trace_to_save or xmax > max_bp_trace_to_save:
        # do not save if the sequence is too short, for test purposes
        continue
    
    
    
    name = file[len(prefix):-len(".csv")]
    print(xmax, name) 
    ax.plot(x,y + y_shift,label = name)
    plt.text(800, 12 + y_shift, name, fontsize = 20)
    y_shift += 6
    n_traces += 1


# plt.ylim(9.8,18.2)
# plt.ylim(0,30)

print("total traces:", n_traces)

plt.xlim(500,5000)
plt.xlabel("extension (nm)", fontsize = 20)
plt.ylabel("force (pN", fontsize = 20)
plt.yticks([])
plt.tick_params(axis='both', which='major', labelsize = 20)
# plt.legend(loc='best')
plt.savefig("result.png", dpi = 600)
plt.show()

