# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 15:16:37 2022

@author: mda153
"""

import pandas as pd
import all_paths as ap

import os
import o3seespy as o3
import numpy as np
import matplotlib.pyplot as plt
import sfsimodels as sm
import loc_o3fibre
import loc_fd_analysis
from bwplot import cbox
import all_paths as ap

name = 'diff_fd_length'

version = 'v1'

version_list = ['v1', 'v2']


out_folder = f'{ap.PM_STUDY_RESULTS}{name}/{version}/'

lfs = [10, 20, 30]

for lf in lfs:
    prefix = f'LF{lf}'.replace('.', 'p')
    
    columns = np.arange(0, 21, 1)
    max_list = []
    
    #outs = ['nx', 'node_disps', 'forces', 'mom', 'shear', 'ndisps', 'spring_forces'] 
    item = "ndisps"
    ndisp_all = pd.read_csv(out_folder + f'{prefix}_{item}.txt', sep=' ', header=None)
    item = 'nx'
    loss_support_length = pd.read_csv(out_folder + f'{prefix}_{item}.txt', sep=' ', header=None)
    item = 'mom'
    bend_moment = pd.read_csv(out_folder + f'{prefix}_{item}.txt', sep=' ', header=None)
    for column in columns: 
        data_list = ndisp_all.iloc[:, column].to_numpy()
        max_list.append(max(abs(data_list)))

print('max list:', max_list)
print('overall max:', max(max_list))



# print(pd.read_csv('//file/Usersm$/mda153/Home/My Documents/2022/Course Work/ENCI682/Parametric Studydiff_fd_length/v3/LF10_ndisps.txt', sep=' ', header=None))