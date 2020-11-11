#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 18:58:24 2020

@author: reyer
"""

import numpy as np
from scipy.stats import norm
from scipy.stats import gamma
from scipy.stats import powernorm
from scipy.stats import lognorm
from scipy.stats import exponnorm
from scipy.stats import expon
from scipy.stats import loggamma
from scipy.stats import poisson
import matplotlib.pyplot as plt
import csv
import pandas as pd

WalkerFolder1 = "/Users/reyer/Data/Python_Simulations/71020/post/sodB130/six/"


#ka = args[0] #GUESS
#kd = args[1] #GUESS
#kx_s_p = args[2] #GUESS
#k_init_wt = args[3]
#k_init_rne = args[4]
#b_e = args[5]
#lnf = args[6]np
#kx_wt = args[7]
#kx_rne = args[8]

dataFile1 = open(WalkerFolder1 + "walkers.csv", "r")

data1 = pd.read_csv(dataFile1, header = None)

data1_values = data1.values

num_steps = 10000
burn_step_k_on = 1000
burn_step_k_off = 1000
burn_step_kx_s_p = 1000
burn_step_b_e = 1000
burn_step_kx_wt = 1000
burn_step_kx_rne = 1000
burn_step_b_ms = 1000
burn_step_c = 1000
num_walkers = 50
num_files = 1
posterior_check = -1000

k_on_data = []
k_off_data = []
kx_s_p_data = []
b_e_data = []
lnf_data = []
kx_wt_data = []
kx_rne_data = []
b_ms_data = []

k_on_list = []
k_off_list = []
kx_s_p_list = []
b_e_list = []
kx_wt_list = []
kx_rne_list = []
b_ms_list = []

dataFile = open(WalkerFolder1 + "walkers.csv", "r")
csvReader = csv.reader(dataFile)
walkers = np.zeros((num_walkers, num_steps, 7))
counter = 0
for row in csvReader:
    first = int(counter / num_steps)
    second = counter % num_steps
    walkers[first][second] = [float(i) for i in row]
    counter+=1
likelihoodDataFile = open(WalkerFolder1 + "likelihood.csv", "r")
csvReaderLikelihoods = csv.reader(likelihoodDataFile)
likelihoods = np.zeros((num_walkers, num_steps))

counter = 0
for row in csvReaderLikelihoods:
    likelihoods[:, counter] = [float(i) for i in row]
    counter+=1
likelihoodsShaped = np.array(likelihoods.reshape((-1)))
likelihoodsShaped2 = np.array(likelihoods.reshape((-1)))
lengthLike = len(likelihoodsShaped)
for i in range(lengthLike):
    if (likelihoodsShaped[i] < -1000):
        likelihoodsShaped[i] = -1000
walkersShaped = walkers.reshape((-1, 7))
k_on_ml,k_off_ml,kx_s_p_ml,b_e_ml,b_ms_ml,kx_wt_ml,kx_rne_ml =walkersShaped[:][np.argmax(likelihoodsShaped)]

bic = 7*np.log(42)-2*np.max(likelihoodsShaped)
print(bic)


