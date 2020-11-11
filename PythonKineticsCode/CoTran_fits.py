#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 11:41:44 2020

@author: reyer
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 09:49:51 2019

@author: reyer
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:04:02 2019

@author: reyer
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 11:44:17 2018

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

WalkerFolder1 = "/Users/reyer/Data/Python_Simulations/8520/old/ptsG/four/"


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
burn_step_d = 1000
num_walkers = 50
num_files = 1
posterior_check = -30000

k_on_data = []
k_off_data = []
kx_s_p_data = []
b_e_data = []
lnf_data = []
kx_wt_data = []
kx_rne_data = []
b_ms_data = []
c_data = []
d_data = []

k_on_list = []
k_off_list = []
kx_s_p_list = []
b_e_list = []
kx_wt_list = []
kx_rne_list = []
b_ms_list = []
c_list = []
d_list = []

dataFile = open(WalkerFolder1 + "walkers.csv", "r")
csvReader = csv.reader(dataFile)
walkers = np.zeros((num_walkers, num_steps, 9))
counter = 0
for row in csvReader:
    first = int(counter / num_steps)
    second = counter % num_steps
    walkers[first][second] = [float(i) for i in row]
    counter+=1
likelihoodDataFile = open(WalkerFolder1 + "likelihood.csv", "r")
csvReaderLikelihoods = csv.reader(likelihoodDataFile)
likelihoods = np.zeros((num_walkers, num_steps))
posteriorDataFile = open(WalkerFolder1 + "posterior.csv", "r")
csvReaderposterior = csv.reader(posteriorDataFile)
posterior = np.zeros((num_walkers, num_steps))
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
walkersShaped = walkers.reshape((-1, 9))
k_on_ml,k_off_ml,kx_s_p_ml,b_e_ml,b_ms_ml,kx_wt_ml,kx_rne_ml,c_ml,d_ml =walkersShaped[:][np.argmax(likelihoodsShaped)]

counter = 0
for row in csvReaderposterior:
    posterior[:, counter] = [float(i) for i in row]
    counter+=1
posteriorShaped = np.array(posterior.reshape((-1)))
posteriorLike = len(posteriorShaped)
for i in range(posteriorLike):
    if (posteriorShaped[i] < -500):
        posteriorShaped[i] = -500

k_on_p = 0
k_off_p = 0
kx_s_p_p = 0
b_e_p = 0
kx_wt_p = 0
kx_rne_p = 0
b_ms_p = 0   
c_p = 0
d_p = 0     

max_count = 0
while max_count < 5:
    
    k_on,k_off,kx_s_p,b_e,b_ms,kx_wt,kx_rne,c,d = walkersShaped[:][np.argmax(posteriorShaped)]
    walkersShaped = np.delete(walkersShaped,np.argmax(posteriorShaped),0)
    posteriorShaped = np.delete(posteriorShaped,np.argmax(posteriorShaped))
    
    
    if k_on == k_on_p and k_off == k_off_p and c == c_p and d == d_p and kx_s_p == kx_s_p_p and b_e == b_e_p and kx_wt == kx_wt_p and kx_rne == kx_rne_p and b_ms == b_ms_p:
        continue
    else:
        k_on_p = k_on
        k_off_p = k_off
        kx_s_p_p = kx_s_p
        b_e_p = b_e
        kx_wt_p = kx_wt
        kx_rne_p = kx_rne
        b_ms_p = b_ms
        c_p = c
        d_p = d
        
        k_on_list.append(k_on) 
        k_off_list.append(k_off)
        kx_s_p_list.append(kx_s_p)
        b_e_list.append(b_e)
        kx_wt_list.append(kx_wt)
        kx_rne_list.append(kx_rne)
        b_ms_list.append(b_ms)
        c_list.append(c)
        d_list.append(d)
        max_count += 1
            
        
    
k_on_max = np.mean(k_on_list)
k_off_max = np.mean(k_off_list)
kx_s_p_max = np.mean(kx_s_p_list)
b_e_max = np.mean(b_e_list)
kx_wt_max = np.mean(kx_wt_list)
kx_rne_max = np.mean(kx_rne_list)
b_ms_wt_max = np.mean(b_ms_list)
c_max = np.mean(c_list)
d_max = np.mean(d_list)
    
   

step_count = 0

for j in range(num_walkers):
    for i in range(burn_step_k_on,num_steps):
        if (data1_values[i + step_count][0] < 1E-12 or data1_values[i + step_count][0] > 300):
            continue
        if posterior[j,i] < posterior_check:
            continue
        k_on_data.append(data1_values[i + step_count][0])
        
    step_count += num_steps  

step_count = 0    
for j in range(num_walkers):
    for i in range(burn_step_k_off,num_steps):
        if (data1_values[i + step_count][1] < 1E-12 or data1_values[i + step_count][1] > 300):
            continue
        if posterior[j,i] < posterior_check:
            continue
        k_off_data.append(data1_values[i + step_count][1])
        
    step_count += num_steps     

     
step_count = 0     
for j in range(num_walkers):
    for i in range(burn_step_kx_s_p,num_steps):  
        if (data1_values[i + step_count][2] > 2):
            continue
        if posterior[j,i] < posterior_check:
            continue
        kx_s_p_data.append(data1_values[i + step_count][2])
    step_count += num_steps

step_count = 0
for j in range(num_walkers):
    for i in range(burn_step_b_e,num_steps):
        if (data1_values[i + step_count][3] < 1E-8 or data1_values[i + step_count][3] > 2):
            continue
        if posterior[j,i] < posterior_check:
            continue
        b_e_data.append(data1_values[i + step_count][3])
    step_count += num_steps

    

    
step_count = 0
for j in range(num_walkers):
    for i in range(burn_step_kx_wt,num_steps):
        if (data1_values[i + step_count][5] < .0001):
            continue
        if posterior[j,i] < posterior_check:
            continue
        kx_wt_data.append(data1_values[i + step_count][5])
    step_count += num_steps

step_count = 0
for j in range(num_walkers):
    for i in range(burn_step_kx_rne,num_steps):
        if data1_values[i + step_count][6] < .0001:
            continue
        if posterior[j,i] < posterior_check:
            continue
        kx_rne_data.append(data1_values[i + step_count][6])
        
    step_count += num_steps
    

step_count = 0     
for j in range(num_walkers):
    for i in range(burn_step_b_ms,num_steps):  
#        if (data1_values[i + step_count][1] > 1):
#            continue
        if posterior[j,i] < posterior_check:
            continue
        b_ms_data.append(data1_values[i + step_count][4])
    step_count += num_steps    
        
step_count = 0       
for j in range(num_walkers):
    for i in range(burn_step_c,num_steps):  
        if data1_values[i + step_count][7] < 0:
            continue
        if posterior[j,i] < posterior_check:
            continue
        c_data.append(data1_values[i + step_count][7])
    step_count += num_steps 
     
step_count = 0       
for j in range(num_walkers):
    for i in range(burn_step_d,num_steps):  
        if data1_values[i + step_count][8] < 0:
            continue
        if posterior[j,i] < posterior_check:
            continue
        d_data.append(data1_values[i + step_count][8])
    step_count += num_steps 
        

# Fit a normal distribution to the data:

fit_alpha, fit_loc, fit_beta = gamma.fit(k_on_data)
#fit_alpha, fit_loc,  = norm.fit(kA_data)

# Plot the histogram.
plt.hist(k_on_data, normed = True, bins=25, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = gamma.pdf(x, fit_alpha, fit_loc, fit_beta)
plt.plot(x, p, 'k', linewidth=2)
title = "k_on : mu = %.3E,  std = %.3E, MLE = %.3E" % (gamma.mean(fit_alpha, fit_loc, fit_beta), gamma.std(fit_alpha, fit_loc, fit_beta), k_on_max)
plt.xscale('log')
plt.title(title)
plt.savefig(WalkerFolder1 + "k_on_fit.png")
plt.close()

fit_alpha, fit_loc, fit_beta = gamma.fit(c_data)
#fit_alpha, fit_loc,  = norm.fit(kA_data)

# Plot the histogram.
plt.hist(c_data, normed = True, bins=25, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = gamma.pdf(x, fit_alpha, fit_loc, fit_beta)
plt.plot(x, p, 'k', linewidth=2)
title = "c : mu = %.3E,  std = %.3E, MLE = %.3E" % (gamma.mean(fit_alpha, fit_loc, fit_beta), gamma.std(fit_alpha, fit_loc, fit_beta), c_max)
plt.xscale('log')
plt.title(title)
plt.savefig(WalkerFolder1 + "c_fit.png")
plt.close()

fit_alpha, fit_loc, fit_beta = gamma.fit(d_data)
#fit_alpha, fit_loc,  = norm.fit(kA_data)

# Plot the histogram.
plt.hist(d_data, normed = True, bins=25, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = gamma.pdf(x, fit_alpha, fit_loc, fit_beta)
plt.plot(x, p, 'k', linewidth=2)
title = "d : mu = %.3E,  std = %.3E, MLE = %.3E" % (gamma.mean(fit_alpha, fit_loc, fit_beta), gamma.std(fit_alpha, fit_loc, fit_beta), d_max)
plt.xscale('log')
plt.title(title)
plt.savefig(WalkerFolder1 + "d_fit.png")
plt.close()

fit_alpha, fit_loc, fit_beta = gamma.fit(k_off_data)
#fit_alpha, fit_loc,  = norm.fit(kA_data)

# Plot the histogram.
plt.hist(k_off_data, normed = True, bins=25, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = gamma.pdf(x, fit_alpha, fit_loc, fit_beta)
plt.plot(x, p, 'k', linewidth=2)
title = "k_off : mu = %.3E,  std = %.3E, MLE = %.3E" % (gamma.mean(fit_alpha, fit_loc, fit_beta), gamma.std(fit_alpha, fit_loc, fit_beta), k_off_max)
plt.xscale('log')
plt.title(title)
plt.savefig(WalkerFolder1 + "k_off_fit.png")
plt.close()


fit_alpha, fit_loc, fit_beta = gamma.fit(kx_s_p_data)

# Plot the histogram.
plt.hist(kx_s_p_data, normed = True, bins=20, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = gamma.pdf(x, fit_alpha, fit_loc, fit_beta)
plt.plot(x, p, 'k', linewidth=2)
title = "kx_s_p : mu = %.3f,  std = %.3f , MLE = %.3f" % (gamma.mean(fit_alpha, fit_loc, fit_beta), gamma.std(fit_alpha, fit_loc,fit_beta), kx_s_p_max)
plt.title(title)
plt.savefig(WalkerFolder1 + "kx_s_p_fit.png")
plt.close()


fit_alpha, fit_loc,fit_beta = gamma.fit(b_e_data)

# Plot the histogram.
plt.hist(b_e_data, normed = True, bins=25, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = gamma.pdf(x, fit_alpha,  fit_loc,fit_beta)
plt.plot(x, p, 'k', linewidth=2)
plt.xscale('log')
title = "b_e : mu = %.3E,  std = %.3E, MLE %.3E" % (gamma.mean(fit_alpha, fit_loc,fit_beta), gamma.std(fit_alpha, fit_loc,fit_beta), b_e_max)
plt.title(title)
plt.savefig(WalkerFolder1 + "b_e_fit.png")
plt.close()


# Plot the histogram.

kx_wt_mu, kx_wt_std = norm.fit(kx_wt_data)

plt.hist(kx_wt_data, normed = True, bins=25, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, kx_wt_mu, kx_wt_std)
plt.plot(x, p, 'k', linewidth=2)
title = "kx_wt : mu = %.3E,  std = %.3E, MLE = %.3E" % (kx_wt_mu, kx_wt_std, kx_wt_max)
plt.title(title)
plt.savefig(WalkerFolder1 + "kx_wt_fit.png")
plt.close()

kx_rne_mu, kx_rne_std = norm.fit(kx_rne_data)

# Plot the histogram.
plt.hist(kx_rne_data, normed = True, bins=25, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, kx_rne_mu, kx_rne_std)
plt.plot(x, p, 'k', linewidth=2)
title = "kx_rne : mu = %.3E,  std = %.3E, MLE = %.3E" % (kx_rne_mu, kx_rne_std, kx_rne_max)
plt.title(title)
plt.savefig(WalkerFolder1 + "kx_rne_fit.png")
plt.close()


fit_alpha, fit_loc = norm.fit(b_ms_data)

# Plot the histogram.
plt.hist(b_ms_data, normed = True, bins=20, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, fit_alpha, fit_loc)
plt.plot(x, p, 'k', linewidth=2)
title = "b_ms_wt : mu = %.6f,  std = %.6f , MLE = %.6f" % (norm.mean(fit_alpha, fit_loc), norm.std(fit_alpha, fit_loc), b_ms_wt_max)
plt.title(title)
plt.savefig(WalkerFolder1 + "b_ms_wt_fit.png")
plt.close()


