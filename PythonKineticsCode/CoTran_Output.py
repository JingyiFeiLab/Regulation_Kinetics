#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 09:47:48 2019

@author: reyer
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:02:21 2019

@author: reyer
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 11:02:21 2018

@author: reyer
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import math
from scipy.integrate import odeint
from multiprocessing import Pool
from scipy import integrate, stats
from scipy.odr import *
import contextlib
import os
import sys
from scipy.optimize import fmin, minimize, basinhopping
import scipy.optimize as op
import csv

WalkerFolder = "/Users/reyer/Data/Python_Simulations/8320/sodB130/four/"

num_walkers = 50
num_steps = 10000

k_on_list = []
k_off_list = []
kx_s_p_list = []
b_e_list = []
b_s_wt_list = []
b_s_rne_list = []
k_init_wt_list = []
k_init_rne_list = []
kx_wt_list = []
kx_rne_list = []
a_s_wt_list = []
a_s_rne_list = []
b_m_list = []
b_ms_list = []
c_list = []
d_list = []


# k_on_max,k_off_max,kx_s_p_max,b_e_max,e_max,b_s_wt_max,b_s_rne_max,k_init_wt_co_max,k_init_rne_co_max,k_init_wt_max,k_init_rne_max,kx_wt_max,kx_rne_max,a_s_wt_max,a_s_rne_max,b_m_max,b_ms_wt_max

def postProcess():
    dataFile = open(WalkerFolder + "walkers.csv", "r")
    csvReader = csv.reader(dataFile)
    walkers = np.zeros((num_walkers, num_steps, 9))
    counter = 0
    for row in csvReader:
        first = int(counter / num_steps)
        second = counter % num_steps
        walkers[first][second] = [float(i) for i in row]
        counter+=1
    for i in range(num_walkers):
        plt.plot(list(range(0, num_steps)), walkers[i, :, 0])
        plt.title("k_on values ")
    plt.yscale('log')
    plt.savefig(WalkerFolder + "k_on.png")
    plt.close()
    for i in range(num_walkers):
        plt.plot(list(range(0, num_steps)), walkers[i, :, 1])
        plt.title("k_off values ")
    plt.yscale('log')
    plt.savefig(WalkerFolder + "k_off.png")
    plt.close()
    for i in range(num_walkers):
        plt.plot(list(range(0, num_steps)), walkers[i, :, 2])
        plt.title("kx_s_p values ")
    plt.savefig(WalkerFolder + "kx_s_p.png")
    plt.close()
    for i in range(num_walkers):
        plt.plot(list(range(0, num_steps)), walkers[i, :, 3])
        plt.title("b_e values ")
    plt.yscale('log')    
    plt.savefig(WalkerFolder + "b_e.png")
    plt.close()
    for i in range(num_walkers):
        plt.plot(list(range(0, num_steps)), walkers[i, :, 5])
        plt.title("kx_wt values ")
    plt.savefig(WalkerFolder + "kx_wt.png")
    plt.close()
    for i in range(num_walkers):
        plt.plot(list(range(0, num_steps)), walkers[i, :, 6])
        plt.title("kx_rne values ")
    plt.savefig(WalkerFolder + "kx_rne.png")
    plt.close()
    
    for i in range(num_walkers):
        plt.plot(list(range(0, num_steps)), walkers[i, :, 4])
        plt.title("b_ms values ")
    plt.savefig(WalkerFolder + "b_ms_wt.png")
    plt.close()
    for i in range(num_walkers):
        plt.plot(list(range(0, num_steps)), walkers[i, :, 7])
        plt.title("c values ")
        
    plt.savefig(WalkerFolder + "c.png")
    plt.close()
    for i in range(num_walkers):
        plt.plot(list(range(0, num_steps)), walkers[i, :, 8])
        plt.title("d values ")
    plt.savefig(WalkerFolder + "d.png")
    plt.close()
#    for i in range(num_walkers):
#        plt.plot(list(range(0, num_steps)), walkers[i, :, 13])
#        plt.title("b_m_rne values ")
#    plt.savefig(WalkerFolder + "b_ms_rne.png")
#    plt.close()
    

def heatMap():
    
    k_on_list = []
    k_off_list = []
    kx_s_p_list = []
    b_e_list = []
    kx_wt_list = []
    kx_rne_list = []
    b_ms_list = []
    c_list = []
    d_list = []
    
    dataFile = open(WalkerFolder + "walkers.csv", "r")
    csvReader = csv.reader(dataFile)
    walkers = np.zeros((num_walkers, num_steps, 9))
    counter = 0
    for row in csvReader:
        first = int(counter / num_steps)
        second = counter % num_steps
        walkers[first][second] = [float(i) for i in row]
        counter+=1
    likelihoodDataFile = open(WalkerFolder + "likelihood.csv", "r")
    csvReaderLikelihoods = csv.reader(likelihoodDataFile)
    likelihoods = np.zeros((num_walkers, num_steps))
    posteriorDataFile = open(WalkerFolder + "posterior.csv", "r")
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
        if (posteriorShaped[i] < -475):
            posteriorShaped[i] = -475
    
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
    while max_count < 10:
        
#        k_on,k_off,kx_s_p,b_e,b_ms,kx_wt,kx_rne,c,d = walkersShaped[:][np.argmax(posteriorShaped)]
        k_on,k_off,kx_s_p,b_e,c,d,kx_wt,kx_rne,b_ms = walkersShaped[:][np.argmax(posteriorShaped)]
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
            kx_rne_list.append(kx_wt)
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
        
    
    plt.scatter(walkersShaped[:,0],posteriorShaped)
    plt.xlim(np.min(walkersShaped[:,0]),np.max(walkersShaped[:,0]))
    plt.xscale('log')
    plt.title("k_on likelihood ")
    plt.savefig(WalkerFolder + "k_on_like.png")
    plt.close()
    
    plt.scatter(walkersShaped[:,1],posteriorShaped)
    plt.xlim(np.min(walkersShaped[:,1]),np.max(walkersShaped[:,1]))
    plt.xscale('log')
    plt.title("k_off likelihood ")
    plt.savefig(WalkerFolder + "k_off_like.png")
    plt.close()
    
    plt.scatter(walkersShaped[:,2],posteriorShaped)
    plt.xlim(np.min(walkersShaped[:,2]),np.max(walkersShaped[:,2]))
    plt.title("kx_s_p likelihood ")
    plt.savefig(WalkerFolder + "kx_s_p_like.png")
    plt.close()
    
    plt.scatter(walkersShaped[:,3],posteriorShaped)
    plt.xlim(np.min(walkersShaped[:,3]),np.max(walkersShaped[:,3]))
    plt.xscale('log')
    plt.title("b_e likelihood ")
    plt.savefig(WalkerFolder + "b_e_like.png")
    plt.close()
    
    
    plt.scatter(walkersShaped[:,5],posteriorShaped)
    plt.xlim(np.min(walkersShaped[:,5]),np.max(walkersShaped[:,5]))
    plt.title("kx_wt likelihood ")
    plt.savefig(WalkerFolder + "kx_wt_like.png")
    plt.close()
    
    plt.scatter(walkersShaped[:,6],posteriorShaped)
    plt.xlim(np.min(walkersShaped[:,6]),np.max(walkersShaped[:,6]))
    plt.title("kx_rne likelihood ")
    plt.savefig(WalkerFolder + "kx_rne_like.png")
    plt.close()
    
    plt.scatter(walkersShaped[:,4],posteriorShaped)
    plt.xlim(np.min(walkersShaped[:,4]),np.max(walkersShaped[:,4]))
    plt.xscale('log')
    plt.title("b_ms likelihood ")
    plt.savefig(WalkerFolder + "b_ms_like.png")
    plt.close()
    
    plt.scatter(walkersShaped[:,7],posteriorShaped)
    plt.xlim(np.min(walkersShaped[:,7]),np.max(walkersShaped[:,7]))
    plt.title("c likelihood ")
    plt.savefig(WalkerFolder + "c_like.png")
    plt.close()
    
    plt.scatter(walkersShaped[:,8],posteriorShaped)
    plt.xlim(np.min(walkersShaped[:,8]),np.max(walkersShaped[:,8]))
    plt.title("d likelihood ")
    plt.savefig(WalkerFolder + "d_like.png")
    plt.close()
    
#    plt.scatter(walkersShaped[:,13],posteriorShaped)
#    plt.xlim(np.min(walkersShaped[:,13]),np.max(walkersShaped[:,13]))
#    plt.title("b_ms_rne likelihood ")
#    plt.savefig(WalkerFolder + "b_ms_rne_like.png")
#    plt.close()
    

    
    walker_kA_likelihood = []
    walker_k_init_rne_start_likelihood = []
    walker_kx_s_p_likelihood = []
    walker_k_init_wt_likelihood = []
    walker_k_init_rne_likelihood = []
    walker_b_e_likelihood = []
    walker_kx_wt_likelihood = []
    walker_kx_rne_likelihood = []
    
    # [ ka, kd, kx_s, k_init_wt, k_init_rne, b_e, kx_wt, kx_rne ]
    
    for i in range(len(walkersShaped)):
        #kA
        if ((walkersShaped[i][1] - k_init_rne_start_ml) < .1 and (walkersShaped[i][2] - kx_s_p_ml) < .05 and (walkersShaped[i][3] - k_init_wt_ml) < 500 and  (walkersShaped[i][4] - k_init_rne_ml) < 500 and  (walkersShaped[5][0] - b_e_ml) < .0005 and  (walkersShaped[i][6] - kx_wt_ml) < .0005 and  (walkersShaped[i][7] - kx_rne_ml) < .0005):
            walker_kA_likelihood.append([walkersShaped[i][0],likelihoodsShaped2[i]])
            
        #k_init_rne_start
        if ((walkersShaped[i][0] - kA_ml) < 1E-9 and (walkersShaped[i][2] - kx_s_p_ml) < .01 and (walkersShaped[i][3] - k_init_wt_ml) < 100 and  (walkersShaped[i][4] - k_init_rne_ml) < 100 and  (walkersShaped[5][0] - b_e_ml) < .0001 and  (walkersShaped[i][6] - kx_wt_ml) < .0005 and  (walkersShaped[i][7] - kx_rne_ml) < .0005):
            walker_k_init_rne_start_likelihood.append([walkersShaped[i][1],likelihoodsShaped2[i]])
        
        #ks_s_p
        if ((walkersShaped[i][1] - k_init_rne_start_ml) < .1 and (walkersShaped[i][0] - kA_ml) < 1E-7 and (walkersShaped[i][3] - k_init_wt_ml) < 100 and  (walkersShaped[i][4] - k_init_rne_ml) < 100 and  (walkersShaped[5][0] - b_e_ml) < .0001 and  (walkersShaped[i][6] - kx_wt_ml) < .0005 and  (walkersShaped[i][7] - kx_rne_ml) < .0005):
            walker_kx_s_p_likelihood.append([walkersShaped[i][2],likelihoodsShaped2[i]])
    
        #k_init_wt
        if ((walkersShaped[i][1] - k_init_rne_start_ml) < .1 and (walkersShaped[i][2] - kx_s_p_ml) < .01 and (walkersShaped[i][0] - kA_ml) < 1E-7 and  (walkersShaped[i][4] - k_init_rne_ml) < 100 and  (walkersShaped[5][0] - b_e_ml) < .0001 and  (walkersShaped[i][6] - kx_wt_ml) < .0005 and  (walkersShaped[i][7] - kx_rne_ml) < .0005):
            walker_k_init_wt_likelihood.append([walkersShaped[i][3],likelihoodsShaped2[i]])
   
        #k_init_rne
        if ((walkersShaped[i][1] - k_init_rne_start_ml) < .1 and (walkersShaped[i][2] - kx_s_p_ml) < .01 and (walkersShaped[i][3] - k_init_wt_ml) < 100 and  (walkersShaped[i][0] - kA_ml) < 1E-7 and  (walkersShaped[5][0] - b_e_ml) < .0001 and  (walkersShaped[i][6] - kx_wt_ml) < .0005 and  (walkersShaped[i][7] - kx_rne_ml) < .0005):
            walker_k_init_rne_likelihood.append([walkersShaped[i][4],likelihoodsShaped2[i]])
    
        #b_e
        if ((walkersShaped[i][1] - k_init_rne_start_ml) < .1 and (walkersShaped[i][2] - kx_s_p_ml) < .01 and (walkersShaped[i][3] - k_init_wt_ml) < 100 and  (walkersShaped[i][4] - k_init_rne_ml) < 100 and  (walkersShaped[i][0] - kA_ml) < 1E-7 and  (walkersShaped[i][6] - kx_wt_ml) < .0005 and  (walkersShaped[i][7] - kx_rne_ml) < .0005):
            walker_b_e_likelihood.append([walkersShaped[i][5],likelihoodsShaped2[i]])
    
        #kx_wt
        if ((walkersShaped[i][1] - k_init_rne_start_ml) < .1 and (walkersShaped[i][2] - kx_s_p_ml) < .01 and (walkersShaped[i][3] - k_init_wt_ml) < 100 and  (walkersShaped[i][4] - k_init_rne_ml) < 100 and  (walkersShaped[5][0] - b_e_ml) < .0001 and  (walkersShaped[i][0] - kA_ml) < 1E-7 and  (walkersShaped[i][7] - kx_rne_ml) < .0005):
            walker_kx_wt_likelihood.append([walkersShaped[i][6],likelihoodsShaped2[i]])
    
        #kx_rne
        if ((walkersShaped[i][1] - k_init_rne_start_ml) < .1 and (walkersShaped[i][2] - kx_s_p_ml) < .01 and (walkersShaped[i][3] - k_init_wt_ml) < 100 and  (walkersShaped[i][4] - k_init_rne_ml) < 100 and  (walkersShaped[5][0] - b_e_ml) < .0001 and  (walkersShaped[i][6] - kx_wt_ml) < .0005 and  (walkersShaped[i][0] - kA_ml) < 1E-7):
            walker_kx_rne_likelihood.append([walkersShaped[i][7],likelihoodsShaped2[i]])
    
    walkersShaped_kA = np.array(sorted(walker_kA_likelihood,key=lambda x: x[0]))
    walkersShaped_k_init_rne_start = np.array(sorted(walker_k_init_rne_start_likelihood,key=lambda x: x[0]))
    walkersShaped_kx_s_p = np.array(sorted(walker_kx_s_p_likelihood,key=lambda x: x[0]))
    walkersShaped_k_init_wt = np.array(sorted(walker_k_init_wt_likelihood,key=lambda x: x[0]))
    walkersShaped_k_init_rne = np.array(sorted(walker_k_init_rne_likelihood,key=lambda x: x[0]))
    walkersShaped_b_e = np.array(sorted(walker_b_e_likelihood,key=lambda x: x[0]))
    walkersShaped_kx_wt = np.array(sorted(walker_kx_wt_likelihood,key=lambda x: x[0]))
    walkersShaped_kx_rne = np.array(sorted(walker_kx_rne_likelihood,key=lambda x: x[0]))
    

    xs_k_init_wt = []
    ys_k_init_wt = []
   
    width_am_wt = 10
    for i in range(int(min(walkersShaped_k_init_wt[:,0])), int(max(walkersShaped_k_init_wt[:,0])), width_am_wt):
        df = []
        for j in range(len(walkersShaped_k_init_wt)):
            if walkersShaped_k_init_wt[j,0] >= i-1 and walkersShaped_k_init_wt[j,0] < i +width_am_wt +1 :
                df.append(walkersShaped_k_init_wt[j,1]) # vertical slice
                
            #if np.isnan(np.min(func(_df.val)):            # ignore nans
            #   continue
        if not(df) :
            continue
        xs_k_init_wt.append(i + width_am_wt)                         
        ys_k_init_wt.append(np.max(df))
        
    xs_k_init_rne = []
    ys_k_init_rne = []
   
    width_am_rne = 10
    for i in range(int(min(walkersShaped_k_init_rne[:,0])), int(max(walkersShaped_k_init_rne[:,0])), width_am_rne):
        df = []
        for j in range(len(walkersShaped_k_init_rne)):
            if walkersShaped_k_init_rne[j,0] >= i-1 and walkersShaped_k_init_rne[j,0] < i +width_am_rne +1 :
                df.append(walkersShaped_k_init_rne[j,1]) # vertical slice
                
            #if np.isnan(np.min(func(_df.val)):            # ignore nans
            #   continue
        if not(df) :
            continue
        xs_k_init_rne.append(i + width_am_rne)                         
        ys_k_init_rne.append(np.max(df))    
        
    xs_kx_wt = []
    ys_kx_wt = []
   
    width_kx_wt = 10
    #kx_val = min(walkersShaped_kx[:,0])
    for i in range(0,len(walkersShaped_kx_wt)-width_kx_wt,width_kx_wt):
        df = []
        
        for j in range(len(walkersShaped_kx_wt)):
            if walkersShaped_kx_wt[j,0] >= walkersShaped_kx_wt[i,0] and walkersShaped_kx_wt[j,0] < walkersShaped_kx_wt[i+width_kx_wt,0] :
                df.append(walkersShaped_kx_wt[j,1]) # vertical slice
                
            #if np.isnan(np.min(func(_df.val)):            # ignore nans
            #   continue
        if not(df) :
            continue
        xs_kx_wt.append(walkersShaped_kx_wt[i+width_kx_wt,0])                         
        ys_kx_wt.append(np.max(df))   
   
    xs_kx_rne = []
    ys_kx_rne = []
   
    width_kx_rne = 10
    #kx_val = min(walkersShaped_kx[:,0])
    for i in range(0,len(walkersShaped_kx_rne)-width_kx_rne,width_kx_rne):
        df = []
        
        for j in range(len(walkersShaped_kx_rne)):
            if walkersShaped_kx_rne[j,0] >= walkersShaped_kx_rne[i,0] and walkersShaped_kx_rne[j,0] < walkersShaped_kx_rne[i+width_kx_rne,0] :
                df.append(walkersShaped_kx_rne[j,1]) # vertical slice
                
            #if np.isnan(np.min(func(_df.val)):            # ignore nans
            #   continue
        if not(df) :
            continue
        xs_kx_rne.append(walkersShaped_kx_rne[i+width_kx_rne,0])                         
        ys_kx_rne.append(np.max(df))
        
    xs_kA = []
    ys_kA = []
   
    width_kA = 1
    #kx_val = min(walkersShaped_kx[:,0])
    for i in range(0,len(walkersShaped_kA)-width_kA,width_kA):
        df = []
        
        for j in range(len(walkersShaped_kA)):
            if walkersShaped_kA[j,0] >= walkersShaped_kA[i,0] and walkersShaped_kA[j,0] < walkersShaped_kA[i+width_kA,0] :
                df.append(walkersShaped_kA[j,1]) # vertical slice
                
            #if np.isnan(np.min(func(_df.val)):            # ignore nans
            #   continue
        if not(df) :
            continue
        xs_kA.append(walkersShaped_kA[i+width_kA,0])                         
        ys_kA.append(np.max(df))
        
    xs_k_init_rne_start = []
    ys_k_init_rne_start = []
   
    width_k_init_rne_start = 10
    #kx_val = min(walkersShaped_kx[:,0])
    for i in range(0,len(walkersShaped_k_init_rne_start)-width_k_init_rne_start,width_k_init_rne_start):
        df = []
        
        for j in range(len(walkersShaped_k_init_rne_start)):
            if walkersShaped_k_init_rne_start[j,0] >= walkersShaped_k_init_rne_start[i,0] and walkersShaped_k_init_rne_start[j,0] < walkersShaped_k_init_rne_start[i+width_k_init_rne_start,0] :
                df.append(walkersShaped_k_init_rne_start[j,1]) # vertical slice
                
            #if np.isnan(np.min(func(_df.val)):            # ignore nans
            #   continue
        if not(df) :
            continue
        xs_k_init_rne_start.append(walkersShaped_k_init_rne_start[i+width_k_init_rne_start,0])                         
        ys_k_init_rne_start.append(np.max(df))
        
    xs_kx_s_p = []
    ys_kx_s_p = []
   
    width_kx_s_p = 10
    #kx_val = min(walkersShaped_kx[:,0])
    for i in range(0,len(walkersShaped_kx_s_p)-width_kx_s_p,width_kx_s_p):
        df = []
        
        for j in range(len(walkersShaped_kx_s_p)):
            if walkersShaped_kx_s_p[j,0] >= walkersShaped_kx_s_p[i,0] and walkersShaped_kx_s_p[j,0] < walkersShaped_kx_s_p[i+width_kx_s_p,0] :
                df.append(walkersShaped_kx_s_p[j,1]) # vertical slice
                
            #if np.isnan(np.min(func(_df.val)):            # ignore nans
            #   continue
        if not(df) :
            continue
        xs_kx_s_p.append(walkersShaped_kx_s_p[i+width_kx_s_p,0])                         
        ys_kx_s_p.append(np.max(df))
        
    xs_b_e = []
    ys_b_e = []
   
    width_b_e = 10
    #kx_val = min(walkersShaped_kx[:,0])
    for i in range(0,len(walkersShaped_b_e)-width_b_e,width_b_e):
        df = []
        
        for j in range(len(walkersShaped_b_e)):
            if walkersShaped_b_e[j,0] >= walkersShaped_b_e[i,0] and walkersShaped_b_e[j,0] < walkersShaped_b_e[i+width_b_e,0] :
                df.append(walkersShaped_b_e[j,1]) # vertical slice
                
            #if np.isnan(np.min(func(_df.val)):            # ignore nans
            #   continue
        if not(df) :
            continue
        xs_b_e.append(walkersShaped_b_e[i+width_b_e,0])                         
        ys_b_e.append(np.max(df))   
    
    # [ ka, kd, kx_s, k_init_wt, k_init_rne, b_e, kx_wt, kx_rne ]
    
    
    one_sd_below_kA = []
    one_sd_below_k_init_rne_start = []
    one_sd_below_kx_s_p = []
    one_sd_below_k_init_wt = []
    one_sd_below_k_init_rne = []
    one_sd_below_b_e = []
    one_sd_below_kx_wt = []
    one_sd_below_kx_rne = []
    
    for i in range(len(xs_kA)):
        if ( (max(walkersShaped_kA[:,1]) -.4945 -.1) <ys_kA[i] < (max(walkersShaped_kA[:,1]) -.4945+.1)):
            one_sd_below_kA.append([xs_kA[i],abs(.4945-(max(walkersShaped_kA[:,1])-ys_kA[i]))])
       
            
    for i in range(len(xs_k_init_rne_start)):
        if ( (max(walkersShaped_k_init_rne_start[:,1]) -.4945 -.1) <ys_k_init_rne_start[i] < (max(walkersShaped_k_init_rne_start[:,1]) -.4945+.1)):
            one_sd_below_k_init_rne_start.append([xs_k_init_rne_start[i],abs(.4945-(max(walkersShaped_k_init_rne_start[:,1])-ys_k_init_rne_start[i]))])
            
    for i in range(len(xs_kx_s_p)):
        if ( (max(walkersShaped_kx_s_p[:,1]) -.4945 -.1) <ys_kx_s_p[i] < (max(walkersShaped_kx_s_p[:,1]) -.4945+.1)):
            one_sd_below_kx_s_p.append([xs_kx_s_p[i],abs(.4945-(max(walkersShaped_kx_s_p[:,1])-ys_kx_s_p[i]))])        
            
    for i in range(len(xs_k_init_wt)):
        if ( (max(walkersShaped_k_init_wt[:,1]) -.4945 -.1) <ys_k_init_wt[i] < (max(walkersShaped_k_init_wt[:,1]) -.4945+.1)):
            one_sd_below_k_init_wt.append([xs_k_init_wt[i],abs(.4945-(max(walkersShaped_k_init_wt[:,1])-ys_k_init_wt[i]))])        
            
    for i in range(len(xs_k_init_rne)):
        if ( (max(walkersShaped_k_init_rne[:,1]) -.4945 -.1) <ys_k_init_rne[i] < (max(walkersShaped_k_init_rne[:,1]) -.4945+.1)):
            one_sd_below_k_init_rne.append([xs_k_init_rne[i],abs(.4945-(max(walkersShaped_k_init_rne[:,1])-ys_k_init_rne[i]))])       
            
    for i in range(len(xs_b_e)):
        if ( (max(walkersShaped_b_e[:,1]) -.4945 -.1) <ys_b_e[i] < (max(walkersShaped_b_e[:,1]) -.4945+.1)):
            one_sd_below_b_e.append([xs_b_e[i],abs(.4945-(max(walkersShaped_b_e[:,1])-ys_b_e[i]))])
            
    for i in range(len(xs_kx_wt)):
        if ( (max(walkersShaped_kx_wt[:,1]) -.4945 -.1) <ys_kx_wt[i] < (max(walkersShaped_kx_wt[:,1]) -.4945+.1)):
            one_sd_below_kx_wt.append([xs_kx_wt[i],abs(.4945-(max(walkersShaped_kx_wt[:,1])-ys_kx_wt[i]))]) 
            
    for i in range(len(xs_kx_rne)):
        if ( (max(walkersShaped_kx_rne[:,1]) -.4945 -.1) <ys_kx_rne[i] < (max(walkersShaped_kx_rne[:,1]) -.4945+.1)):
            one_sd_below_kx_rne.append([xs_kx_rne[i],abs(.4945-(max(walkersShaped_kx_rne[:,1])-ys_kx_rne[i]))])        
            
            
    # [ ka, kd, kx_s_p, k_init_wt, k_init_rne, b_e, kx_wt, kx_rne ]        
    
    
    k_init_wt_low_candidates = []
    k_init_wt_high_candidates = []
    
    for i in range(len(one_sd_below_k_init_wt)):
        if one_sd_below_k_init_wt[i][0] < k_init_wt_ml:
            k_init_wt_low_candidates.append([one_sd_below_k_init_wt[i][0],one_sd_below_k_init_wt[i][1]])
        elif one_sd_below_k_init_wt[i][0] > k_init_wt_ml:
            k_init_wt_high_candidates.append([one_sd_below_k_init_wt[i][0],one_sd_below_k_init_wt[i][1]])
    
    kA_low_candidates = []
    kA_high_candidates = []
    
    
    for i in range(len(one_sd_below_kA)):
        if one_sd_below_kA[i][0] < kA_ml:
            kA_low_candidates.append([one_sd_below_kA[i][0],one_sd_below_kA[i][1]])
        elif one_sd_below_kA[i][0] > kA_ml:
            kA_high_candidates.append([one_sd_below_kA[i][0],one_sd_below_kA[i][1]])
            
            
    k_init_rne_start_low_candidates = []
    k_init_rne_start_high_candidates = []
    
    
    for i in range(len(one_sd_below_k_init_rne_start)):
        if one_sd_below_k_init_rne_start[i][0] < k_init_rne_start_ml:
            k_init_rne_start_low_candidates.append([one_sd_below_k_init_rne_start[i][0],one_sd_below_k_init_rne_start[i][1]])
        elif one_sd_below_k_init_rne_start[i][0] > k_init_rne_start_ml:
            k_init_rne_start_high_candidates.append([one_sd_below_k_init_rne_start[i][0],one_sd_below_k_init_rne_start[i][1]])      
            
    kx_s_p_low_candidates = []
    kx_s_p_high_candidates = []
    
    
    for i in range(len(one_sd_below_kx_s_p)):
        if one_sd_below_kx_s_p[i][0] < kx_s_p_ml:
            kx_s_p_low_candidates.append([one_sd_below_kx_s_p[i][0],one_sd_below_kx_s_p[i][1]])
        elif one_sd_below_kx_s_p[i][0] > kx_s_p_ml:
            kx_s_p_high_candidates.append([one_sd_below_kx_s_p[i][0],one_sd_below_kx_s_p[i][1]])        
            
    k_init_rne_low_candidates = []
    k_init_rne_high_candidates = []
    
    
    for i in range(len(one_sd_below_k_init_rne)):
        if one_sd_below_k_init_rne[i][0] < k_init_rne_ml:
            k_init_rne_low_candidates.append([one_sd_below_k_init_rne[i][0],one_sd_below_k_init_rne[i][1]])
        elif one_sd_below_k_init_rne[i][0] > k_init_rne_ml:
            k_init_rne_high_candidates.append([one_sd_below_k_init_rne[i][0],one_sd_below_k_init_rne[i][1]])
            
    b_e_low_candidates = []
    b_e_high_candidates = []
    
    
    for i in range(len(one_sd_below_b_e)):
        if one_sd_below_b_e[i][0] < b_e_ml:
            b_e_low_candidates.append([one_sd_below_b_e[i][0],one_sd_below_b_e[i][1]])
        elif one_sd_below_b_e[i][0] > b_e_ml:
            b_e_high_candidates.append([one_sd_below_b_e[i][0],one_sd_below_b_e[i][1]])         
            
    kx_wt_low_candidates = []
    kx_wt_high_candidates = []
    
    
    for i in range(len(one_sd_below_kx_wt)):
        if one_sd_below_kx_wt[i][0] < kx_wt_ml:
            kx_wt_low_candidates.append([one_sd_below_kx_wt[i][0],one_sd_below_kx_wt[i][1]])
        elif one_sd_below_kx_wt[i][0] > kx_wt_ml:
            kx_wt_high_candidates.append([one_sd_below_kx_wt[i][0],one_sd_below_kx_wt[i][1]])       
    
    kx_rne_low_candidates = []
    kx_rne_high_candidates = []
    
    
    for i in range(len(one_sd_below_kx_rne)):
        if one_sd_below_kx_rne[i][0] < kx_rne_ml:
            kx_rne_low_candidates.append([one_sd_below_kx_rne[i][0],one_sd_below_kx_rne[i][1]])
        elif one_sd_below_kx_rne[i][0] > kx_rne_ml:
            kx_rne_high_candidates.append([one_sd_below_kx_rne[i][0],one_sd_below_kx_rne[i][1]])       
    
    # [ ka, kd, kx_s_p, k_init_wt, k_init_rne, b_e, kx_wt, kx_rne ] 
    kA_low_candidates = np.array(kA_low_candidates)  
    kA_high_candidates = np.array(kA_high_candidates)   

    kA_low = kA_low_candidates[np.argmin(kA_low_candidates[:,1]),0]
    kA_high = kA_high_candidates[np.argmin(kA_high_candidates[:,1]),0]
    
    kA_sd = np.mean([kA_ml-kA_low,kA_high-kA_ml])
    
    k_init_rne_start_low_candidates = np.array(k_init_rne_start_low_candidates)  
    k_init_rne_start_high_candidates = np.array(k_init_rne_start_high_candidates)   

    k_init_rne_start_low = k_init_rne_start_low_candidates[np.argmin(k_init_rne_start_low_candidates[:,1]),0]
    k_init_rne_start_high = k_init_rne_start_high_candidates[np.argmin(k_init_rne_start_high_candidates[:,1]),0]
    
    k_init_rne_start_sd = np.mean([k_init_rne_start_ml-k_init_rne_start_low,k_init_rne_start_high-k_init_rne_start_ml])
    
    
    kx_s_p_low_candidates = np.array(kx_s_p_low_candidates)  
    kx_s_p_high_candidates = np.array(kx_s_p_high_candidates)   

    kx_s_p_low = kx_s_p_low_candidates[np.argmin(kx_s_p_low_candidates[:,1]),0]
    kx_s_p_high = kx_s_p_high_candidates[np.argmin(kx_s_p_high_candidates[:,1]),0]
    
    kx_s_p_sd = np.mean([kx_s_p_ml-kx_s_p_low,kx_s_p_high-kx_s_p_ml])
    
    k_init_wt_low_candidates = np.array(k_init_wt_low_candidates)  
    k_init_wt_high_candidates = np.array(k_init_wt_high_candidates)   

    k_init_wt_low = k_init_wt_low_candidates[np.argmin(k_init_wt_low_candidates[:,1]),0]
    k_init_wt_high = k_init_wt_high_candidates[np.argmin(k_init_wt_high_candidates[:,1]),0]
    
    k_init_wt_sd = np.mean([k_init_wt_ml-k_init_wt_low,k_init_wt_high-k_init_wt_ml])
    
    k_init_rne_low_candidates = np.array(k_init_rne_low_candidates)  
    k_init_rne_high_candidates = np.array(k_init_rne_high_candidates)   

    k_init_rne_low = k_init_rne_low_candidates[np.argmin(k_init_rne_low_candidates[:,1]),0]
    k_init_rne_high = k_init_rne_high_candidates[np.argmin(k_init_rne_high_candidates[:,1]),0]
    
    k_init_rne_sd = np.mean([k_init_rne_ml-k_init_rne_low,k_init_rne_high-k_init_rne_ml])
    
    b_e_low_candidates = np.array(b_e_low_candidates)  
    b_e_high_candidates = np.array(b_e_high_candidates)   

    b_e_low = b_e_low_candidates[np.argmin(b_e_low_candidates[:,1]),0]
    b_e_high = b_e_high_candidates[np.argmin(b_e_high_candidates[:,1]),0]
    
    b_e_sd = np.mean([b_e_ml-b_e_low,b_e_high-b_e_ml])
    
    kx_wt_low_candidates = np.array(kx_wt_low_candidates)  
    kx_wt_high_candidates = np.array(kx_wt_high_candidates)   

    kx_wt_low = kx_wt_low_candidates[np.argmin(kx_wt_low_candidates[:,1]),0]
    kx_wt_high = kx_wt_high_candidates[np.argmin(kx_wt_high_candidates[:,1]),0]
    
    kx_wt_sd = np.mean([kx_wt_ml-kx_wt_low,kx_wt_high-kx_wt_ml])
    
    kx_rne_low_candidates = np.array(kx_rne_low_candidates)  
    kx_rne_high_candidates = np.array(kx_rne_high_candidates)   

    kx_rne_low = kx_rne_low_candidates[np.argmin(kx_rne_low_candidates[:,1]),0]
    kx_rne_high = kx_rne_high_candidates[np.argmin(kx_rne_high_candidates[:,1]),0]
    
    kx_rne_sd = np.mean([kx_rne_ml-kx_rne_low,kx_rne_high-kx_rne_ml])
    
    
    labels=["kA", "k_init_rne_start", "kx_s_p", "k_init_wt", "k_init_rne", "b_e", "kx_wt", "kx_rne"]
    ml_labels = [kA_ml, k_init_rne_start_ml, kx_s_p_ml, k_init_wt_ml, k_init_rne_ml, b_e_ml, kx_wt_ml, kx_rne_ml]
    sd_labels = [kA_sd,k_init_rne_start_sd,kx_s_p_sd,k_init_wt_sd,k_init_rne_sd,b_e_sd,kx_wt_sd,kx_rne_sd]
    
    for par1 in range(7):
        for par2 in range(par1 + 1, 8):
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["blue", "violet","red"])
            plt.xlabel(labels[par1])
            plt.ylabel(labels[par2])
            #title = (labels[par1] + "_ml = %.5f, " + labels[par1] +"_sd = %.5f, " + labels[par2] + "_ml = %.5f, " + labels[par2] + "_sd = %.5f")% (ml_labels[par1],sd_labels[par1],ml_labels[par2],sd_labels[par2])
            #print(title)
            #plt.title(title)
            plt.scatter(walkersShaped[:, par1], walkersShaped[:, par2], c = likelihoodsShaped, cmap = cmap)
            plt.axis([min(walkersShaped[:, par1]), np.percentile(walkersShaped[:, par1], 75), min(walkersShaped[:, par2]), max(walkersShaped[:, par2])])
            plt.colorbar()
            plt.savefig(WalkerFolder + labels[par1] + "_" + labels[par2] + ".png")
            plt.close()
 
  
postProcess()
heatMap()