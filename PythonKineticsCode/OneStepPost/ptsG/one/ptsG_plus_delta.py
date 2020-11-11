# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 18:52:48 2018

@author: reyer
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 19:30:07 2018

@author: reyer
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 17:29:43 2018

@author: reyer
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May  7 17:40:37 2018

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
import emcee
import corner
import scipy.optimize as op
from multiprocessing import cpu_count
import csv
from schwimmbad import MPIPool


xData = [0, 60, 180, 360, 720, 1080, 1440]

result = []


numReplicates = 2
timePoints = 8
numMeasurable = 6
#0 is p_wt, 1 is m_wt, 2 is s_wt, 3 p_rne, 4 m_rne, 5 s_rine
yData = np.zeros(shape=(numReplicates, numMeasurable, timePoints))

yData[0, 0, :] = [134917.1029, 184452.1771, 180579.3257,216245.1063, 432250.9255, 839879.6454, 1170873.494, 1903248.354]
yData[0, 1, :] = [22443.71086, 339380.1687, 774276.8132, 989754.3786, 1178364.294, 1168843.171, 1039977.762, 1152377.24]
yData[0, 2, :] = [1299487.304, 2158748.53, 1461797.642, 1479650.818, 1213052.738, 1132760.79, 930608.6544, 743058.2862]
yData[0, 3, :] = [257171.8156, 242054.4163, 248115.1836, 246027.831, 383643.0025, 802367.1203, 1131430.97, 1382649.746]
yData[0, 4, :] = [16225.77785, 7646.443305, 194166.3494, 714018.5826, 1418654.816, 1533096.423, 1617435.188, 2046367.339]
yData[0, 5, :] = [2188990.237, 2050903.278, 2390768.473, 1788277.574, 1296136.946, 1018016.032, 973140.0654, 973140.0654]

yData[1, 0, :] = [144202.1426, 333919.142, 149161.7576, 304584.7737, 738063.1426, 808832.9779, 1392651.226, 1681171.266]
yData[1, 1, :] = [33318.96697, 262752.4716, 627298.9573, 851325.6287, 1159651.072, 1151829.933, 1016580.917, 1118273.902]
yData[1, 2, :] = [1788327.436, 1787030.495, 1406293.534, 1188055.269, 1088649.355, 825880.4722, 785990.9093, 744898.0836]
yData[1, 3, :] = [64338.27187, 36900.77477, 47182.78534, 52802.52599, 298986.6166, 731084.9614, 1055270.669, 1940635.146]
yData[1, 4, :] = [28424.46802, 39375.74472, 237406.821, 1027314.934, 2003857.809, 2042997.469, 1758470.575, 2012571.706]
yData[1, 5, :] = [2547665.889, 2079076.002, 2324980.48, 2178616.805, 1944435.717, 956261.5957, 1099876.149, 1250870.797]

yData[:,5,:] = ((yData[:,5,:]+5.2879e+04)/7.3302e+03)
yData[:,2,:] = ((yData[:,2,:]+5.2879e+04)/7.3302e+03) #RyhB has 4 probes, SgrS has 9

yData[:,1,:] = ((yData[:,1,:]+1.3257e+05)/6.0368e+03)+1
yData[:,4,:] = ((yData[:,4,:]+1.3257e+05)/6.0368e+03)+1 #RyhB has 4 probes, SgrS has 9


#SRNA, MRNA, M_S, Protein
# initial_cond1 = [yS1_wt[0],yM1_wt[0],0,yP1_wt[3],yS1_rne[2],yM1_rne[2],0, np.mean([yP1_rne[4],yP1_rne[5]])]
# initial_cond2 = [yS2_wt[0],yM2_wt[0],0,yP2_wt[3],yS2_rne[2],yM2_rne[2],0, np.mean([yP2_rne[4],yP2_rne[5]])]
# print(initial_cond1)
# print(initial_cond2)

initial_cond = np.zeros(shape=(1, 8))
for i in range(1):
    initial_cond[i, :] = [yData[i, 2, 0], yData[i, 1, 0],0, yData[i, 0, 3], yData[i, 5, 0], yData[i, 4, 0],0, yData[i, 3, 3]]

a_m_wt = 1.905
a_m_wt_minus =  1.905
a_m_wt_co = 1.905
sigma_a_m_wt = 0.305

a_m_rne = 1.431
a_m_rne_minus = 1.431
a_m_rne_co = 1.431
sigma_a_m_rne = 0.031

b_m = 0.003229666667
b_m_ml = 0.003229666667
sigma_b_m = 0.00003108590249

k_on = 1.29E-4
k_off = 1E-2
k_on_nuc = 1E-4
k_off_nuc = 0.11

a_s_wt = 0.3646
a_s_wt_ml = 0.3646
sigma_a_s_wt = 0.02771858582
b_s_wt = 0.00147425
b_s_wt_ml = 0.00147425
sigma_b_s_wt = 0.0002935200249

a_s_rne = 0.37649
a_s_rne_ml = 0.37649
sigma_a_s_rne = 0.05064298767
b_s_rne = 0.001022905
b_s_rne_ml = 0.001022905
sigma_b_s_rne = 0.0001349089028

kx_wt = 13.7033
kx_rne = 5.45035
kx_wt_ml = 13.7033
kx_rne_ml = 5.45035
sigma_kx_wt = 1.6967999999999996
sigma_kx_rne = 0.437349



kx_s_p = .358

b_e = 2.82E-03#GUESS
b_ms = .0025
b_ms_ml = .0025

b_p = 0.00018
lnf = -1.05


def eq(initial_conditions):
    #  #-time-grid-----------------------------------
    t = xData
    t_start = 0
    t_end = 1440
    t_step = 0.1
     # differential-eq-system----------------------
     # USES LSODA From fortran package and determines t dynamically

    

    def funct(t, y):
         # print(t)
         # print(y)

         s_wt=y[0]
         m_wt=y[1]
         ms_wt=y[2]
         p_wt=y[3]

         s_rne=y[4]
         m_rne=y[5]
         ms_rne=y[6]
         p_rne=y[7]
         
         


         f0 = a_s_wt - b_s_wt * s_wt - k_on * s_wt * m_wt + k_off * ms_wt 
         f1 = a_m_wt  - b_m * m_wt - k_on * s_wt * m_wt + k_off * ms_wt #dm/dt
         f2 = k_on * s_wt * m_wt - b_ms * ms_wt - b_e * ms_wt - k_off * ms_wt
         f3 = m_wt * kx_wt + (kx_s_p*kx_wt) * ms_wt  - p_wt * b_p

         f4 = a_s_rne - b_s_rne * s_rne - k_on * s_rne * m_rne + k_off * ms_rne
         f5 = a_m_rne  - b_m * m_rne - k_on * s_rne * m_rne + k_off * ms_rne#dm/dt
         f6 = k_on * s_rne * m_rne - b_ms * ms_rne - k_off * ms_rne
         f7 = m_rne * kx_rne + (kx_s_p*kx_rne) * ms_rne - p_rne * b_p

         #print(a_m(a_m_wt,kA,c,s_wt,ms_wt))

         return [f0, f1, f2, f3, f4, f5, f6, f7];
    ode = integrate.ode(funct)
    # BDF method suited to stiff systems of ODEs
    ode.set_integrator('lsoda', nsteps=50000, min_step=1E-25)
    filteredInitialConditions = initial_conditions
    #Set to some arbitraty starting conditions if the beginnings is nan

    if (np.isnan(initial_conditions[0]) or np.isnan(initial_conditions[2]) or np.isnan(initial_conditions[3])):
        filteredInitialConditions[0] = 1060613.113
        filteredInitialConditions[1] = 85078.09744
        filteredInitialConditions[2] = 0
        filteredInitialConditions[3] = 466524.1913
    if (np.isnan(initial_conditions[4]) or np.isnan(initial_conditions[6]) or np.isnan(initial_conditions[7])):
        filteredInitialConditions[0] = 2071941.54
        filteredInitialConditions[1] = 159857.6577
        filteredInitialConditions[2] = 0
        filteredInitialConditions[3] = 538447.8646

    ode.set_initial_value(initial_conditions,t_start)
    t = []
    ys = []

    while ode.successful() and ode.t < t_end:
        ys.append(ode.y)
        t.append(ode.t)
        ode.integrate(ode.t + t_step)
    ys = np.array(ys)
    # return 0
    return (ys[:,0],ys[:,1],ys[:,2],ys[:,3],ys[:,4],ys[:,5],ys[:,6],ys[:,7], t)



#In the original code, they pass x, y, yerr in the place of temp1-temp3
#We don't really care about these 3 parameters because our x, y and yerr are global variables
def lnlike(args, arg1, arg2, arg3):
    global k_on , k_off, kx_s_p, b_e, kx_wt, kx_rne, b_ms, a_m_wt, a_m_rne
    k_on = args[0] #GUESS
    k_off = args[1]
    kx_s_p = args[2] #GUESS
    b_e = args[3]
    a_m_wt = args[4]
    a_m_rne = args[5]

    kx_wt = arg1
    kx_rne = arg2
    b_ms = arg3
    
    s_wt = []
    m_wt = []
    ms_wt = []
    p_wt = []
    s_rne = []
    m_rne = []
    ms_rne = []
    p_rne = []
    t = []
    
    distance_mRNA = []
    distance_protein = []
    distance_sRNA = []
    inv_s2_mRNA = []
    inv_s2_protein = []
    inv_s2_sRNA = []
    distance = []
    for i in range(1):
        distance_mRNA.append(0)
        distance_protein.append(0)
        distance_sRNA.append(0)
        inv_s2_mRNA.append(0)
        inv_s2_protein.append(0)
        inv_s2_sRNA.append(0)
        distance.append(0)

        s_wt.append(0)
        m_wt.append(0)
        ms_wt.append(0)
        p_wt.append(0)
        s_rne.append(0)
        m_rne.append(0)
        ms_rne.append(0)
        p_rne.append(0)
        t.append(0)

        s_wt[i], m_wt[i], ms_wt[i], p_wt[i], s_rne[i], m_rne[i], ms_rne[i], p_rne[i], t[i] = eq(initial_cond[i])
        if (len(t[i])< 14400):
            return -np.inf



    x_RNA = [0,600,1800,3600,7200,10800,14400]

    x_protein_wt = [0,3600,7200,10800]
    x_protein_rne = [0,3600,7200,10800]
    real_RNA_time = [0,1,2,3,4,5,6]
    real_protein_time_wt = [3,4,5,6]
    real_protein_time_rne = [3,4,5,6]
    likeAll = [0, 0, 0, 0]
    
    for e,i in enumerate([0,3,4,5]):
        temp = yData[:, 1, :]
        mRNAData = temp[np.isfinite(temp)]
        temp = yData[:, 2, :]
        sRNAData = temp[np.isfinite(temp)]
        temp = yData[:, 0, :]
        proteinData = temp[np.isfinite(temp)]

        error_mRNA = np.std(mRNAData)
        error_sRNA = np.std(sRNAData)
        error_protein = np.std(proteinData)
        for r in range(1):
            #0 is p_wt, 1 is m_wt, 2 is s_wt, 3 p_rne, 4 m_rne, 5 s_rine
            if (np.isnan(yData[r, 1, real_RNA_time[i]]) or np.isnan(yData[r, 2, real_RNA_time[i]]) or np.isnan(yData[r, 0, real_protein_time_wt[e]])):
                continue
            distance_mRNA[r] = (m_wt[r][x_RNA[i]] + ms_wt[r][x_RNA[i]] - yData[r, 1, real_RNA_time[i]])
            distance_sRNA[r] = (s_wt[r][x_RNA[i]] + ms_wt[r][x_RNA[i]] - yData[r, 2, real_RNA_time[i]])
            distance_protein[r] = (p_wt[r][x_protein_wt[e]] - yData[r, 0, real_protein_time_wt[e]])
            inv_s2_mRNA[r] = 1/(error_mRNA + (m_wt[r][x_RNA[i]] + ms_wt[r][x_RNA[i]])**2*np.exp(2*lnf))
            inv_s2_protein[r] = 1/(error_protein + p_wt[r][x_protein_wt[e]]**2*np.exp(2*lnf))
            inv_s2_sRNA[r] = 1/(error_sRNA + (s_wt[r][x_RNA[i]] + ms_wt[r][x_RNA[i]])**2*np.exp(2*lnf))
            distance[r] = ((distance_mRNA[r]**2)*inv_s2_mRNA[r] - np.log(inv_s2_mRNA[r])) + ((distance_protein[r]**2)*inv_s2_protein[r] - np.log(inv_s2_protein[r])) + ((distance_sRNA[r]**2)*inv_s2_sRNA[r] - np.log(inv_s2_sRNA[r]))

        likeAll[0] += -0.5*(sum(distance))

    for i in [1,2,6]:

        temp = yData[:, 1, :]
        mRNAData = temp[np.isfinite(temp)]
        temp = yData[:, 2, :]
        sRNAData = temp[np.isfinite(temp)]

        error_mRNA = np.std(mRNAData)
        error_sRNA = np.std(sRNAData)
        for r in range(1):
            #0 is p_wt, 1 is m_wt, 2 is s_wt, 3 p_rne, 4 m_rne, 5 s_rine
            if (np.isnan(yData[r, 1, real_RNA_time[i]]) or np.isnan(yData[r, 2, real_RNA_time[i]])):
                continue
            distance_mRNA[r] = (m_wt[r][x_RNA[i]] + ms_wt[r][x_RNA[i]] - yData[r, 1, real_RNA_time[i]])
            distance_sRNA[r] = (s_wt[r][x_RNA[i]] + ms_wt[r][x_RNA[i]] - yData[r, 2, real_RNA_time[i]])
            distance_protein[r] = (p_wt[r][x_protein_wt[e]] - yData[r, 0, real_protein_time_wt[e]])

            inv_s2_mRNA[r] = 1/(error_mRNA + (m_wt[r][x_RNA[i]] + ms_wt[r][x_RNA[i]])**2*np.exp(2*lnf))
            inv_s2_sRNA[r] = 1/(error_sRNA + (s_wt[r][x_RNA[i]] + ms_wt[r][x_RNA[i]])**2*np.exp(2*lnf))
            distance[r] = ((distance_mRNA[r]**2)*inv_s2_mRNA[r] - np.log(inv_s2_mRNA[r])) + ((distance_sRNA[r]**2)*inv_s2_sRNA[r] - np.log(inv_s2_sRNA[r]))

        likeAll[1] += -0.5*(sum(distance))
    #RNE
    for e,i in enumerate([0,3,4,5]):
        temp = yData[:, 4, :]
        mRNAData = temp[np.isfinite(temp)]
        temp = yData[:, 5, :]
        sRNAData = temp[np.isfinite(temp)]
        temp = yData[:, 3, :]
        proteinData = temp[np.isfinite(temp)]

        error_mRNA = np.std(mRNAData)
        error_sRNA = np.std(sRNAData)
        error_protein = np.std(proteinData)

        for r in range(1):
            if (np.isnan(yData[r, 4, real_RNA_time[i]]) or np.isnan(yData[r, 5, real_RNA_time[i]]) or np.isnan(yData[r, 3, real_protein_time_rne[e]])):
                continue
            distance_mRNA[r] = (m_rne[r][x_RNA[i]] + ms_rne[r][x_RNA[i]] - yData[r, 4, real_RNA_time[i]])
            distance_sRNA[r] = (s_rne[r][x_RNA[i]] + ms_rne[r][x_RNA[i]] - yData[r, 5, real_RNA_time[i]])
            distance_protein[r] = (p_rne[r][x_protein_rne[e]] - yData[r, 3, real_protein_time_rne[e]])
            # print(m_rne[r][x_RNA[i]])
            # print(yData[r, 4, real_RNA_time[i]])

            inv_s2_mRNA[r] = 1/(error_mRNA + (m_rne[r][x_RNA[i]] + ms_rne[r][x_RNA[i]])**2*np.exp(2*lnf))
            inv_s2_protein[r] = 1/(error_protein + p_rne[r][x_protein_rne[e]]**2*np.exp(2*lnf))
            inv_s2_sRNA[r] = 1/(error_sRNA + (s_rne[r][x_RNA[i]] + ms_rne[r][x_RNA[i]])**2*np.exp(2*lnf))
            distance[r] = ((distance_mRNA[r]**2)*inv_s2_mRNA[r] - np.log(inv_s2_mRNA[r])) + ((distance_protein[r]**2)*inv_s2_protein[r] - np.log(inv_s2_protein[r])) + ((distance_sRNA[r]**2)*inv_s2_sRNA[r] - np.log(inv_s2_sRNA[r]))

        likeAll[2] += -0.5*(sum(distance))

    for i in [1,2,6]:
        temp = yData[:, 4, :]
        mRNAData = temp[np.isfinite(temp)]
        temp = yData[:, 5, :]
        sRNAData = temp[np.isfinite(temp)]
        error_mRNA = np.std(mRNAData)
        error_sRNA = np.std(sRNAData)
        for r in range(1):
            if (np.isnan(yData[r, 4, real_RNA_time[i]]) or np.isnan(yData[r, 5, real_RNA_time[i]])):
                continue
            distance_mRNA[r] = (m_rne[r][x_RNA[i]] + ms_rne[r][x_RNA[i]] - yData[r, 4, real_RNA_time[i]])
            distance_sRNA[r] = (s_rne[r][x_RNA[i]] + ms_rne[r][x_RNA[i]] - yData[r, 5, real_RNA_time[i]])
            inv_s2_mRNA[r] = 1/(error_mRNA + (m_rne[r][x_RNA[i]] + ms_rne[r][x_RNA[i]])**2*np.exp(2*lnf))
            inv_s2_sRNA[r] = 1/(error_sRNA + (s_rne[r][x_RNA[i]] + ms_rne[r][x_RNA[i]])**2*np.exp(2*lnf))
            distance[r] = ((distance_mRNA[r]**2)*inv_s2_mRNA[r] - np.log(inv_s2_mRNA[r])) + ((distance_sRNA[r]**2)*inv_s2_sRNA[r] - np.log(inv_s2_sRNA[r]))

        likeAll[3] += -0.5*(sum(distance))
    return sum(likeAll)


#lnprior. Function of likelihood od parameters
def lnprior(args):
    #Set conditions for args

    k_on = args[0] #GUESS
    k_off = args[1]
    kx_s_p = args[2] #GUESS
    b_e = args[3]
    a_m_wt = args[4]
    a_m_rne = args[5]

    kx_wt = args[6]
    kx_rne = args[7]
    b_ms = args[8]
    

        
    
    if not (1E-9 <= k_on <= 1E-1 and 1E-5 <= k_off <= 1 and .01 <= kx_s_p <= 1.10 and 1E-5 <= b_e <= 0.1 and kx_wt_ml - 1.01*sigma_kx_wt <= kx_wt <= kx_wt_ml + 1.01*sigma_kx_wt and kx_rne_ml - 1.01*sigma_kx_rne <= kx_rne <= kx_rne_ml + 1.01*sigma_kx_rne and a_m_rne_minus - sigma_a_m_rne - .010 <= a_m_rne <= a_m_rne_minus + sigma_a_m_rne + .010 and a_m_wt_minus - sigma_a_m_wt <= a_m_wt <= a_m_wt_minus + sigma_a_m_wt + .010 and b_m <= b_ms <= 0.8):
        return -np.inf
    
    def GaussError(mu,sigma,curr_val):
        return (np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5*(curr_val-mu)**2/sigma**2)


    return np.sum([GaussError(kx_wt_ml,sigma_kx_wt,kx_wt),GaussError(kx_rne_ml,sigma_kx_rne,kx_rne),GaussError(a_m_wt_minus,sigma_a_m_wt,a_m_wt),GaussError(a_m_rne_minus,sigma_a_m_rne,a_m_rne)])

#In the original code, they pass x, y, yerr in the place of temp1-temp3
#We don't really care about these 3 parameters because our x, y and yerr are global variables
def lnprob(args):
    # print(args)
    lp = lnprior(args)
    if not np.isfinite(lp):
        return -np.inf, -np.inf, -np.inf
    lnl = lnlike(args, args[6], args[7], args[8])
    if not np.isfinite(lnl):
        return -np.inf, -np.inf, -np.inf
    return lp + lnl, lnl, lp + lnl


def fixPos():
    for i in range(len(pos)):
        pos[i][0] = np.clip(pos[i][0], 1E-9, 1E-1) #k_on
        pos[i][1] = np.clip(pos[i][1], 1E-5, 1) #k_off
        pos[i][2] = np.clip(pos[i][2], .01, 1.10) #Kx_s_p
        pos[i][3] = np.clip(pos[i][3], 1E-5, .1)#b_e
        pos[i][4] = np.clip(pos[i][4], a_m_wt_minus - sigma_a_m_wt, a_m_wt_minus + sigma_a_m_wt)#a_m_wt
        pos[i][5] = np.clip(pos[i][5], a_m_rne_minus - sigma_a_m_rne, a_m_rne_minus + sigma_a_m_rne)#a_m_rne
        pos[i][6] = np.clip(pos[i][6], kx_wt_ml - sigma_kx_wt, kx_wt_ml + sigma_kx_wt)#kx_wt
        pos[i][7] = np.clip(pos[i][7], kx_rne_ml - sigma_kx_rne, kx_rne_ml + sigma_kx_rne)#kx_rne
        pos[i][8] = np.clip(pos[i][8], b_m_ml, .8)#b_m_s




def graph():
    global result, pos
    ndim, nwalkers = 9, 50
#    nll = lambda *args: -lnlike(*args)
#    bnds = ((1E-8, 1E-1), (.2, .99), (.4*a_m_wt_minus, a_m_wt_minus), (.4*a_m_rne_minus, a_m_rne_minus), (5E-4, 0.01), (.0012, .0016), (.0007, .0011))
#    result = op.minimize(nll, [kA, kx_s_p, a_m_wt, a_m_rne, b_e, b_s_wt, b_s_rne], args = (kx_wt_ml, kx_rne_ml), method ='TNC', bounds = bnds)
#    kA_ml, kx_s_p_ml, a_m_wt_ml, a_m_rne_ml, b_e_ml, b_s_wt_ml, b_s_rne_ml = result["x"]
    k_on_ml = 2.284E-4
    k_off_ml = 2.653E-2
    kx_s_p_ml = 0.411
    b_e_ml = 2.82E-03
    b_ms_ml = .00420
    pos = [[k_on_ml, k_off_ml, kx_s_p_ml, b_e_ml, a_m_wt_minus, a_m_rne_minus, kx_wt_ml, kx_rne_ml, b_ms_ml]
            + np.array([2E-3*np.random.randn(), 1E-1*np.random.randn(), .1*np.random.randn(), 1e-2*np.random.randn(), 1.0*np.random.randn(), 1.0*np.random.randn(), 1.0*np.random.randn(), 1.0*np.random.randn(), .001*np.random.randn()])
            for i in range(nwalkers)]
    dtype = [("likelihood", float), ("posterior", float)]
    fixPos()
    with MPIPool() as pool:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool = pool, blobs_dtype = dtype)
        sampler.run_mcmc(pos, 10000)
    samples = sampler.chain[:, :, :].reshape((-1, ndim))
    np.savetxt("walkers.csv", samples, delimiter = ',', fmt='%10.12f')
    blobs = sampler.get_blobs()
    likelihood_samps = blobs["likelihood"]
    posterior_samps = blobs["posterior"]
    np.savetxt("likelihood.csv", likelihood_samps, delimiter = ',', fmt='%10.12f') #200 x 20
    np.savetxt("posterior.csv", posterior_samps, delimiter = ',', fmt='%10.12f') #200 x 20
def postProcess():
    dataFile = open("walkers.csv", "r")
    csvReader = csv.reader(dataFile)
    walkers = np.zeros((100, 300, 9))
    counter = 0
    for row in csvReader:
        first = int(counter / 300)
        second = counter % 300
        walkers[first][second] = [float(i) for i in row]
        counter+=1
    for i in range(100):
        plt.plot(list(range(0, 300)), walkers[i, :, 0])
        plt.title("kx values ")
    plt.savefig("kx.png")
    plt.close()
    for i in range(100):
        plt.plot(list(range(0, 300)), walkers[i, :, 1])
        plt.title("a_m values ")
    plt.savefig("a_m.png")
    plt.close()


graph()
# postProcess()
# heatMap()
