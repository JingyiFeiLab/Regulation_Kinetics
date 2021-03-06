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

yData[1, 0, :] = [1130948.676,892617.3287,1074207.349,1375535.856,2229638.274,2729098.729,3567882.043,5582180.8]
yData[1, 1, :] = [10050.51839,305731.8887,803354.734,860573.0497,284796.5291,881811.9806,918839.9995,920764.8788]
yData[1, 2, :] = [430758.4504,399694.8762,385861.0661,425008.5325,117141.0959,317523.1091,392639.7222,311909.1325]
yData[0, 3, :] = [537540.133,489861.7415,495303.3963,636534.158,1067081.43,1469663.453,2574256.314,2543087.016]
yData[0, 4, :] = [15372.71674,28146.29857,440707.2556,1226991.127,1666440.507,1320825.78,1559866.142,344620.4189]
yData[0, 5, :] = [593375.0032,455116.0876,568228.2777,511043.6822,508376.8531,461330.1373,395811.775,230814.0326]


yData[0, 0, :] = [1077656.5,1109629.183,1219469.886,1555835.48,2574231.765,3197925.634,4967996.936,5582180.8]
yData[0, 1, :] = [4809.786906,405355.402,570542.6576,739668.6922,788078.7262,859802.2183,782590.881,739851.4549]
yData[0, 2, :] = [231573.46,373192.0521,272348.6771,291664.8198,234679.4617,264677.7101,288519.334,210534.048]
yData[1, 3, :] = [536904.3621,461661.7087,495223.7801,678235.2065,1227608.695,1859059.735,2832320.07,3746419.851]
yData[1, 4, :] = [19898.01018,123286.8915,428240.4877,1243484.439,1724017.649,1690148.639,1684720.668,1434147.074]
yData[1, 5, :] = [461429.9867,478100.4001,485829.6041,455485.2018,518848.1603,573561.3758,465289.356,276181.5914]



yData[:,5,:] = yData[:,5,:] *(9/4)
yData[:,2,:] = yData[:,2,:] *(9/4) #RyhB has 4 probes, SgrS has 9

yData[:,5,:] = ((yData[:,5,:]+5.2879e+04)/7.3302e+03)
yData[:,2,:] = ((yData[:,2,:]+5.2879e+04)/7.3302e+03) #RyhB has 4 probes, SgrS has 9

yData[:,1,:] = ((yData[:,1,:]+1.3257e+05)/6.0368e+03)+1
yData[:,4,:] = ((yData[:,4,:]+1.3257e+05)/6.0368e+03)+1 #RyhB has 4 probes, SgrS has 9


#SRNA, MRNA, M_S, Protein
# initial_cond1 = [yS1_wt[0],yM1_wt[0],0,yP1_wt[3],yS1_rne[2],yM1_rne[2],0, np.mean([yP1_rne[4],yP1_rne[5]])]
# initial_cond2 = [yS2_wt[0],yM2_wt[0],0,yP2_wt[3],yS2_rne[2],yM2_rne[2],0, np.mean([yP2_rne[4],yP2_rne[5]])]
# print(initial_cond1)
# print(initial_cond2)

initial_cond = np.zeros(shape=(1, 12))
for i in range(1):
    initial_cond[i, :] = [yData[i, 2, 0],0,0, yData[i, 1, 0],0, yData[i, 0, 3], yData[i, 5, 0],0,0, yData[i, 4, 0],0, yData[i, 3, 3]]


DNA = 20

k_init_wt = 0.10843
k_init_wt_ml = 0.10843
sigma_k_init_wt = 0.0191216082

k_init_rne = 0.1067575
k_init_rne_ml = 0.1067575
sigma_k_init_rne = 0.01219053834


k_elon = 0.0592
k_elon_prime = 0.0592

b_m = 0.005336285714
b_m_ml = 0.005336285714
sigma_b_m = 0.0002118173067

k_on = 3.23E-05
k_off = 2.69E-01
k_on_nuc = 0.0
k_off_nuc = 0.0

a_s_wt = 0.263325
a_s_wt_ml = 0.263325
sigma_a_s_wt = 0.06681451975/2
b_s_wt = 0.00279305
b_s_wt_ml = 0.00279305
sigma_b_s_wt = 0.0002329916844

a_s_rne = 0.204975
a_s_rne_ml = 0.204975
sigma_a_s_rne = 0.07343303923/2
b_s_rne = 0.0019644
b_s_rne_ml = 0.0019644
sigma_b_s_rne = 0.0003454923733

kx_wt = 26.60793333
kx_rne = 10.355625
kx_wt_ml = 26.60793333
kx_rne_ml = 10.355625
sigma_kx_wt = 3.952737716
sigma_kx_rne = 1.757157171


kx_s_p = .358

b_e = 1.64E-03#GUESS
b_ms = .0065
b_ms_ml = .0065

b_p = 0.00018
c = 0.5
d = 0.5




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
         mi_wt = y[1]
         mis_wt = y[2]
         m_wt=y[3]
         ms_wt=y[4]
         p_wt=y[5]

         s_rne=y[6]
         mi_rne = y[7]
         mis_rne = y[8]
         m_rne=y[9]
         ms_rne=y[10]
         p_rne=y[11]
         
         


         f0 = a_s_wt - b_s_wt * s_wt - k_on * s_wt * (m_wt) - k_on_nuc * s_wt*mi_wt + k_off * (ms_wt) + k_off_nuc*mis_wt #ds/dt contributes sRNA signal
         f1 = k_init_wt*DNA - k_elon*mi_wt - k_on_nuc * s_wt * mi_wt + k_off_nuc*mis_wt#dm_i/dt 
         f2 = k_on_nuc * s_wt * mi_wt - k_elon_prime * mis_wt - (b_e + b_ms) * mis_wt - k_off_nuc*mis_wt #dm_i_s/dt contributes sRNA signal
         f3 = k_elon * mi_wt  - b_m * m_wt - k_on * s_wt * m_wt + k_off * ms_wt #dm/dt contributes mRNA signal
         f4 = k_elon_prime * mis_wt + k_on * s_wt * m_wt - b_ms * ms_wt - b_e * ms_wt - k_off * ms_wt #dms/dt contributes mRNA signal and sRNA signal
         f5 = m_wt * kx_wt + (kx_s_p*kx_wt) * ms_wt  - p_wt * b_p #dp/dt protein signal

         f6 = a_s_rne - b_s_rne * s_rne - k_on * s_rne * (m_rne) - k_on_nuc * s_rne*mi_rne+ k_off * (ms_rne) + k_off_nuc*mis_rne#ds/dt contributes sRNA signal
         f7 = k_init_rne*DNA - k_elon*mi_rne - k_on_nuc * s_rne * mi_rne + k_off_nuc*mis_rne#dm_i/dt 
         f8 = k_on_nuc * s_rne * mi_rne - k_elon_prime * mis_rne - (b_ms) * mis_rne - k_off_nuc*mis_rne #dm_i_s/dt contributes sRNA signal
         f9 = k_elon * mi_rne  - b_m * m_rne - k_on * s_rne * m_rne + k_off * ms_rne #dm/dt contributes mRNA signal
         f10 = k_elon_prime * mis_rne + k_on * s_rne * m_rne - b_ms * ms_rne - k_off * ms_rne #dms/dt contributes mRNA signal and sRNA signal
         f11 = m_rne * kx_rne + (kx_s_p*kx_rne) * ms_rne  - p_rne * b_p #dp/dt protein signal

         #print(a_m(a_m_wt,kA,c,s_wt,ms_wt))

         return [f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11];
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
    return (ys[:,0],ys[:,1],ys[:,2],ys[:,3],ys[:,4],ys[:,5],ys[:,6],ys[:,7],ys[:,8],ys[:,9],ys[:,10],ys[:,11], t)



#In the original code, they pass x, y, yerr in the place of temp1-temp3
#We don't really care about these 3 parameters because our x, y and yerr are global variables
def lnlike(args, arg1, arg2):
    global k_on , kx_s_p, b_e, kx_wt, kx_rne, b_ms, k_elon, k_off
    k_on = args[0] 
    k_off = args[1]
    kx_s_p = args[2] #GUESS
    b_e = args[3]
    b_ms = args[4]

    
    kx_wt = arg1
    kx_rne = arg2
    
    s_wt = []
    mi_wt = []
    mis_wt = []
    m_wt = []
    ms_wt = []
    p_wt = []
    
    s_rne = []
    mi_rne = []
    mis_rne = []
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
        mi_wt.append(0)
        mis_wt.append(0)
        m_wt.append(0)
        ms_wt.append(0)
        p_wt.append(0)
        s_rne.append(0)
        mi_rne.append(0)
        mis_rne.append(0)
        m_rne.append(0)
        ms_rne.append(0)
        p_rne.append(0)
        t.append(0)

        s_wt[i], mi_wt[i], mis_wt[i], m_wt[i], ms_wt[i], p_wt[i], s_rne[i], mi_rne[i], mis_rne[i], m_rne[i], ms_rne[i], p_rne[i], t[i] = eq(initial_cond[i])
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
            distance_sRNA[r] = (s_wt[r][x_RNA[i]] + ms_wt[r][x_RNA[i]] + mis_wt[r][x_RNA[i]] - yData[r, 2, real_RNA_time[i]])
            distance_protein[r] = (p_wt[r][x_protein_wt[e]] - yData[r, 0, real_protein_time_wt[e]])
            inv_s2_mRNA[r] = 1/(error_mRNA**2)
            inv_s2_protein[r] = 1/(error_protein**2)
            inv_s2_sRNA[r] = 1/(error_sRNA**2)
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
            distance_sRNA[r] = (s_wt[r][x_RNA[i]] + ms_wt[r][x_RNA[i]] + mis_wt[r][x_RNA[i]] - yData[r, 2, real_RNA_time[i]])
            
            inv_s2_mRNA[r] = 1/(error_mRNA**2)
            inv_s2_sRNA[r] = 1/(error_sRNA**2)
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
            distance_sRNA[r] = (s_rne[r][x_RNA[i]] + ms_rne[r][x_RNA[i]] + mis_rne[r][x_RNA[i]] - yData[r, 5, real_RNA_time[i]])
            distance_protein[r] = (p_rne[r][x_protein_rne[e]] - yData[r, 3, real_protein_time_rne[e]])
            # print(m_rne[r][x_RNA[i]])
            # print(yData[r, 4, real_RNA_time[i]])

            inv_s2_mRNA[r] = 1/(error_mRNA**2)
            inv_s2_protein[r] = 1/(error_protein**2)
            inv_s2_sRNA[r] = 1/(error_sRNA**2)
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
            distance_sRNA[r] = (s_rne[r][x_RNA[i]] + ms_rne[r][x_RNA[i]] + mis_rne[r][x_RNA[i]] - yData[r, 5, real_RNA_time[i]])
            inv_s2_mRNA[r] = 1/(error_mRNA**2)
            inv_s2_sRNA[r] = 1/(error_sRNA**2)
            distance[r] = ((distance_mRNA[r]**2)*inv_s2_mRNA[r] - np.log(inv_s2_mRNA[r])) + ((distance_sRNA[r]**2)*inv_s2_sRNA[r] - np.log(inv_s2_sRNA[r]))

        likeAll[3] += -0.5*(sum(distance))
    return sum(likeAll)


#lnprior. Function of likelihood od parameters
def lnprior(args):
    #Set conditions for args

    k_on = args[0] 
    k_off = args[1]
    kx_s_p = args[2] #GUESS
    b_e = args[3]
    b_ms = args[4]

    kx_wt = args[5]
    kx_rne = args[6]
    

        
    
    if not (1E-9 <= k_on <= 1E-2 and 1E-5 <= k_off <= 1.0 and .01 <= kx_s_p <= 1.10 and 1E-6 <= b_e <= 1.0 and kx_wt_ml - 1.01*sigma_kx_wt <= kx_wt <= kx_wt_ml + 1.01*sigma_kx_wt and kx_rne_ml - 1.01*sigma_kx_rne <= kx_rne <= kx_rne_ml + 1.01*sigma_kx_rne and b_m_ml - 1.01*sigma_b_m <= b_m <= b_m_ml + 1.01*sigma_b_m and b_m <= b_ms <= 0.8):
        return -np.inf
    
    def GaussError(mu,sigma,curr_val):
        return (np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5*(curr_val-mu)**2/sigma**2)


    return np.sum([GaussError(kx_wt_ml,sigma_kx_wt,kx_wt),GaussError(kx_rne_ml,sigma_kx_rne,kx_rne),GaussError(b_m_ml,sigma_b_m,b_m)])

#In the original code, they pass x, y, yerr in the place of temp1-temp3
#We don't really care about these 3 parameters because our x, y and yerr are global variables
def lnprob(args):
    # print(args)
    lp = lnprior(args)
    if not np.isfinite(lp):
        return -np.inf, -np.inf, -np.inf
    lnl = lnlike(args, args[5], args[6])
    if not np.isfinite(lnl):
        return -np.inf, -np.inf, -np.inf
    return lp + lnl, lnl, lp + lnl


def fixPos():
    for i in range(len(pos)):
        pos[i][0] = np.clip(pos[i][0], 1E-9, 1E-2) #k_on
        pos[i][1] = np.clip(pos[i][1], 1E-5, 1.0) #k_off
        pos[i][2] = np.clip(pos[i][2], .01, 1.10) #Kx_s_p
        pos[i][3] = np.clip(pos[i][3], 1E-6, 1.0)#b_e
        
        pos[i][5] = np.clip(pos[i][5], kx_wt_ml - sigma_kx_wt, kx_wt_ml + sigma_kx_wt)#kx_wt
        pos[i][6] = np.clip(pos[i][6], kx_rne_ml - sigma_kx_rne, kx_rne_ml + sigma_kx_rne)#kx_rne
        
        pos[i][4] = np.clip(pos[i][4], b_m_ml,.8)#b_m_s
        



def graph():
    global result, pos
    ndim, nwalkers = 7, 50
#    nll = lambda *args: -lnlike(*args)
#    bnds = ((1E-8, 1E-1), (.2, .99), (.4*a_m_wt_minus, a_m_wt_minus), (.4*a_m_rne_minus, a_m_rne_minus), (5E-4, 0.01), (.0012, .0016), (.0007, .0011))
#    result = op.minimize(nll, [kA, kx_s_p, a_m_wt, a_m_rne, b_e, b_s_wt, b_s_rne], args = (kx_wt_ml, kx_rne_ml), method ='TNC', bounds = bnds)
#    kA_ml, kx_s_p_ml, a_m_wt_ml, a_m_rne_ml, b_e_ml, b_s_wt_ml, b_s_rne_ml = result["x"]
    k_on_ml = 1.591E-4
    k_off_ml = 1.591E-6
    kx_s_p_ml = 0.411
    b_e_ml = 1.243E-2
    b_ms_ml = b_m_ml + sigma_b_m
    pos = [[k_on_ml, k_off_ml, kx_s_p_ml, b_e_ml, b_ms_ml, kx_wt_ml, kx_rne_ml]
            + np.array([1E-2*np.random.randn(), 1E-2*np.random.randn(), .1*np.random.randn(), 0.01*np.random.randn(), .01*np.random.randn(), 0.1*np.random.randn(), 0.1*np.random.randn()])
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
