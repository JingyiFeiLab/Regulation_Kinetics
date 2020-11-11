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
from scipy.interpolate import griddata


xData = [0, 60, 180, 360, 720, 1080, 1440]

result = []


numReplicates = 3
timePoints = 8
numMeasurable = 6
#0 is p_wt, 1 is m_wt, 2 is s_wt, 3 p_rne, 4 m_rne, 5 s_rine
yData = np.zeros(shape=(numReplicates, numMeasurable, timePoints))

yData[1, 0, :] = [347513.2853, 367107.5708, 367229.4779, 526899.8252, 765383.2836, 893353.6218, 1194864.312, 1611808.456]
yData[1, 1, :] = [126759.8746, 174201.1888, 861566.6491, 829035.7939, 1255372.267, 1307704.871, 1447041.427, 1047774.86]
yData[1, 2, :] = [1846825.091, 2026286.474, 1517518.905, 1095534.472, 1302512.952, 1455675.663, 1423814.456, 1391623.404]
yData[0, 3, :] = [30740.06404, 7611.392039, 28438.54689, 23048.9751, 167085.4364, 296431.9934, 591670.9207, 712535.0848]
yData[0, 4, :] = [22969.08289, 13337.21322, 162812.0361, 435224.8164, 905086.9004, 1067163.977, 1082526.701, 1055862.227]
yData[0, 5, :] = [2085525.097, 2072528.385, 2038819.034, 2209784.252, 1869579.711, 2355808.953, 1676761.391, 1123533.284]


yData[0, 0, :] = [435500.661, 412698.7449, 429219.6278, 663807.8444, 1035433.729, 2144265.792, 2486745.291, 2684783.042]
yData[0, 1, :] = [55235.83625, 123960.7343, 1164299.635, 1652964.543, 1516875.402, 1577627.268, 1183959.582, 1500131.281]
yData[0, 2, :] = [1371638.973, 1479740.334, 2300059.658, 1877851.897, 1420127.483, 1090447.296, 1096737.351, 1324216.114]
yData[1, 3, :] = [25053.85663, 27674.565, 30010.12027, 67783.35453, 159716.7395, 416788.5436, 627153.735, 323987.3492]
yData[1, 4, :] = [21402.5558, 470096.1962, 114914.2647, 502208.8933, 797456.0465, 1006640.656, 859554.7947, 666903.2247]
yData[1, 5, :] = [2056380.255, 2030448.259, 2268022.237, 2531172.333, 2419731.03, 2262055.679, 1935417.938, 831907.646]

yData[2, 0, :] = [103121.4438,119740.8521,178169.0566,173472.7545,344063.0791,520752.0862,820034.7537,1184189.091]
yData[2, 1, :] = [18532.05393,42313.70869,650846.2755,1082916.612,1025421.802,1344767.177,1548645.684,1403606.94]
yData[2, 2, :] = [1142440.97,949293.0491,1536117.023,1213753.668,664036.2236,1182033.928,1183425.272,1304065.765]
yData[2, 3, :] = [25053.85663, 27674.565, 30010.12027, 67783.35453, 159716.7395, 416788.5436, 627153.735, 323987.3492]
yData[2, 4, :] = [21402.5558, 470096.1962, 114914.2647, 502208.8933, 797456.0465, 1006640.656, 859554.7947, 666903.2247]
yData[2, 5, :] = [2056380.255, 2030448.259, 2268022.237, 2531172.333, 2419731.03, 2262055.679, 1935417.938, 831907.646]

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

k_init_wt = 0.0582535
k_init_wt_ml = 0.0582535
sigma_k_init_wt = 0.003884137549

k_init_rne = 0.0432
k_init_rne_ml = 0.0432
sigma_k_init_rne = 0.002404163056


k_elon = 0.0641
k_elon_prime = 0.0641

b_m = 0.003248
b_m_ml = 0.003248
sigma_b_m = .0002

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

kx_wt = 10.5897
kx_rne = 5.4598
kx_wt_ml = 10.5897
kx_rne_ml = 5.4598
sigma_kx_wt = 0.9456846092
sigma_kx_rne = 0.2870853532


kx_s_p = .358

b_e = 2.82E-03#GUESS
b_ms = .0025
b_ms_ml = .0025

b_nuc_wt = 10
b_nuc_rne = 10

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
         
         


         f0 = a_s_wt - b_s_wt * s_wt - k_on * s_wt * (m_wt) - k_on * s_wt*mi_wt + k_off * (ms_wt + mis_wt) #ds/dt contributes sRNA signal
         f1 = k_init_wt*DNA - k_elon*mi_wt - k_on * s_wt * mi_wt + c*(k_off*mis_wt)#dm_i/dt 
         f2 = k_on * s_wt * mi_wt - k_off*mis_wt #dm_i_s/dt contributes sRNA signal
         f3 = k_elon * mi_wt  - b_m * m_wt - k_on * s_wt * m_wt + k_off * ms_wt #dm/dt contributes mRNA signal
         f4 = k_on * s_wt * m_wt - b_ms * ms_wt - b_e * ms_wt - k_off * ms_wt #dms/dt contributes mRNA signal and sRNA signal
         f5 = m_wt * kx_wt + (kx_s_p*kx_wt) * ms_wt  - p_wt * b_p #dp/dt protein signal

         f6 = a_s_rne - b_s_rne * s_rne - k_on * s_rne * (m_rne) - k_on * s_rne*mi_rne + k_off * (ms_rne + mis_rne) #ds/dt contributes sRNA signal
         f7 = k_init_rne*DNA - k_elon*mi_rne - k_on * s_rne * mi_rne + d*(k_off*mis_rne)#dm_i/dt 
         f8 = k_on * s_rne * mi_rne - k_off*mis_rne #dm_i_s/dt contributes sRNA signal
         f9 = k_elon * mi_rne  - b_m * m_rne - k_on * s_rne * m_rne + k_off * ms_rne #dm/dt contributes mRNA signal
         f10 = k_on * s_rne * m_rne - b_ms * ms_rne - k_off * ms_rne #dms/dt contributes mRNA signal and sRNA signal
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
def lnlike(args, arg1, arg2, arg3, arg4):
    global k_on , kx_s_p, b_e, kx_wt, kx_rne, b_ms, k_elon, k_off, c, d
    k_on = args[0] 
    k_off = args[1]
    kx_s_p = args[2] #GUESS
    b_e = args[3]
    b_ms = args[4]

    kx_wt = arg1
    kx_rne = arg2
    c = arg3
    d = arg4
    
    
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
    c = args[7]
    d = args[8]
    

        
    
    if not (1E-6 <= k_on <= 1E-2 and 1E-3 <= k_off <= 10.0 and .01 <= kx_s_p <= 1.10 and .01 <= c <= 1.00 and .01 <= d <= 1.00 and 1E-6 <= b_e <= 1.0 and kx_wt_ml - 1.01*sigma_kx_wt <= kx_wt <= kx_wt_ml + 1.01*sigma_kx_wt and kx_rne_ml - 1.01*sigma_kx_rne <= kx_rne <= kx_rne_ml + 1.01*sigma_kx_rne and b_m <= b_ms <= 0.8):
        return -np.inf
    
    def GaussError(mu,sigma,curr_val):
        return (np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5*(curr_val-mu)**2/sigma**2)


    return np.sum([GaussError(kx_wt_ml,sigma_kx_wt,kx_wt),GaussError(kx_rne_ml,sigma_kx_rne,kx_rne)])

#In the original code, they pass x, y, yerr in the place of temp1-temp3
#We don't really care about these 3 parameters because our x, y and yerr are global variables
def lnprob(args):
    # print(args)
    lp = lnprior(args)
    if not np.isfinite(lp):
        return -np.inf, -np.inf, -np.inf
    lnl = lnlike(args, args[5], args[6], args[7], args[8])
    if not np.isfinite(lnl):
        return -np.inf, -np.inf, -np.inf
    return lp + lnl, lnl, lp + lnl


def fixPos():
    for i in range(len(pos)):
        pos[i][0] = np.clip(pos[i][0], 1E-6, 1E-2) #k_on
        pos[i][1] = np.clip(pos[i][1], 1E-3, 10.0) #k_off
        pos[i][2] = np.clip(pos[i][2], .01, 1.10) #Kx_s_p
        pos[i][3] = np.clip(pos[i][3], 1E-6, 1.0)#b_e
        
        pos[i][5] = np.clip(pos[i][5], kx_wt_ml - sigma_kx_wt, kx_wt_ml + sigma_kx_wt)#kx_wt
        pos[i][6] = np.clip(pos[i][6], kx_rne_ml - sigma_kx_rne, kx_rne_ml + sigma_kx_rne)#kx_rne
        pos[i][7] = np.clip(pos[i][7], 0.01, 1.0)#b_e
        pos[i][8] = np.clip(pos[i][8], 0.01, 1.0)#b_e
        
        pos[i][4] = np.clip(pos[i][4], b_m_ml,.8)#b_m_s



def graph():
    global result, pos
    ndim, nwalkers = 9, 50
#    nll = lambda *args: -lnlike(*args)
#    bnds = ((1E-8, 1E-1), (.2, .99), (.4*a_m_wt_minus, a_m_wt_minus), (.4*a_m_rne_minus, a_m_rne_minus), (5E-4, 0.01), (.0012, .0016), (.0007, .0011))
#    result = op.minimize(nll, [kA, kx_s_p, a_m_wt, a_m_rne, b_e, b_s_wt, b_s_rne], args = (kx_wt_ml, kx_rne_ml), method ='TNC', bounds = bnds)
#    kA_ml, kx_s_p_ml, a_m_wt_ml, a_m_rne_ml, b_e_ml, b_s_wt_ml, b_s_rne_ml = result["x"]
    k_on_ml = 1.591E-4
    k_off_ml = 1.591E-1
    kx_s_p_ml = 0.411
    b_e_ml = 1.243E-2
    b_ms_ml = b_m_ml + sigma_b_m
    c_ml = 0.5
    d_ml = 0.5
    pos = [[k_on_ml, k_off_ml, kx_s_p_ml, b_e_ml, b_ms_ml, kx_wt_ml, kx_rne_ml, c_ml, d_ml]
            + np.array([1E-2*np.random.randn(), 1E-2*np.random.randn(), .1*np.random.randn(), 0.01*np.random.randn(), 0.01*np.random.randn(), 0.1*np.random.randn(), 0.1*np.random.randn(), 0.1*np.random.randn(), 0.1*np.random.randn()])
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


graph()
# postProcess()
# heatMap()
