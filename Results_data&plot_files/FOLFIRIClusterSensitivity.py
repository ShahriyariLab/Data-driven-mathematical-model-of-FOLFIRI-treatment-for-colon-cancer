'''
FOLFIRIClusterSensitivity: script for analysis of  parameter sensitivity  of 
                   Data driven mathematical model of FOLFIRI treatment for colon cancer*

Author: Arkadz Kirshtein, https://sites.google.com/site/akirshtein/
(c) Shahriyari Lab https://sites.google.com/site/leilishahriyari/

*part of https://github.com/ShahriyariLab/Data-driven-mathematical-model-of-FOLFIRI-treatment-for-colon-cancer
 If using this or related code please cite 
  Budithi, A.; Su, S.; Kirshtein,A.; Shahriyari L.
    Data driven mathematical model of FOLFIRI treatment for colon cancer. //Cancers.
     (Manuscript submitted for publication).
'''
from qspmodel import *
from folfiri_qspmodel import *
import pandas as pd
import csv
import os
import scipy as sp

if not os.path.exists('Results/'):
    os.makedirs('Results/Sensitivity/')
elif not os.path.exists('Data/Sensitivity/'):
    os.makedirs('Results/Sensitivity/')
        

# some global parameters
lmod=[0, 1, 2, 3, 4, 5, 6, 8, 9] #indices of variables in cell data (excluding M0)
clusters=5 #number of clusters

number_of_cycles=6
cycle_time=14
injection_time=2

# time discretization with more points during the infusion than the rest of the cycle
# i.e. 513 points during the infusion and 128 point in the rest of the cycle
T_treatment=number_of_cycles*cycle_time
start_treatment = 10
t=np.linspace(0,start_treatment,(2*start_treatment)+1)
cycle_t=np.concatenate((np.linspace(0,injection_time,1025),np.linspace(injection_time,cycle_time,257)[1:]))

for cycle in range(number_of_cycles):
    t=np.concatenate((t,(start_treatment+cycle_t+cycle*cycle_time)[1:]))

# adds extra time after the treatment discretizing it with (2*time_in_days+1) points
obs_time = 365*3
T=start_treatment+T_treatment+obs_time
t=np.concatenate((t,np.linspace(start_treatment+T_treatment,T,(2*obs_time)+1)[1:]))

# T=3500
# t=np.linspace(0, T, 35001)

# pre-defined parameters for treatment deltas, cell_steady, cell_extreme, drug_ratio, drug_eff, 5FU_killrate
drug_deltas=[71.3, 1.45, 2.56]
cell_steady_ind=[5, 6, 13]
cell_extreme_ind=[5, 6, 12]
drug_levels=np.array([770, 300, 725])/cycle_time
drug_min=np.array([598, 208, 75])/cycle_time
drug_ratio=drug_levels/drug_min
drug_eff=[0.2, 0.4, 0.05]
FU_killrate=100

# create parameters for step function 
infusion_indices=np.concatenate([15*np.ones(number_of_cycles), 16*np.ones(number_of_cycles), 17*np.ones(number_of_cycles)]).astype(int)
FUinfusion=np.array([2/24, 2])
IRinfusion=np.array([0, 1/24])
LVinfusion=np.array([1/24, 2/24])
infusion_intervals=np.concatenate([[start_treatment+n*cycle_time+FUinfusion for n in range(number_of_cycles)],
                                    [start_treatment+n*cycle_time+IRinfusion for n in range(number_of_cycles)],
                                    [start_treatment+n*cycle_time+LVinfusion for n in range(number_of_cycles)]])
infusion_values=np.concatenate([((cycle_time/(FUinfusion[1]-FUinfusion[0]))*drug_deltas[0]*np.ones(number_of_cycles) if drug_ratio[0]>0 else np.zeros(number_of_cycles)),
                                ((cycle_time/(IRinfusion[1]-IRinfusion[0]))*drug_deltas[1]*np.ones(number_of_cycles) if drug_ratio[1]>0 else np.zeros(number_of_cycles)),
                                ((cycle_time/(LVinfusion[1]-LVinfusion[0]))*drug_deltas[2]*np.ones(number_of_cycles) if drug_ratio[2]>0 else np.zeros(number_of_cycles))])

#other parameters
nvar=Colon_5fu_Functions().nvar # number of variables
nparam=Colon_5fu_Functions().nparam # number of parameters

# infusion function
r= step_vector(nvar, indices=infusion_indices, intervals=infusion_intervals, values=infusion_values)


# define cell data for non-treatment parameter derivation
# rates acquired from bio research [delTn delTh delTc delTr delD delM lamgCmax lamgCmin delmu1 delmu2 delmu3 delH delIg delGb]
globvals=np.array([9.4951e-4, 0.231, 0.406, 0.231, 0.277, 0.02, 7.5e-3, 6.7152e-04, 1.07, 4.62, 58.7, 33.27, 499])
# average values of each variable across all patients [Tc mu1 Ig Gb]
meanvals=np.array([2920.305826, 132.316278, 6.603464, 20017.990343])
extremevals=np.array([2.67309712e+04, 1.23113649e+04, 1.91073391e+04, 7.21023854e+03,
                              3.41730382e+03, 4.32746936e+03, 2.31601361e+04, 3.24717521e+05,
                              1.62358760e+05, 1.05768591e+03, 1.39828734e+04, 2.46968870e+04,
                              1.02214368e+02, 1.63407223e+05])
celldata=[globvals, extremevals, meanvals]


clustercells=pd.read_csv('input/Large_Tumor_cell_data.csv').to_numpy()

#Computations of Sensitivity
print('Starting Sensitivity computations')

IC= np.ones((clusters,nvar))
IC[:,-3:]=0;


#IC=np.array(pd.read_csv('input/Initial_Conditions.csv'))
#IC=np.block([IC,np.zeros((clusters,3))])

# coefficients for variable sensitivity
lambda0=np.zeros((9,2))
lambda0[7,0]=1 # just cancer

for cluster in range(clusters):
    print('Starting computations for cluster '+str(cluster+1)+' of 5') 
    filename='V63-cluster-'+str(cluster+1)+'-of-5-results-'
    
    lambda0[0:9,1]=np.copy(clustercells[cluster,lmod]/np.sum(clustercells[cluster,lmod])) # all cells   

    QSP0=QSP.from_data([clustercells[cluster]]+celldata)
    # QSP0=QSP.from_cell_data(clustercells[cluster])
    qspcore=Colon_5fu_Functions(parameters=QSP0.par)
    #QSP_=QSP(parameters=qspcore.parameters_from_assumptions([drug_deltas, (clustercells[cluster,cell_steady_ind]/extremevals[cell_extreme_ind]), drug_ratio, drug_eff, FU_killrate]), qspcore=qspcore)
    QSP_=QSP.from_data([drug_deltas, (clustercells[cluster,cell_steady_ind]/extremevals[cell_extreme_ind]), drug_ratio, drug_eff, FU_killrate],  qspcore=qspcore)
    
    print(' Parameters set. Computing the solution')
    u, _ = QSP_.solve_ode(t, IC[cluster], 'given', inhomogeneity=r, jumps=True)
    
    
    print(' Computing local sensitivity')
    S=QSP_.Sensitivity(method="time-full", t=t, IC=IC[cluster], inhomogeneity=r, variables=np.arange(9), jumps=True)
    c=csv.writer(open('Results/Sensitivity/'+filename+'sensitivity_local.csv',"w"))
    c.writerows(np.dot(np.mean(S, axis=0),lambda0) )
    del c
    c=csv.writer(open('Results/Sensitivity/'+filename+'sensitivity_local_relative.csv',"w"))
    c.writerows(np.dot(np.mean(S*QSP_.par[np.newaxis,:,np.newaxis]/u[:,np.newaxis,np.arange(9)], axis=0),lambda0) )
    del c
    
    
