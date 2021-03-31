#!/usr/bin/env python
# coding: utf-8

# In[ ]:



from qspmodel import *
import pandas as pd
import csv
import os
import scipy as sp

class Colon_5fu_Functions(object):
    def __init__(self,parameters=None,inject=True):
        self.nparam=19
        self.nvar=17
        self.variable_names=['Naive T-cells', 'helper T-cells', 'cytotoxic cells', 'Treg-cells', 'Naive Dendritic cells', 'Dendritic cells', 'Macrophages', 
                                     'Cancer', 'Necrotic cells', '\mu_1', '\mu_2', 'HMGB1', 'IFN-\gamma', 'TNF-\beta', '5-FU', 'IR', 'LV']

        self.parameter_names=['\lambda_{ThD5fu}', '\lambda_{5fuDTc}', '\delta_{C5fuIg}', '\delta_{5fuM}', '\alpha_{5fu}', '\delta_{5fuD}', '\delta_{5fu}',                                 '\delta_{CFU}', '\delta_{CFULV}', '\alpha_{IRC}', '\alpha_{IRTr}', '\delta_{CIR}', '\delta_{TrIR}', '\delta_{IR}', '\delta_{LV}',                                 '\alpha_{LV}', 'A_{injFU}', 'A_{injIR}', 'A_{injLV}']
        if parameters is None:
            self.par0 = [2.13992019e+00, 3.91714541e-02, 1.18709540e-02, 3.25025264e+00, 6.00530110e-01, 8.28979547e-01, 5.11203137e-02, 6.41203692e-02, 2.61342952e-01, 1.01287620e-01, 6.03885790e-04, 5.45360147e-04, 2.44819003e-02, 5.33231239e-03, 3.64525256e-03, 1.50000000e+00, 9.51544836e-02, 9.68665712e-01, 6.17980430e-03, 3.68555643e+00, 1.88102188e-01, 7.46341386e-01, 5.66450410e+01, 1.46026677e+00, 2.86891398e-01, 1.59945661e-01, 1.47855221e-01, 8.89798679e+00, 1.98429703e+01, 4.52904287e+00, 4.14967297e+02, 8.40327027e+01, 9.49510000e-04, 4.29109010e-01, 1.53085359e+00, 2.31000000e-01, 7.54191593e-01, 2.69059116e+00, 4.06000000e-01, 7.13220230e-01, 2.31000000e-01, 3.35771449e-02, 5.20534274e-02, 2.77000000e-01, 2.00000000e-02, 9.15995362e-04, 6.59506923e-04, 2.12702213e-03, 7.86258064e-04, 6.73317372e-03, 1.07000000e+00, 4.62000000e+00, 5.87000000e+01, 3.32700000e+01, 4.99000000e+02, 1.50058799e+00, 1.02638785e+00, 1.78030065e+0, 2.13992019e+00, 3.91714541e-02, 1.18709540e-02, 3.25025264e+00, 6.00530110e-01]
        else:
            self.par0 = parameters
        self.inject=inject
    def __call__(self,t,x,par):
        # ODE right-hand side
        x=x*(x>0) #Q
        return np.array([self.par0[55] - self.par0[32]*x[0] - 
                       self.par0[61]*x[0]*(self.par0[5]*x[1] + 
                       self.par0[6]*x[10] + self.par0[7]*x[13]) - 
                       self.par0[59]*x[0]*(self.par0[0]*x[5] + 
                       self.par0[1]*x[6] + self.par0[2]*x[9] + 
                       par[0]*x[5]*x[14]) - self.par0[60]*x[0]*
                       (self.par0[3]*x[1] + self.par0[4]*x[5] + 
                       par[1]*x[5]*x[14]), 
                       (-x[1])*(self.par0[35] + self.par0[34]*x[3] + 
                       self.par0[33]*x[10]) + x[0]*(self.par0[0]*x[5] + 
                       self.par0[1]*x[6] + self.par0[2]*x[9] + 
                       par[0]*x[5]*x[14]), 
                       (-x[2])*(self.par0[38] + self.par0[37]*x[3] + 
                       self.par0[36]*x[10]) + x[0]*(self.par0[3]*x[1] + 
                       self.par0[4]*x[5] + par[1]*x[5]*x[14]), 
                       x[0]*(self.par0[5]*x[1] + self.par0[6]*x[10] + 
                       self.par0[7]*x[13]) - x[3]*(self.par0[40] + 
                       self.par0[39]*x[9] + par[12]*x[15]), 
                       self.par0[56] - self.par0[62]*x[4]*(self.par0[9]*x[7] + 
                       self.par0[8]*x[11]) - x[4]*(self.par0[43] + 
                       self.par0[41]*x[11]), 
                       x[4]*(self.par0[9]*x[7] + self.par0[8]*x[11]) - 
                       x[5]*(self.par0[43] + self.par0[42]*x[7] + 
                       self.par0[41]*x[11]), (-self.par0[44])*x[6] + 
                       (self.par0[57] - x[6])*(self.par0[12]*x[1] + 
                       self.par0[10]*x[10] + self.par0[11]*x[12]), 
                       x[7]*(1 - x[7]/self.par0[58])*(self.par0[13] + 
                       self.par0[14]*x[9]) - x[7]*(self.par0[48] + 
                       self.par0[47]*x[2] + self.par0[46]*x[12] + 
                       self.par0[45]*x[13] + par[7]*x[14] - par[3]*x[6]*x[14] + 
                       par[2]*x[12]*x[14] + par[11]*x[15] + par[8]*x[14]*x[16]), 
                       (-self.par0[49])*x[8] + self.par0[15]*x[7]*
                       (self.par0[48] + self.par0[47]*x[2] + 
                       self.par0[46]*x[12] + self.par0[45]*x[13] + 
                       par[7]*x[14] - par[3]*x[6]*x[14] + par[2]*x[12]*x[14] + 
                       par[11]*x[15] + par[8]*x[14]*x[16]), 
                       self.par0[16]*x[1] + self.par0[18]*x[5] + 
                       self.par0[17]*x[6] - self.par0[50]*x[9], 
                       self.par0[21]*x[3] + self.par0[20]*x[5] + 
                       self.par0[19]*x[6] - self.par0[51]*x[10], 
                       self.par0[24]*x[1] + self.par0[25]*x[2] + 
                       self.par0[26]*x[3] + self.par0[23]*x[6] + 
                       self.par0[22]*x[8] - self.par0[52]*x[11], 
                       self.par0[27]*x[1] + self.par0[28]*x[2] + 
                       self.par0[29]*x[6] - self.par0[53]*x[12], 
                       self.par0[31]*x[3] + self.par0[30]*x[6] - 
                       self.par0[54]*x[13], par[16] - par[6]*x[14] - 
                       par[5]*x[5]*x[14] - par[4]*x[7]*(par[7]*x[14] - 
                       par[3]*x[6]*x[14] + par[2]*x[12]*x[14] + 
                       par[8]*x[14]*x[16]), par[17] - par[13]*x[15] - 
                       par[10]*par[12]*x[3]*x[15] - par[9]*par[11]*x[7]*x[15], 
                       par[18] - par[14]*x[16] - par[8]*par[15]*x[7]*x[14]*x[16]])
    def Ju(self,t,x,par):
        # Jacobian with respect to variables
        x=x*(x>0)
        return np.array([[-self.par0[32] - self.par0[61]*(self.par0[5]*x[1] + 
                          self.par0[6]*x[10] + self.par0[7]*x[13]) - 
                          self.par0[59]*(self.par0[0]*x[5] + self.par0[1]*x[6] + 
                          self.par0[2]*x[9] + par[0]*x[5]*x[14]) - 
                          self.par0[60]*(self.par0[3]*x[1] + self.par0[4]*x[5] + 
                          par[1]*x[5]*x[14]), (-self.par0[3])*self.par0[60]*
                          x[0] - self.par0[5]*self.par0[61]*x[0], 0, 0, 0, 
                          (-self.par0[59])*x[0]*(self.par0[0] + par[0]*x[14]) - 
                          self.par0[60]*x[0]*(self.par0[4] + par[1]*x[14]), 
                          (-self.par0[1])*self.par0[59]*x[0], 0, 0, 
                          (-self.par0[2])*self.par0[59]*x[0], 
                          (-self.par0[6])*self.par0[61]*x[0], 0, 0, 
                          (-self.par0[7])*self.par0[61]*x[0], 
                          (-self.par0[59])*par[0]*x[0]*x[5] - self.par0[60]*par[1]*
                          x[0]*x[5], 0, 0], [self.par0[0]*x[5] + 
                          self.par0[1]*x[6] + self.par0[2]*x[9] + 
                          par[0]*x[5]*x[14], -self.par0[35] - self.par0[34]*x[3] - 
                          self.par0[33]*x[10], 0, (-self.par0[34])*x[1], 0, 
                          x[0]*(self.par0[0] + par[0]*x[14]), self.par0[1]*x[0], 0, 
                          0, self.par0[2]*x[0], (-self.par0[33])*x[1], 0, 0, 0, 
                          par[0]*x[0]*x[5], 0, 0], 
                          [self.par0[3]*x[1] + self.par0[4]*x[5] + par[1]*x[5]*x[14], 
                          self.par0[3]*x[0], -self.par0[38] - self.par0[37]*x[3] - 
                          self.par0[36]*x[10], (-self.par0[37])*x[2], 0, 
                          x[0]*(self.par0[4] + par[1]*x[14]), 0, 0, 0, 0, 
                          (-self.par0[36])*x[2], 0, 0, 0, par[1]*x[0]*x[5], 0, 0], 
                          [self.par0[5]*x[1] + self.par0[6]*x[10] + 
                          self.par0[7]*x[13], self.par0[5]*x[0], 0, 
                          -self.par0[40] - self.par0[39]*x[9] - par[12]*x[15], 0, 0, 
                          0, 0, 0, (-self.par0[39])*x[3], self.par0[6]*x[0], 0, 0, 
                          self.par0[7]*x[0], 0, (-par[12])*x[3], 0], 
                          [0, 0, 0, 0, -self.par0[43] - self.par0[41]*x[11] - 
                          self.par0[62]*(self.par0[9]*x[7] + self.par0[8]*x[11]), 
                          0, 0, (-self.par0[9])*self.par0[62]*x[4], 0, 0, 0, 
                          (-self.par0[41])*x[4] - self.par0[8]*self.par0[62]*x[4], 
                          0, 0, 0, 0, 0], [0, 0, 0, 0, self.par0[9]*x[7] + 
                          self.par0[8]*x[11], -self.par0[43] - 
                          self.par0[42]*x[7] - self.par0[41]*x[11], 0, 
                          self.par0[9]*x[4] - self.par0[42]*x[5], 0, 0, 0, 
                          self.par0[8]*x[4] - self.par0[41]*x[5], 0, 0, 0, 0, 0], 
                          [0, self.par0[12]*(self.par0[57] - x[6]), 0, 0, 0, 0, 
                          -self.par0[44] - self.par0[12]*x[1] - 
                          self.par0[10]*x[10] - self.par0[11]*x[12], 0, 0, 0, 
                          self.par0[10]*(self.par0[57] - x[6]), 0, 
                          self.par0[11]*(self.par0[57] - x[6]), 0, 0, 0, 0], 
                          [0, 0, (-self.par0[47])*x[7], 0, 0, 0, par[3]*x[7]*x[14], 
                          -self.par0[48] - self.par0[47]*x[2] - 
                          (x[7]*(self.par0[13] + self.par0[14]*x[9]))/
                          self.par0[58] + (1 - x[7]/self.par0[58])*
                          (self.par0[13] + self.par0[14]*x[9]) - 
                          self.par0[46]*x[12] - self.par0[45]*x[13] - 
                          par[7]*x[14] + par[3]*x[6]*x[14] - par[2]*x[12]*x[14] - 
                          par[11]*x[15] - par[8]*x[14]*x[16], 0, 
                          self.par0[14]*x[7]*(1 - x[7]/self.par0[58]), 0, 0, 
                          (-x[7])*(self.par0[46] + par[2]*x[14]), 
                          (-self.par0[45])*x[7], (-x[7])*(par[7] - par[3]*x[6] + 
                          par[2]*x[12] + par[8]*x[16]), (-par[11])*x[7], 
                          (-par[8])*x[7]*x[14]], [0, 0, self.par0[15]*self.par0[47]*
                          x[7], 0, 0, 0, (-self.par0[15])*par[3]*x[7]*x[14], 
                          self.par0[15]*(self.par0[48] + self.par0[47]*x[2] + 
                          self.par0[46]*x[12] + self.par0[45]*x[13] + 
                          par[7]*x[14] - par[3]*x[6]*x[14] + par[2]*x[12]*x[14] + 
                          par[11]*x[15] + par[8]*x[14]*x[16]), -self.par0[49], 0, 0, 
                          0, self.par0[15]*x[7]*(self.par0[46] + par[2]*x[14]), 
                          self.par0[15]*self.par0[45]*x[7], self.par0[15]*x[7]*
                          (par[7] - par[3]*x[6] + par[2]*x[12] + par[8]*x[16]), 
                          self.par0[15]*par[11]*x[7], self.par0[15]*par[8]*x[7]*
                          x[14]], [0, self.par0[16], 0, 0, 0, self.par0[18], 
                          self.par0[17], 0, 0, -self.par0[50], 0, 0, 0, 0, 0, 0, 0], 
                          [0, 0, 0, self.par0[21], 0, self.par0[20], self.par0[19], 
                          0, 0, 0, -self.par0[51], 0, 0, 0, 0, 0, 0], 
                          [0, self.par0[24], self.par0[25], self.par0[26], 0, 0, 
                          self.par0[23], 0, self.par0[22], 0, 0, -self.par0[52], 
                          0, 0, 0, 0, 0], [0, self.par0[27], self.par0[28], 0, 0, 0, 
                          self.par0[29], 0, 0, 0, 0, 0, -self.par0[53], 0, 0, 0, 0], 
                          [0, 0, 0, self.par0[31], 0, 0, self.par0[30], 0, 0, 0, 0, 
                          0, 0, -self.par0[54], 0, 0, 0], [0, 0, 0, 0, 0, 
                          (-par[5])*x[14], par[3]*par[4]*x[7]*x[14], 
                          (-par[4])*(par[7]*x[14] - par[3]*x[6]*x[14] + 
                          par[2]*x[12]*x[14] + par[8]*x[14]*x[16]), 0, 0, 0, 0, 
                          (-par[2])*par[4]*x[7]*x[14], 0, -par[6] - par[5]*x[5] - 
                          par[4]*x[7]*(par[7] - par[3]*x[6] + par[2]*x[12] + 
                          par[8]*x[16]), 0, (-par[4])*par[8]*x[7]*x[14]], 
                          [0, 0, 0, (-par[10])*par[12]*x[15], 0, 0, 0, 
                          (-par[9])*par[11]*x[15], 0, 0, 0, 0, 0, 0, 0, 
                          -par[13] - par[10]*par[12]*x[3] - par[9]*par[11]*x[7], 0], 
                          [0, 0, 0, 0, 0, 0, 0, (-par[8])*par[15]*x[14]*x[16], 0, 0, 0, 
                          0, 0, 0, (-par[8])*par[15]*x[7]*x[16], 0, 
                          -par[14] - par[8]*par[15]*x[7]*x[14]]] )
    def Jp(self,t,x,par):
        # Jacobian with respect to the parameters
        x=x*(x>0)
        return np.array([[(-self.par0[59])*x[0]*x[5]*x[14], x[0]*x[5]*x[14], 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                           [(-self.par0[60])*x[0]*x[5]*x[14], 0, x[0]*x[5]*x[14], 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 
                           (-x[7])*x[12]*x[14], self.par0[15]*x[7]*x[12]*x[14], 0, 0, 
                           0, 0, 0, (-par[4])*x[7]*x[12]*x[14], 0, 0], 
                           [0, 0, 0, 0, 0, 0, 0, x[6]*x[7]*x[14], (-self.par0[15])*x[6]*
                           x[7]*x[14], 0, 0, 0, 0, 0, par[4]*x[6]*x[7]*x[14], 0, 0], 
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                           (-x[7])*(par[7]*x[14] - par[3]*x[6]*x[14] + 
                           par[2]*x[12]*x[14] + par[8]*x[14]*x[16]), 0, 0], 
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-x[5])*x[14], 0, 
                           0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x[14], 0, 0], 
                           [0, 0, 0, 0, 0, 0, 0, (-x[7])*x[14], self.par0[15]*x[7]*
                           x[14], 0, 0, 0, 0, 0, (-par[4])*x[7]*x[14], 0, 0], 
                           [0, 0, 0, 0, 0, 0, 0, (-x[7])*x[14]*x[16], 
                           self.par0[15]*x[7]*x[14]*x[16], 0, 0, 0, 0, 0, 
                           (-par[4])*x[7]*x[14]*x[16], 0, (-par[15])*x[7]*x[14]*x[16]], 
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                           (-par[11])*x[7]*x[15], 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, (-par[12])*x[3]*x[15], 0], 
                           [0, 0, 0, 0, 0, 0, 0, (-x[7])*x[15], self.par0[15]*x[7]*
                           x[15], 0, 0, 0, 0, 0, 0, (-par[9])*x[7]*x[15], 0], 
                           [0, 0, 0, (-x[3])*x[15], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                           (-par[10])*x[3]*x[15], 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, -x[15], 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, -x[16]], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, (-par[8])*x[7]*x[14]*x[16]], 
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], 
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], 
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]] )
                            
    def parameters_from_assumptions(self,data):
        #Function returning parameter values based on assumptions and given reference data:
        # data is a list of arrays in the following order 
        # [deltas, cell_steady, cell_extreme, drug_ratio, drug_eff, 5FU_killrate]
        #   cell_ratio are given in order [D, M, Ig] and are equal to cell_steady/cell_extreme
        #   deltas=[\delta_{5fu}, \delta_{IR}, \delta_{LV}] - are the drug decay rates
        #   drug_ratio - !daily! average injection rates divided by minimal injection rate
        #   drug_eff - drug efficiencies [5FU_eff, IR_eff, LV_eff]
        #   5FU_killrate - coefficient relating \delta_{CIR} with \delta_{C}
        
        
        deltas, cell_ratio, drug_ratio, drug_eff, FU_killrate= data
        par=np.empty(19)
        par[0]=0.5*self.par0[0]*drug_ratio[0]/deltas[0] #\lambda_{ThD5fu}
        par[1]=0.5*self.par0[4]*drug_ratio[0]/deltas[0] #\lambda_{TcD5fu}
        par[3]=0.5*self.par0[48]*cell_ratio[1]*drug_ratio[0]/deltas[0] #\delta_{5fuM}
        par[5]=0.5*self.par0[48]*cell_ratio[0]*drug_ratio[0]/deltas[0] #\delta_{5fuD}
        par[6]=deltas[0] #\delta_{5fu}
        par[7]=FU_killrate*drug_ratio[0]*self.par0[48]/deltas[0] #\delta_{C5fu}
        par[2]=1.5*par[7]*cell_ratio[2] #\delta_{C5fuIg}
        par[8]=0.1*drug_ratio[2]*par[7]/deltas[2] # \delta_{C5fuLV}
        par[11]=0.5*par[7]*drug_ratio[1]*deltas[0]/(deltas[1]*drug_ratio[0]) # \delta_{CIR}
        par[12]=par[11] # \delta_{TrIR}
        par[13]=deltas[1] # \delta_{IR}
        par[14]=deltas[2] # \delta_{LV}
        par[4]=(( (drug_eff[0]/(1-drug_eff[0]))*deltas[0] - par[5] )/
                ( par[7]+par[2]-par[3]+par[8]*(1-drug_eff[2]) )) #\alpha_{5fu}
        par[9]=drug_eff[1]*deltas[1]/((1-drug_eff[1])*(par[11]+par[12])) # \alpha_{IRC}
        par[10]=par[9] # \alpha_{IRTr}
        par[15]=drug_eff[2]*deltas[2]/((1-drug_eff[2])*par[8]*(1-drug_eff[0])) # \alpha_{LV}
        if self.inject:
            par[16]=deltas[0] #A_{injfu}
            par[17]=deltas[1] #A_{injir}
            par[18]=deltas[2] #A_{injlv}
        else:
            par[16]=0
            par[17]=0
            par[18]=0
            
        return par

