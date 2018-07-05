# Orbital fitting program
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import matplotlib.cm as cm
import lmfit as lm
import os as os
import param as op
import sys
import time as ti
plt.ion() # Activate update





def orbital(flname):
    """
    Function that performs orbital fitting to file
    """

    if os.path.isfile(flname) == False:
        print ('Error: No input file')
        print ('Check that file is in current folder')
        sys.exit()

    data=pd.io(flname,usecols=(0,1))

    hjd,vel,err=[],[],[]
    if len(data[0])==3:
        for i in data:
            hjd.append(i[0]),vel.append(i[1]),err.append(i[2])
    if len(data[0])==2:
        for i in data:
            hjd.append(i[0]),vel.append(i[1])
    hjd,vel=np.array(hjd),np.array(vel)












def res_sin3(pars, x, data=None,sigma=errors_data):
    porb = pars['porb'].value
    hjd0 = pars['hjd0'].value
    gama = pars['gama'].value
    k1 = pars['k1'].value
    #phaseoff=pars['phaseoff'].value
    model=gama+k1*np.sin(2*np.pi*((x-hjd0)/porb))
    if data is None:
        return model
    if sigma is  None:
        return (model - data)
    return (model - data)/sigma
def gaussian(params2, x, data=None,sigma=None):
    amp= params2['amp'].value              # Dilution Factor
    lam0= params2['lam0'].value              # Dilution Factor
    sig= params2['sig'].value              # Dilution Factor
    y0= params2['y0'].value              # Dilution Factor
    model=amp/np.sqrt(2.0*np.pi)/sig*np.exp(-np.power((lam0-x),2)/(2*np.power(sig,2)))+y0
    if data is None:
        return model
    if sigma is  None:
        return (model - data)
    return (model - data)/sigma

def phaser(hjd,hjd0,period):
    phase=(hjd-hjd0)/period-np.fix((hjd-hjd0)/period)
    #print phase
    if np.size(hjd) == 1:
        if phase <0.0:
            return phase+1.0
        else: return phase
    else:
        ss= phase < 0.0
        phase[ss]=phase[ss]+1.0
        return phase
