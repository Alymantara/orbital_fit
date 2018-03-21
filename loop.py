# Performs a circular orbital fit to data
# looping over one of the
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
#import cataclysmic as cv
import lmfit as lm
import os as os
import param as op
import sys
import time as ti
plt.ion() # Activate update
from importlib import reload
reload(op)

def loop(variable,initial,delta,num,):
    import param as op
    reload(op)
    runner=initial+np.arange(num)*delta
    chisqr=[]
    for kk in runner:
        if os.path.isfile(op.file) == False:
            print('Error: No input file')
            print('Check that file is in current folder')
            sys.exit()

        data=np.loadtxt(op.file,usecols=(0,1))

        hjd,vel,err=[],[],[]
        if len(data[0])==3:
            for i in data:
                hjd.append(i[0]),vel.append(i[1]),err.append(i[2])
        if len(data[0])==2:
            for i in data:
                hjd.append(i[0]),vel.append(i[1])
        hjd,vel=np.array(hjd),np.array(vel)

        # %%%%%%%%%%%%%%%%%%%%%%%% ORBIT PARAMETERS

        params = lm.Parameters()
        if variable == 'porb':
            params.add('porb', value= kk, vary=False)
            var_label = 'P$_{orb}$ / d'
        else:
            params.add('porb', value= op.porb, vary=op.fix_porb)
        if variable == 'hjd0':
            params.add('hjd0', value= kk, vary=False)
            var_label = 'HJD$_{0}$ / d'
        else:
            params.add('hjd0',value=op.hjd0, vary=op.fix_hjd0)
        if variable == 'gama':
            params.add('gama', value= kk, vary=False)
            var_label = '$\gamma$ / km s$^{-1}$'
        else:
            params.add('gama',value=op.gama, vary=op.fix_gama)
        if variable == 'k1':
            params.add('k1', value= kk, vary=False)
            var_label = 'K$_{1}$ / km s$^{-1}$'
        else:
            params.add('k1',value=op.k1, vary=op.fix_k1)

        if op.errors:
            err = np.loadtxt(op.file)[:,2]
        else:
            err = op.sigma+np.zeros(len(hjd))

        def res_sin3(pars, x, data=None,sigma=err):
            porb = pars['porb'].value
            hjd0 = pars['hjd0'].value
            gama = pars['gama'].value
            k1 = pars['k1'].value
            model=gama+k1*np.sin(2*np.pi*((x-hjd0)/porb))
            if data is None:
                return model
            if sigma is  None:
                return (model - data)
            return (model - data)/sigma

        results=lm.minimize(res_sin3,params, args=(hjd,vel))
        print('Looping '+str(variable)+': '+str(kk)+'\n')
        lm.report_errors(params,show_correl=False)

        print( 'DoF = ',results.nfree)
        print( 'Chi-squared = ',results.chisqr)
        print( '--------------')
        chisqr.append(results.chisqr)
    chisqr = np.array(chisqr)
    print("Minimum Chi^2: ", runner[chisqr == np.min(chisqr)][0])
    plt.figure(num = 'Loop',facecolor='w')
    plt.clf()
    plt.plot(runner,chisqr/results.nfree,'rs',markersize=6)
    plt.plot(runner,chisqr/results.nfree,'k-')
    plt.axvline(x=runner[chisqr == np.min(chisqr)][0],ls='--',color='b',alpha=0.4)
    plt.xlabel(var_label)
    plt.ylabel('$\chi^2$ / '+str(results.nfree)+' dof')
    plt.show()
    plt.tight_layout()
loop(op.variable,op.initial,op.delta,op.num)
