# Orbital fitting program
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import matplotlib.cm as cm
import lmfit as lm
import os as os
import param as op
import sys
import time as ti
plt.ion() # Activate update
from importlib import reload
reload(op)
import corner

if os.path.isfile(op.file) == False:
    print ('Error: No input file')
    print ('Check that file is in current folder')
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
params.add('porb', value= op.porb, vary=op.fix_porb)
params.add('hjd0',value=op.hjd0-int(op.hjd0), vary=op.fix_hjd0,min=op.hjd0-int(op.hjd0)-op.hjd0/2.,max=op.hjd0-int(op.hjd0)+op.hjd0/2.)
params.add('gama', value= op.gama, vary=op.fix_gama)
params.add('k1', value= op.k1, vary=op.fix_k1)
#params.add('phaseoff', value= op.phaseoff, vary=op.fix_phaseoff)

if op.sigma == None:
    errors_data = np.loadtxt(op.file)
    errors_data = errors_data[:,2]
else:
    errors_data = op.sigma*np.ones(len(hjd))

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

results=lm.minimize(res_sin3,params, args=(hjd-int(op.hjd0),vel))
print ('Best fit for '+op.object)
print ('--------------')
lm.report_errors(results.params,show_correl=False)
print ('--------------')
print ('DoF = ',results.nfree)
print ('Chi-squared = ',results.chisqr)

if op.scale_errors:
    rescale = np.sqrt(results.chisqr/results.nfree)
else:
    rescale = 1.0

print ('Rescale = ',rescale)
phase1=[]
phase = phaser(hjd,results.params['hjd0'].value+int(op.hjd0),results.params['porb'].value)

fig=plt.figure(num='Orbital',figsize=(12,8),facecolor='w')
plt.clf()
ax1=plt.subplot2grid((4, 1), (0, 0),rowspan=3)



plt.errorbar(phase,vel,yerr=errors_data*rescale,marker='o',ms=3,mec='k',ecolor='gray',mfc='k',linestyle='None',capsize=0)
plt.errorbar(phase+1.0,vel,yerr=errors_data*rescale,marker='o',ms=3,mec='k',ecolor='gray',mfc='k',linestyle='None',capsize=0)
plt.errorbar(phase-1.0 ,vel,yerr=errors_data*rescale,marker='o',ms=3,mec='k',ecolor='gray',mfc='k',linestyle='None',capsize=0)

#plt.plot(np.arange(0.,1.,0.01),ffpar[0][0]+ffpar[0][1]*np.sin((np.arange(0.,1.,0.01)-0.05)*2*np.pi))
xx = np.arange(-1,3,0.01)
plt.plot(xx,results.params['gama'].value+results.params['k1'].value*np.sin(2.0*np.pi*xx),'b-')
plt.axis([0,2,min(vel)*op.plotlim,max(vel)*op.plotlim])
plt.ylabel('Radial Velocity, km s$^{-1}$')

ax1.yaxis.set_major_locator(plt.MaxNLocator(9))
ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
plt.minorticks_on()
plt.setp(ax1.get_xticklabels(), visible=False)
import datetime
plt.title('Object: '+op.object+'\n Filename: '+op.file+'  '+datetime.datetime.fromtimestamp(ti.time()).strftime('%Y-%m-%d %H:%M:%S'))

size=18
plt.xticks(size=size)
plt.yticks(size=size)


ax2=plt.subplot2grid((4, 1), (3, 0))
res=vel-res_sin3(results.params,hjd-np.int(op.hjd0))

plt.errorbar(phase,res,yerr=errors_data*rescale,marker='o',ms=3,mec='k',ecolor='gray',mfc='k',linestyle='None',capsize=0)
plt.errorbar(phase+1.0,res,yerr=errors_data*rescale,marker='o',ms=3,mec='k',ecolor='gray',mfc='k',linestyle='None',capsize=0)
plt.errorbar(phase-1.0,res,yerr=errors_data*rescale,marker='o',ms=3,mec='k',ecolor='gray',mfc='k',linestyle='None',capsize=0)
ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
ax2.xaxis.set_major_locator(plt.MaxNLocator(4))
plt.minorticks_on()
plt.xlabel('Orbital Phase, $\phi$')
plt.axhline(y=0.0,linestyle='--',color='k')
plt.ylabel('Residuals')
plt.axis([0,2,min(res)*op.lim_res,max(res)*op.lim_res])
plt.tight_layout()
plt.subplots_adjust(hspace = .000)

plt.savefig(op.psname+'.'+op.output,dpi=300,facecolor='w',format=op.output)

'''
def res_sin3(pars, x, data=None,sigma=op.sigma+np.zeros(len(phase))):
    porb = pars['porb'].value
    gama = pars['gama'].value
    k1 = pars['k1'].value
    phaseoff=pars['phaseoff'].value
    model=gama+k1*np.sin(2.*np.pi*(x+phaseoff))
    if data is None:
        return model
    if sigma is  None:
        return (model - data)
    return (model - data)/sigma
'''

print( '---------')
if op.do_bootstrap:
    grid_boot = np.zeros(0,dtype=[('k1', 'f8'),('phase', 'f8'), ('hjd0', 'f8'),('gama','f8'), ('porb','f8')])
    num1 = len(phase)
    for ii in np.arange(op.boot_iter):
        tt1 = np.sort(np.array((np.random.rand(1,num1)[0]*num1)).astype('int'))
        params3 = lm.Parameters()
        params3.add('porb', value= op.porb, vary=op.fix_porb)
        params3.add('gama', value= op.gama, vary=op.fix_gama)
        params3.add('k1', value= op.k1, vary=op.fix_k1)
        params3.add('hjd0',value=op.hjd0 - int(op.hjd0), vary=op.fix_hjd0,min=op.hjd0-int(op.hjd0)-op.hjd0/2.,max=op.hjd0-int(op.hjd0)+op.hjd0/2.)
        #params3.add('phaseoff', value= op.phaseoff, vary=op.fix_phaseoff)
        res_boot=lm.minimize(res_sin3,params3, args=(hjd[tt1]-int(op.hjd0),vel[tt1]))


        if op.fix_porb:
            phaoff_boot = (res_boot.params['hjd0'].value+int(op.hjd0) - op.hjd0)/res_boot.params['porb'].value
        else:
            phaoff_boot = (res_boot.params['hjd0'].value+int(op.hjd0) - op.hjd0)/op.porb

        grid_boot = np.append(grid_boot, np.array([(res_boot.params['k1'].value,
                     phaoff_boot,res_boot.params['hjd0'].value+int(op.hjd0),res_boot.params['porb'].value,
                     res_boot.params['gama'].value)], dtype=grid_boot.dtype))
    #######################
    ##### CORRELATION PLOTS
    #######################
    data = np.array(grid_boot.tolist())
    if op.fix_porb:
        #print("yes hjd0, yes per")
        labels = [r"$K_1$", r"$\phi_0$",r"HJD$_0$",r"$\gamma$"]
        data = np.delete(data, 4, axis=1)
        data[:,2] -= int(op.hjd0)
    else:
        labels = [r"$K_1$", r"$\phi_0$",r"HJD$_0$",r"$\gamma$",r"P$_{orb}$"]
        #data[:,3] -= hd_fac
    corner.corner(data,labels=labels,show_titles=True,quantiles=[0.15,0.5, .86],top_ticks=False)
    #plt.savefig(output_root+'_bootstrap.'+output_ext,dpi=300)
print ('Errors Least Squares:  K1 = ',results.params['k1'].value,' +/-',results.params['k1'].stderr)
if op.do_bootstrap: rint ('Errors Bootstrap:  K1 = ',np.mean(grid_boot['k1']),' +/-',np.std(grid_boot['k1']))
print ('------------')
print ('Errors Least Squares:  Gamma = ',results.params['gama'].value,' +/-',results.params['gama'].stderr)
if op.do_bootstrap: print ('Errors Bootstrap:  Gamma = ',np.mean(grid_boot['gama']),' +/-',np.std(grid_boot['gama']))
print ('------------')
print ('Errors Least Squares:  Porb = ',results.params['porb'].value,' +/-',results.params['porb'].stderr)
if op.do_bootstrap:print ('Errors Bootstrap:  Porb = ',np.mean(grid_boot['porb']),' +/-',np.std(grid_boot['porb']))
print ('------------')
print ('Errors Least Squares:  HJD0 = ',results.params['hjd0'].value+int(op.hjd0),' +/-',results.params['hjd0'].stderr)
if op.do_bootstrap:print ('Errors Bootstrap:  HJD0 = ',np.mean(grid_boot['hjd0'])+int(op.hjd0),' +/-',np.std(grid_boot['hjd0']))
print ('------------')
lori = np.sqrt( (1./results.params['porb'].value)**2*(results.params['hjd0'].stderr)**2 + ( (results.params['hjd0'].value -op.hjd0)**2/(results.params['porb'].value)**2 )* results.params['porb'].stderr**2 )
print ('Errors Least Squares:  PhaseOff = ',((results.params['hjd0'].value - op.hjd0)/results.params['porb'].value)%1,' +/-',lori)
if op.do_bootstrap: lori = np.sqrt( (1./results.params['porb'].value)**2*(np.std(grid_boot['hjd0']))**2 + ( (results.params['hjd0'].value -op.hjd0)**2/(results.params['porb'].value)**2 )* np.std(grid_boot['porb'])**2 )
if op.do_bootstrap:print ('Errors Bootstrap:  PhaseOff = ',((results.params['hjd0'].value - op.hjd0)/results.params['porb'].value)%1,' +/-',lori)
print ('------------')
