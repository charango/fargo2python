import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import fnmatch
import os
import re

def plottqwk():

    # first import global variables
    import par
    
    # first prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.19, right=0.94, top=0.94, bottom=0.12)
    ax = fig.gca()
    ax.set_xlabel('time [orbit]')
    if par.plot_tqwk == 'torque':
        ytitle = 'Specific torque on planet'
    if par.plot_tqwk == 'rtatorque':
        ytitle = 'Specific r.t.a. torque on planet'
    if par.plot_tqwk == 'power':
        ytitle = 'Specific power on planet'
    if par.plot_tqwk == 'rtapower':
        ytitle = 'Specific r.t.a. power on planet'
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
    plt.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    
    # several directories are possible
    directory = par.directory
    if isinstance(par.directory, str) == True:
        directory = [par.directory]

    # loop over directories
    for j in range(len(directory)): 
       # start by reading planet0.dat file to get the initial radial position of the planet
       if par.fargo3d == 'Yes':
           f1, xpla, ypla, f4, f5, f6, f7, f8, date, omega = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
       else:
           f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
       rpla_0 = np.sqrt( xpla[0]*xpla[0] + ypla[0]*ypla[0] )
       # count how many planets 
       nbplanets = len(fnmatch.filter(os.listdir(directory[j]), 'orbit*.dat'))    
       # now, read tqwk0.dat file
       for k in range(nbplanets): 
           f1, it, ot, f4, f5, ip, op, f8, f9, time = np.loadtxt(directory[j]+"/tqwk"+str(k)+".dat",unpack=True)
           time /= (2.0*np.pi*rpla_0*np.sqrt(rpla_0))  # time in orbital periods at inner planet's initial location
           tq = it+ot
           pw = ip+op
           if par.plot_tqwk == 'torque' or par.plot_tqwk == 'rtatorque':
               y = tq
           if par.plot_tqwk == 'power' or par.plot_tqwk == 'rtapower':
               y = pw
           if par.plot_tqwk == 'rtatorque' or par.plot_tqwk == 'rtapower':
               for i in range(1,len(y)):
                   y[i] = (i*y[i-1] + tq[i])/(i+1.0)
           if (par.mytmin == '#'):
               xmin = time.min()
           else:
               xmin = par.mytmin
           if (par.mytmax == '#'):
               xmax = time.max()
           else:
               xmax = par.mytmax
           ax.set_xlim(xmin,xmax)

           ymin = 0.0
           ymax = 0.0
           if ('myymin' in open('paramsf2p.dat').read()) and (par.myymin != '#'):
               ymin = par.myymin
           if ('myymax' in open('paramsf2p.dat').read()) and (par.myymax != '#'):
               ymax = par.myymax
           if ymin != 0.0 or ymax != 0.0:
               ax.set_ylim(ymin,ymax)
               
           
           ax.plot(time[::par.take_one_point_every], y[::par.take_one_point_every], color=par.c20[k*len(directory)+j], lw=2., linestyle = 'solid', label=directory[j])
           ax.legend(frameon=False,fontsize=15)
           fig.add_subplot(ax)
              
    # save file
    if len(directory) == 1:           
       outfile = par.plot_tqwk+'_'+str(directory[0])
    else:
       outfile = par.plot_tqwk
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)
