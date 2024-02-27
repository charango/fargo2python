import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import fnmatch
import os
import subprocess
import sys
import re

def plotplanet():

    # first import global variables
    import par

    # first prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.17, right=0.94, top=0.94, bottom=0.12)
    ax = fig.gca()
    if par.plot_planet[0] == 't':
        xtitle = 'Time [orbits]'
    if par.plot_planet[0] == 'a':
        xtitle = 'Semi-major axis'
        if par.physical_units == 'Yes':
            xtitle += ' [au]'
    if par.plot_planet[1] == 'a':
        ytitle = 'Semi-major axis'
        if par.physical_units == 'Yes':
            ytitle += ' [au]'
    if par.plot_planet[1] == 'r':
        ytitle = 'Orbital radius'
        if par.physical_units == 'Yes':
            ytitle += ' [au]'
    if par.plot_planet[1] == 'e':
        ytitle = 'Eccentricity'
    if par.plot_planet[1] == 'i':
        ytitle = 'Inclination [deg]'
    if par.plot_planet[1] == 'm':
        ytitle = 'Planet mass'
        if par.physical_units == 'Yes':
            ytitle += ' [Earth mass]'
    if par.plot_planet[1] == 'p':
        ytitle = 'Orbital period ratio'
    if par.plot_planet[1] == 'mmr':
        plt.subplots_adjust(left=0.12, right=0.90, top=0.95, bottom=0.12)
        ax2 = ax.twinx()
        if par.mmr_integers[0] == 1:
            str0 = ''
        else:
            str0 = str(par.mmr_integers[0])
        if par.mmr_integers[1] == 1:
            str1 = ''
        else:
            str1 = str(par.mmr_integers[1])
        if par.mmr_integers[0]-par.mmr_integers[1] == 1:
            strdiff = ''
        else:
            strdiff = str(par.mmr_integers[0]-par.mmr_integers[1])
        ytitle = r'$\Psi_{\rm i}$ = '+str0+r'$\lambda_{\rm o} - $'+str1+r'$\lambda_{\rm i} - $'+strdiff+r'$\varpi_{\rm i}$'
        ytitle_o = r'$\Psi_{\rm o}$ = '+str0+r'$\lambda_{\rm o} - $'+str1+r'$\lambda_{\rm i} - $'+strdiff+r'$\varpi_{\rm o}$'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=False, useMathText=True)

    # several directories are possible
    directory = par.directory
    if isinstance(par.directory, str) == True:
        directory = [par.directory]
    
    # loop over directories
    for j in range(len(directory)):

        # work out physical units
        if par.physical_units == 'Yes':
            if par.fargo3d == 'No':
            # units.dat contains physical units of mass [kg], length [m], time [s], and temperature [k] 
                cumass, culength, cutime, cutemp = np.loadtxt(directory[j]+"/units.dat",unpack=True)
            else:
            # get units via variable.par file
                command = 'awk " /^UNITOFLENGTHAU/ " '+directory[j]+'/variables.par'
                # check which version of python we're using
                if sys.version_info[0] < 3:   # python 2.X
                    buf = subprocess.check_output(command, shell=True)
                else:                         # python 3.X
                    buf = subprocess.getoutput(command)
                culength = float(buf.split()[1])*1.5e11  #from au to meters
                command = 'awk " /^UNITOFMASSMSUN/ " '+directory[j]+'/variables.par'
                # check which version of python we're using
                if sys.version_info[0] < 3:   # python 2.X
                    buf = subprocess.check_output(command, shell=True)
                else:                         # python 3.X
                    buf = subprocess.getoutput(command)
                cumass = float(buf.split()[1])*2e30  #from Msol to kg
                # unit of time = sqrt( pow(L,3.) / 6.673e-11 / M );
                cutime = np.sqrt( culength**3.0 / 6.673e-11 / cumass)
        
        if par.fargo3d == 'No':
            time, e, a, M, V, PA, mylambda, varpi = np.loadtxt(directory[j]+"/orbit0.dat",unpack=True)
        else:
            command = 'awk "{print NF; exit}" '+directory[j]+'/orbit0.dat'
            # check which version of python we're using
            if sys.version_info[0] < 3:   # python 2.X
                buf = subprocess.check_output(command, shell=True)
            else:                         # python 3.X
                buf = subprocess.getoutput(command)
            nbcol = float(buf.split()[0])
            if nbcol == 11:
                time, e, a, M, V, argPA, phiangle, incl, longAN, PA, mylambda = np.loadtxt(directory[j]+"/orbit0.dat",unpack=True)
            else:
                time, e, a, M, V, argPA, phiangle, incl, longAN, PA = np.loadtxt(directory[j]+"/orbit0.dat",unpack=True)
           
        rpla_0 = a[0]
       
        # count how many planets 
        nbplanets = len(fnmatch.filter(os.listdir(directory[j]), 'orbit*.dat'))
        if (nbplanets <= 1 and par.plot_planet[1] == 'p'):
            sys.exit('ERROR: you requested to plot orbital period ratio, but there is only one planet simulated in directory: ', directory[j])
           
        # now, read orbitN.dat or bigplanetN.dat files
        for k in range(nbplanets):
            
            if par.fargo3d == 'No':
                time, e, a, M, V, PA, mylambda, varpi = np.loadtxt(directory[j]+"/orbit"+str(k)+".dat",unpack=True)
            else:
                if nbcol == 11:
                    time, e, a, M, V, argPA, phiangle, incl, longAN, PA, mylambda = np.loadtxt(directory[j]+"/orbit"+str(k)+".dat",unpack=True)
                else:
                    time, e, a, M, V, argPA, phiangle, incl, longAN, PA = np.loadtxt(directory[j]+"/orbit"+str(k)+".dat",unpack=True)
            time /= (2.0*np.pi*rpla_0*np.sqrt(rpla_0))  # time in orbital periods at inner planet's initial location
            
            if (par.mytmin == '#'):
                xmin = time.min()
                umin = 0
            else:
                xmin = par.mytmin
                umin = np.argmin(np.abs(time-xmin))
            if (par.mytmax == '#'):
                xmax = time.max()
                umax = len(time)-1
            else:
                xmax = par.mytmax
                umax = np.argmin(np.abs(time-xmax))
                
            if par.plot_planet[0] == 't':                              
                #x = time[umin:umax]
                x = time[umin:umax+1]
                ax.set_xlim(xmin,xmax)
                
            if par.plot_planet[0] == 'a':
                x = a[umin:umax+1]
                
            if par.plot_planet[1] == 'e':
                y = e[umin:umax+1]
                
            if par.plot_planet[1] == 'i':
                y = incl[umin:umax+1]*180.0/np.pi
                
            if par.plot_planet[1] == 'a':
                y = a[umin:umax+1]
                if par.physical_units == 'Yes':
                    y *= (culength / 1.5e11) # in au
                    
            if par.plot_planet[1] == 'r' or par.plot_planet[1] == 'm':
                if par.fargo3d == 'No':
                    f1, xpla, ypla, f4, f5, mpla, f7, date, f9, f10, f11 = np.loadtxt(directory[j]+"/bigplanet"+str(k)+".dat",unpack=True)
                else:
                    f1, xpla, ypla, zpla, f5, f6, f7, mpla, date, f10 = np.loadtxt(directory[j]+"/bigplanet"+str(k)+".dat",unpack=True)
                mytime = date/(2.0*np.pi*rpla_0*np.sqrt(rpla_0))
                # repeat stuff above just in case the bigplanet and orbit.dat files have different lengths...
                if (par.mytmin == '#'):
                    umin = 0
                else:
                    umin = np.argmin(np.abs(mytime-xmin))
                if (par.mytmax == '#'):
                    umax = len(mytime)-1
                else:
                    umax = np.argmin(np.abs(mytime-xmax))
                x = mytime[umin:umax+1]
                if par.plot_planet[1] == 'r':
                    if par.fargo3d == 'No':
                        y = np.sqrt( xpla*xpla + ypla*ypla )
                    else:
                        y = np.sqrt( xpla*xpla + ypla*ypla + zpla*zpla )
                    y = y[umin:umax+1]
                    if par.physical_units == 'Yes':
                        y *= (culength / 1.5e11) # in au
                if par.plot_planet[1] == 'm':
                    y = mpla[umin:umax+1]
                    if par.physical_units == 'Yes':
                        y *= (cumass / 2e30 / 3e-6) # in Earth masses
                        
            if par.plot_planet[1] == 'p':
                if k == 0:
                    p0 = a**(1.5)
                if k == 1:
                    p1 = a**(1.5)
                    
            if par.plot_planet[1] == 'mmr':
                if k == 0:
                    lambda0 = mylambda
                    varpi0  = PA
                if k == 1:
                    lambda1 = mylambda
                    varpi1  = PA
                    
            if par.plot_planet[1] != 'p' and par.plot_planet[1] != 'mmr':
                if len(directory) == 1:
                    ax.legend_ = None
                else:
                    ax.legend(frameon=False,fontsize=15)
                ax.plot(x[::par.take_one_point_every], y[::par.take_one_point_every], color=par.c20[k*len(directory)+j], lw=2., linestyle = 'solid', label=directory[j])
                fig.add_subplot(ax)
                
        if par.plot_planet[1] == 'p':
            y = p1/p0
            if len(directory) == 1:
                ax.legend_ = None
            ax.plot(x[::par.take_one_point_every], y[::par.take_one_point_every], color=par.c20[j], lw=2., linestyle = 'solid', label=directory[j])
            fig.add_subplot(ax)
            
        if par.plot_planet[1] == 'mmr':
            y0 = par.mmr_integers[0]*lambda1 - par.mmr_integers[1]*lambda0 - (par.mmr_integers[0]-par.mmr_integers[1])*varpi0
            y1 = par.mmr_integers[0]*lambda1 - par.mmr_integers[1]*lambda0 - (par.mmr_integers[0]-par.mmr_integers[1])*varpi1
            for l1 in range(len(y0)):
                while y0[l1] >= 2.0*np.pi:
                    y0[l1] -= 2.0*np.pi
                while y0[l1] < 0:
                    y0[l1] += 2.0*np.pi
            for l2 in range(len(y1)):
                while y1[l2] >= 2.0*np.pi:
                    y1[l2] -= 2.0*np.pi
                while y1[l2] < 0:
                    y1[l2] += 2.0*np.pi
            ax.yaxis.label.set_color(color=par.c20[2*j])
            ax.scatter(x[::par.take_one_point_every], y0[::par.take_one_point_every], s=5, c=par.c20[2*j],   alpha=1.0)
            ax2.scatter(x[::par.take_one_point_every], y1[::par.take_one_point_every], s=5, c=par.c20[2*j+1], alpha=1.0)
            ax.legend(frameon=False,fontsize=15)
            ax2.yaxis.label.set_color(color=par.c20[2*j+1])
            ax2.set_ylabel(ytitle_o)
            fig.add_subplot(ax)
            
    # save file
    if len(directory) == 1:           
        outfile = par.plot_planet[0]+'_'+par.plot_planet[1]+'_'+str(directory[0])
    else:
        outfile = par.plot_planet[0]+'_'+par.plot_planet[1]
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)
