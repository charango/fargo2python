import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
import fnmatch
import os
import re

from field import *

def plottqwk():

    # first import global variables
    import par

    # first prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.17, right=0.94, top=0.94, bottom=0.12)
    ax = fig.gca()
    ax.set_xlabel('time [orbits]')
    if par.plot_tqwk == 'torque':
        ytitle = 'Specific torque on planet'
    if par.plot_tqwk == 'normtorque':
        ytitle = r'$\Gamma / \Gamma_0$'
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

        if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
            if len(directory) == 1:
                mylabel = str(par.use_legend)
            else:
                mylabel = str(par.use_legend[j])
        else:
            mylabel = str(directory[j])

        # Locally check if simulations were carried out with Fargo3D
        summary0_file = directory[j]+'/summary0.dat'
        if os.path.isfile(summary0_file) == True:
            fargo3d = 'Yes'
        else:
            fargo3d = 'No'

        # start by reading planet0.dat file to get the initial radial position of the planet
        if fargo3d == 'Yes':
            f1, xpla, ypla, f4, f5, f6, f7, mpla, date, omega = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
        else:
            f1, xpla, ypla, f4, f5, mpla, f7, date, omega, f10, f11 = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
        rpla_0 = np.sqrt( xpla[0]*xpla[0] + ypla[0]*ypla[0] )
        rpla_0 = 1.0  # CUIDADIN!
        # count how many planets 
        nbplanets = len(fnmatch.filter(os.listdir(directory[j]), 'orbit*.dat'))  
        print('nbplanets = ',nbplanets)

        # Normalized torque by Gamma_0 = q/h^2 x Sigma(r_p) r_p^4 Omega^2(r_p)
        if par.plot_tqwk == 'normtorque':
            # get planet-to-star mass ratio q
            q = mpla[len(mpla)-1]   # time-varying array
            # get local disc's aspect ratio
            if fargo3d == 'Yes':
                command  = par.awk_command+' " /^ASPECTRATIO/ " '+directory[j]+'/*.par'
                command2 = par.awk_command+' " /^FLARINGINDEX/ " '+directory[j]+'/*.par'
            else:
                command  = par.awk_command+' " /^AspectRatio/ " '+directory[j]+'/*.par'
                command2 = par.awk_command+' " /^FlaringIndex/ " '+directory[j]+'/*.par'
            buf = subprocess.getoutput(command)
            aspectratio = float(buf.split()[1])
            buf2 = subprocess.getoutput(command2)
            fli = float(buf2.split()[1])
            rpla0_normtq = np.sqrt( xpla[0]*xpla[0] + ypla[0]*ypla[0] )
            h = aspectratio*(rpla0_normtq**fli)  # constant in time
            # get local azimuthally averaged surface density
            myfield0 = Field(field='dens', fluid='gas', on=0, directory=directory[j], physical_units='No', nodiff='Yes', fieldofview=par.fieldofview, slice=par.slice, onedprofile='Yes', override_units=par.override_units)
            dens = np.sum(myfield0.data,axis=1) / myfield0.nsec
            imin = np.argmin(np.abs(myfield0.rmed-rpla0_normtq))
            sigmap = dens[imin]
            # Finally infer Gamma_0
            Gamma_0 = (q/h/h)*sigmap*rpla0_normtq
            print('q = ', q)
            print('h = ', h)
            print('rpla0_normtq = ', rpla0_normtq)
            print('sigmap = ', sigmap)
            print('Gamma_0 = ', Gamma_0)


        # now, read tqwk0.dat file
        for k in range(nbplanets): 
            f1, it, ot, f4, f5, ip, op, f8, f9, time = np.loadtxt(directory[j]+"/tqwk"+str(k)+".dat",unpack=True)
            time /= (2.0*np.pi*rpla_0*np.sqrt(rpla_0))  # time in orbital periods at inner planet's initial location
            tq = it+ot
            pw = ip+op
            if par.plot_tqwk == 'torque' or par.plot_tqwk == 'rtatorque' or par.plot_tqwk == 'normtorque':
                y = tq
            if par.plot_tqwk == 'power' or par.plot_tqwk == 'rtapower':
                y = pw
            if par.plot_tqwk == 'rtatorque' or par.plot_tqwk == 'rtapower':
                for i in range(1,len(y)):
                    y[i] = (i*y[i-1] + tq[i])/(i+1.0)
            if (par.mytmin == '#'):
                xmin = time.min()
                imin = 0.0
            else:
                xmin = par.mytmin
                imin = np.argmin(np.abs(time-xmin))
            if (par.mytmax == '#'):
                xmax = time.max()
                imax = len(time)-1
            else:
                xmax = par.mytmax
                imax = np.argmin(np.abs(time-xmax))
            ax.set_xlim(xmin,xmax)

            ymin = 0.0
            ymax = 0.0
            if ('myymin' in open('paramsf2p.dat').read()) and (par.myymin != '#'):
                ymin = par.myymin
            #else:
            #    ymin = y[imin:imax].min()
            if ('myymax' in open('paramsf2p.dat').read()) and (par.myymax != '#'):
                ymax = par.myymax
            #else:
            #    ymax = y[imin:imax].max()
            if ymin != 0.0 or ymax != 0.0:
                ax.set_ylim(ymin,ymax)
               
            # new (Nov. 2023): display in y-axis log scale (indirect term
            # project)
            if par.log_xyplots_y == 'Yes':
                y = np.abs(y)
                ax.set_yscale('log')
                ytitle = str('|')+ytitle+str('|')

            if par.plot_tqwk == 'normtorque':
                y /= Gamma_0
               
            ax.plot(time[::par.take_one_point_every], y[::par.take_one_point_every], color=par.c20[k*len(directory)+j], lw=2., linestyle = 'solid', label=mylabel)
            ax.legend(frameon=False,fontsize=15)
            fig.add_subplot(ax)

    ax.set_axisbelow(False)
    ax.grid(axis='both', which='major', ls='-', alpha=0.8)
    legend = plt.legend(loc='upper right',fontsize=15,facecolor='white',edgecolor='white',framealpha=0.85,numpoints=1,bbox_transform=plt.gcf().transFigure)
    for line, text in zip(legend.get_lines(), legend.get_texts()):
        text.set_color(line.get_color())
              
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
