import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import fnmatch
import os
import re

from mesh import *
from field import *

# ============
# Function that plots time evolution of the disc eccentricity (primarily for binary disc simulations)
# ============
def plotdiscecc():

    # first import global variables
    import par
    
    # first prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.18, right=0.94, top=0.94, bottom=0.12)
    ax = fig.gca()
    
    # if par.physical_units == 'No':
    #     xtitle = r'time [$T_0$]'
    #     ytitle = r'$e_{\rm disc}$'
    # else:
    #     xtitle = 'time [years]'
    #     ytitle = r'$e_{\rm disc}$'        

    xtitle = r'time [$T_0$]'
    ytitle = r'$e_{\rm disc}$'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
    plt.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))


    # several directories can be considered
    directory = par.directory
    if isinstance(par.directory, str) == True:
        directory = [par.directory]

    if par.take_one_point_every == '#':
        take_one_point_every = 1
    else:
        take_one_point_every = par.take_one_point_every

    mytmin = 1e4
    mytmax = 0.0

    # loop over directories
    for j in range(len(directory)):

        # find how many output numbers were produced for each directory
        if par.fargo3d == 'No':
            nboutputs = len(fnmatch.filter(os.listdir(directory[j]), 'gasdens*.dat'))
        else:
            nboutputs = len(fnmatch.filter(os.listdir(directory[j]), 'summary*.dat'))
        print('number of outputs for directory ',directory[j],': ',nboutputs)

        on = range(0,nboutputs-1,take_one_point_every)

        disc_ecc = np.zeros(len(on))
        mytime    = np.zeros(len(on))

        if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
           if len(directory) == 1:
               mylabel = str(par.use_legend)
           else:
               mylabel = str(par.use_legend[j])
        else:
           mylabel = str(directory[j])
        
        first_time = 0

        # loop over output numbers
        for k in range(len(on)):     

            print('output number =',str(k),'out of', str(len(on)),end='\r')

            # get 2D gas surface density field (not compatible with 3D yet...)
            dens = Field(field='dens', fluid='gas', on=on[k], directory=directory[j], physical_units='No', nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile='No', override_units=par.override_units)
            ecc  = Field(field='ecc',  fluid='gas', on=on[k], directory=directory[j], physical_units='No', nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile='No', override_units=par.override_units)

            # things we do only when entering for loop
            if first_time == 0:
            
                first_time = 1
            
                # get surface area of every cell
                surface = np.zeros((dens.nrad,dens.nsec))
                Rinf = dens.redge[0:len(dens.redge)-1]
                Rsup = dens.redge[1:len(dens.redge)]
                surf = np.pi * (Rsup*Rsup - Rinf*Rinf) / dens.nsec
                for th in range(dens.nsec):
                    surface[:,th] = surf

                # get time
                if dens.fargo3d == 'Yes':
                    f1, xpla, ypla, f4, f5, f6, f7, f8, date, omega = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
                else:
                    f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
                with open(directory[j]+"/orbit0.dat") as f_in:
                    firstline_orbitfile = np.genfromtxt(itertools.islice(f_in, 0, 1, None), dtype=float)


            # get total disc mass and time at current output number
            disc_ecc[k] = np.sum(ecc.data*dens.data*surface)/np.sum(dens.data*surface)
            # get time
            mytime[k] = date[take_one_point_every*k]/2.0/np.pi/(1.0**1.5)  # orbital periods at apla=1

            # if par.physical_units == 'Yes':
            #     mytime[k] *= dens.cutime/3.15e7   # in year

        # find minimum anx maximum
        ymin = disc_ecc.min()
        ymax = disc_ecc.max()
        mytmin = mytime.min()
        mytmax = mytime.max()

        # display data as scatter plot for each directory
        ax.scatter(mytime, disc_ecc, s=20, c=par.c20[j], alpha=1.0, label=mylabel)

        
    # set x-range
    if par.mytmin != '#':
        mytmin = par.mytmin
    if par.mytmax != '#':
        mytmax = par.mytmax
    ax.set_xlim(mytmin,mytmax)

    # set y-range
    if ( ('myymin' in open('paramsf2p.dat').read()) and (par.myymin != '#') and (('myymax' in open('paramsf2p.dat').read()) and (par.myymax != '#')) ):
        ymin = par.myymin
        ymax = par.myymax
    ax.set_ylim(ymin,ymax)
    
    # finally add legend
    ax.set_axisbelow(False)
    ax.grid(axis='both', which='major', ls='-', alpha=0.8)

    if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != 'None'):
        legend = plt.legend(loc='upper right',fontsize=15,facecolor='white',edgecolor='white',framealpha=0.85,numpoints=1,bbox_transform=plt.gcf().transFigure)
        for line, text in zip(legend.get_lines(), legend.get_texts()):
            text.set_color(line.get_color())
    
    # And save file
    if len(directory) == 1:           
       outfile = 'discecc_'+str(directory[0])
    else:
       outfile = 'discecc'
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)



# ============
# Function that plots time evolution of the disc's pericenter argument varpi (primarily for binary disc simulations)
# ============
def plotdiscperarg():

    # first import global variables
    import par
    
    # first prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.18, right=0.94, top=0.94, bottom=0.12)
    ax = fig.gca()
    
    # if par.physical_units == 'No':
    #     xtitle = r'time [$T_0$]'
    #     ytitle = r'$\varpi_{\rm disc}$'
    # else:
    #     xtitle = 'time [years]'
    #     ytitle = r'$\varpi_{\rm disc}$'        

    xtitle = r'time [$T_0$]'
    ytitle = r'$\varpi_{\rm disc}$ [rad]'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
    plt.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))


    # several directories can be considered
    directory = par.directory
    if isinstance(par.directory, str) == True:
        directory = [par.directory]

    if par.take_one_point_every == '#':
        take_one_point_every = 1
    else:
        take_one_point_every = par.take_one_point_every

    mytmin = 1e4
    mytmax = 0.0

    # loop over directories
    for j in range(len(directory)):

        # find how many output numbers were produced for each directory
        if par.fargo3d == 'No':
            nboutputs = len(fnmatch.filter(os.listdir(directory[j]), 'gasdens*.dat'))
        else:
            nboutputs = len(fnmatch.filter(os.listdir(directory[j]), 'summary*.dat'))
        print('number of outputs for directory ',directory[j],': ',nboutputs)

        on = range(0,nboutputs-1,take_one_point_every)

        disc_varpi = np.zeros(len(on))
        mytime    = np.zeros(len(on))

        if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
           if len(directory) == 1:
               mylabel = str(par.use_legend)
           else:
               mylabel = str(par.use_legend[j])
        else:
           mylabel = str(directory[j])
        
        first_time = 0

        # loop over output numbers
        for k in range(len(on)):     

            print('output number =',str(k),'out of', str(len(on)),end='\r')

            # get 2D gas surface density field (not compatible with 3D yet...)
            dens = Field(field='dens', fluid='gas', on=on[k], directory=directory[j], physical_units='No', nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile='No', override_units=par.override_units)
            varpi= Field(field='varpi',fluid='gas', on=on[k], directory=directory[j], physical_units='No', nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile='No', override_units=par.override_units)


            # things we do only when entering for loop
            if first_time == 0:
            
                first_time = 1
            
                # get surface area of every cell
                surface = np.zeros((dens.nrad,dens.nsec))
                Rinf = dens.redge[0:len(dens.redge)-1]
                Rsup = dens.redge[1:len(dens.redge)]
                surf = np.pi * (Rsup*Rsup - Rinf*Rinf) / dens.nsec
                for th in range(dens.nsec):
                    surface[:,th] = surf

                # get time
                if dens.fargo3d == 'Yes':
                    f1, xpla, ypla, f4, f5, f6, f7, f8, date, omega = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
                else:
                    f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
                with open(directory[j]+"/orbit0.dat") as f_in:
                    firstline_orbitfile = np.genfromtxt(itertools.islice(f_in, 0, 1, None), dtype=float)


            # get total disc mass and time at current output number
            disc_varpi[k] = np.sum(varpi.data*dens.data*surface)/np.sum(dens.data*surface)
            # get time
            mytime[k] = date[take_one_point_every*k]/2.0/np.pi/(1.0**1.5)  # orbital periods at apla=1.0

            # if par.physical_units == 'Yes':
            #     mytime[k] *= dens.cutime/3.15e7   # in year

        # find minimum anx maximum
        ymin = disc_varpi.min()
        ymax = disc_varpi.max()
        mytmin = mytime.min()
        mytmax = mytime.max()

        # display data as scatter plot for each directory
        ax.scatter(mytime, disc_varpi, s=20, c=par.c20[j], alpha=1.0, label=mylabel)

        
    # set x-range
    if par.mytmin != '#':
        mytmin = par.mytmin
    if par.mytmax != '#':
        mytmax = par.mytmax
    ax.set_xlim(mytmin,mytmax)

    # set y-range
    if ( ('myymin' in open('paramsf2p.dat').read()) and (par.myymin != '#') and (('myymax' in open('paramsf2p.dat').read()) and (par.myymax != '#')) ):
        ymin = par.myymin
        ymax = par.myymax
    ax.set_ylim(ymin,ymax)
    
    # finally add legend
    ax.set_axisbelow(False)
    ax.grid(axis='both', which='major', ls='-', alpha=0.8)

    if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != 'None'):
        legend = plt.legend(loc='upper right',fontsize=15,facecolor='white',edgecolor='white',framealpha=0.85,numpoints=1,bbox_transform=plt.gcf().transFigure)
        for line, text in zip(legend.get_lines(), legend.get_texts()):
            text.set_color(line.get_color())
    
    # And save file
    if len(directory) == 1:           
       outfile = 'discvarpi_'+str(directory[0])
    else:
       outfile = 'discvarpi_'
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)