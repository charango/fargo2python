import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import fnmatch
import re

from mesh import *
from field import *

def plotlibcross():

    # first import global variables
    import par
    
    # first prepare figures
    '''
    font = {'family' : 'Arial',
            'size'   : 16}
    matplotlib.rc('font', **font)
    '''

    fig, axs = plt.subplots(3,1,figsize=(8, 16))

    plt.subplots_adjust(left=0.16, right=0.95, top=0.97, bottom=0.06)
    #ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
    plt.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))

    axs[0].set_ylabel(r'$\omega_{\rm lib} / \omega_{\rm cross}$')
    axs[1].set_ylabel(r'$\omega^{-1}_{\rm lib}$')
    axs[2].set_ylabel(r'$\omega^{-1}_{\rm cross}$')

    for ax in axs.flat:
        ax.set(xlabel='time [orbits]')

    # several directories can be considered
    directory = par.directory
    if isinstance(par.directory, str) == True:
        directory = [par.directory]

    if par.take_one_point_every == '#':
        take_one_point_every = 1
    else:
        take_one_point_every = par.take_one_point_every


    # loop over directories
    for j in range(len(directory)):

        # get density at on=0
        dens = Field(field='dens', fluid='gas', on=0, directory=directory[j], physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile='No', override_units=par.override_units)

        # find how many output numbers were produced for each directory
        if dens.fargo3d == 'No':
            nboutputs = len(fnmatch.filter(os.listdir(directory[j]), 'gasdens*.dat'))
        else:
            nboutputs = len(fnmatch.filter(os.listdir(directory[j]), 'summary*.dat'))
        print('number of outputs for directory ',directory[j],': ',nboutputs)

        on = range(0,nboutputs,take_one_point_every)

        omega_lib   = np.zeros(len(on))
        omega_cross = np.zeros(len(on))
        ratio       = np.zeros(len(on))
        mytime      = np.zeros(len(on))

        # get time
        if dens.fargo3d == 'Yes':
            f1, xpla, ypla, f4, f5, f6, f7, f8, date, omega = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
        else:
            f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
        with open(directory[j]+"/orbit0.dat") as f_in:
            firstline_orbitfile = np.genfromtxt(itertools.islice(f_in, 0, 1, None), dtype=float)
        apla = firstline_orbitfile[2]
        
        # handle labels
        if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
           if len(directory) == 1:
               mylabel = str(par.use_legend)
           else:
               mylabel = str(par.use_legend[j])
        else:
           mylabel = str(directory[j])
        

        # azimuthal index of sector where planet is
        jpla = dens.nsec//2  # cuidadin, only true if frame is corotating with planet!

        # azimuthal index of sector about where librating gas is (about 1 radian ahead of planet)
        jlib = jpla + int(1.0*dens.nsec/(2.0*np.pi))
        jlib_range = int(np.abs(0.1*(jlib-jpla)))


        # loop over output numbers
        for k in range(len(on)):     

            print('output number =',str(k),'out of', str(len(on)),end='\r')

            # get 2D gas surface density field (not compatible with 3D yet...)  (nrad,nsec)
            vortensity = Field(field='vortensity', fluid='gas', on=on[k], directory=directory[j], physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile='No', override_units=par.override_units).data

            # time
            mytime[k] = date[take_one_point_every*k]/2.0/np.pi/(apla**1.5)  # orbital periods at apla

            # find planet position and radial index of ring where planet is
            rpla = np.sqrt( xpla[take_one_point_every*k]*xpla[take_one_point_every*k] + ypla[take_one_point_every*k]*ypla[take_one_point_every*k] )
            ipla = np.argmin(np.abs(dens.rmed-rpla))

            # average inverse vortensity of librating gas, estimated at planet's orbital radius 
            # and averaged over a small range of azimuthal angles around theta_pla + 1 radian (by eye inspection!)
            omega_lib[k] = np.mean(vortensity[ipla-1:ipla+1,jlib-jlib_range:jlib+jlib_range])
            omega_lib[k] = 1./omega_lib[k]

            # average inverse vortensity of orbit-crossing flow, estimated at planet's orbital 
            # radius and in a range of azimuths between that of the planet and 0.2 radians below 
            # the planet (by-eye inspection!)
            jcross = jpla - int(0.2*dens.nsec/(2.0*np.pi))
            omega_cross[k] = np.mean(vortensity[ipla-1:ipla+1,jcross:jpla])
            omega_cross[k] = 1./omega_cross[k]

            # ratio of librating and orbit-crossing inverse vortensities
            ratio[k] = omega_lib[k] / omega_cross[k]

            print(on[k], mytime[k], rpla, ipla, 1./omega_lib[k], 1./omega_cross[k], ratio[k])


        # display data as scatter plot for each directory
        axs[0].scatter(mytime, ratio, s=20, c=par.c20[j], alpha=1.0, label=mylabel)
        axs[1].scatter(mytime, 1./omega_lib, s=20, c=par.c20[j], alpha=1.0, label=mylabel)
        axs[2].scatter(mytime, 1./omega_cross, s=20, c=par.c20[j], alpha=1.0, label=mylabel)


    # finally add legend
    '''
    legend = plt.legend(loc='upper right',fontsize=15,facecolor='white',edgecolor='white',framealpha=0.85,numpoints=1,bbox_transform=plt.gcf().transFigure)
    for line, text in zip(legend.get_lines(), legend.get_texts()):
        text.set_color(line.get_color())
    '''
    for ax in axs.flat:
        ax.set(xlim=(0,mytime.max()))

    # And save file
    if len(directory) == 1:           
       outfile = 'wlib_over_wcross_'+str(directory[0])
    else:
       outfile = 'wlib_over_wcross'
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)
