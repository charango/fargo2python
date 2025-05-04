import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import fnmatch
import re

from mesh import *
from field import *


def plotpowerspectrum():

    # first import global variables
    import par
    
    # Define range of output numbers to consider (in case a time-averaged spectrum is required)
    if par.take_one_point_every == '#':
        take_one_point_every = 1
    else:
        take_one_point_every = par.take_one_point_every

    if np.isscalar(par.on) == False:
        on = range(par.on[0],par.on[1]+1,par.take_one_point_every)
    else:
        on = [par.on]
        #nboutputs = len(fnmatch.filter(os.listdir(par.directory), 'summary*.dat'))
        #on = range(0,nboutputs,take_one_point_every)
    #print('output numbers = ', on)

    # 2D arrays with radius and azimuth
    dens = Field(field='dens', fluid='gas', on=0, directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', override_units=par.override_units)
    pmed2d  = np.zeros((dens.nrad,dens.nsec))
    surface = np.zeros((dens.nrad,dens.nsec))

    for r  in range(dens.nrad):
        pmed2d[r,:] = dens.pmed

    Rinf = dens.redge[0:len(dens.redge)-1]
    Rsup = dens.redge[1:len(dens.redge)]
    surf = np.pi * (Rsup*Rsup - Rinf*Rinf) / dens.nsec
    for th in range(dens.nsec):
        surface[:,th] = surf

    # Range of azimuthal wavenumbers    
    m_min = 1
    m_max = int(dens.nsec/8)
    azi_wavenb = range(m_min,m_max,1)
    # if grid's azimuthal extent is pi: only odd values of m are relevant
    if np.abs(dens.pmed[-1]-dens.pmed[0]-3.14) < 0.1:
        azi_wavenb = range(2*m_min,m_max,2)

    # allocate arrays for Fourier-decomposition
    an = np.zeros(len(azi_wavenb))
    bn = np.zeros(len(azi_wavenb))
    cn = np.zeros(len(azi_wavenb))


    # ========================
    # loop over output numbers
    # ========================
    for k in range(len(on)):

        # get disc midplane density: array of size (nrad, nsec)
        dens = Field(field='dens', fluid='gas', on=on[k], directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', override_units=par.override_units).data

        # total mass
        mass = np.sum(dens*surface)

        # -------------------------------
        # loop over azimuthal wavenumbers
        # -------------------------------
        for m in range(len(azi_wavenb)):

            # real part of Fourier decomposition
            an[m] = (np.sum(dens*surface*np.cos(azi_wavenb[m]*pmed2d))/mass)

            # imaginary part of Fourier decomposition
            bn[m] = (np.sum(dens*surface*np.sin(azi_wavenb[m]*pmed2d))/mass)

            # amplitude
            cn[m] += np.sqrt( an[m]*an[m] + bn[m]*bn[m] )
            if azi_wavenb[m] == 2:
                print('cn[2] = ',  np.sqrt( an[m]*an[m] + bn[m]*bn[m] ) )
    
    # final amplitude - divide by len(on) in case of average over multiple output numbers
    for m in range(len(azi_wavenb)):
        cn[m] /= len(on)
        #print(azi_wavenb[m],cn[m])

    # prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.16, right=0.96, top=0.95, bottom=0.12)
    ax = fig.gca()
    xtitle = 'Azimuthal wavenumber m'
    ytitle = r'Fourier amplitude coefficient $c_m$'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')

    # handle labels
    if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
        mylabel = str(par.use_legend)
    else:
        mylabel = str(par.directory)

    ax.set_yscale('log')
    ax.scatter(azi_wavenb, cn, color=par.c20[0], s=10, label=mylabel)

    # And save file
    outfile = 'power_spectrum_'+str(par.directory)+'_'
    if np.isscalar(par.on) == False:
        outfile += str(par.on[0])+'_'+str(par.on[1])
    else:
        outfile += str(par.on)
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)
