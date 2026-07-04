import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import fnmatch
import re
import os

from mesh import *
from field import *


# Function that plots time evolution of mode with azimuthal wavenumber 
# m set in paramsf2p.dat parameter file. Based on density field.

def plotfourier():

    # first import global variables
    import par
    
    # Azimuthal wavenumber
    # if isinstance(par.plot_fourier, int) == 'True':
    #     azi_wavenb = par.plot_fourier
    # else:
    #     azi_wavenb = 1
    #     print('You have entered a non-integer value for plot_fourier in paramsf2p.dat. Defaulting to m=1!')
    azi_wavenb = int(par.plot_fourier)

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
    dens = Field(field='dens', fluid='gas', on=0, directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units)
    pmed2d  = np.zeros((dens.nrad,dens.nsec))
    surface = np.zeros((dens.nrad,dens.nsec))

    for r  in range(dens.nrad):
        pmed2d[r,:] = dens.pmed

    Rinf = dens.redge[0:len(dens.redge)-1]
    Rsup = dens.redge[1:len(dens.redge)]
    surf = np.pi * (Rsup*Rsup - Rinf*Rinf) / dens.nsec
    for th in range(dens.nsec):
        surface[:,th] = surf

    # get time
    if dens.fargo3d == 'Yes':
        f1, xpla, ypla, f4, f5, f6, f7, f8, date, omega = np.loadtxt(par.directory+"/planet0.dat",unpack=True)
    else:
        f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(par.directory+"/planet0.dat",unpack=True)

    mytime = np.zeros(len(on))
    cn = np.zeros(len(on))

    # ========================
    # loop over output numbers
    # ========================
    for k in range(len(on)):

        print('output number =',str(k),'out of', str(len(on)),end='\r')
        # get disc midplane density: array of size (nrad, nsec)
        dens = Field(field='dens', fluid='gas', on=on[k], directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units).data

        # total mass
        mass = np.sum(dens*surface)

        # get time
        mytime[k] = date[take_one_point_every*k]/2.0/np.pi  # orbital periods at apla

        # ---------------------
        # Fourier decomposition
        # ---------------------
        # real part of Fourier decomposition
        an = np.sum(dens*surface*np.cos(azi_wavenb*pmed2d)) / mass
        # an[m] = np.sum(dens*np.cos(azi_wavenb[m]*pmed2d)) / np.sum(dens)

        # imaginary part of Fourier decomposition
        bn = np.sum(dens*surface*np.sin(azi_wavenb*pmed2d)) / mass
        # bn[m] = np.sum(dens*np.sin(azi_wavenb[m]*pmed2d)) / np.sum(dens)

        # amplitude (the += arises when averaging over mutliple outputs)
        cn[k] = np.sqrt( an*an + bn*bn )


    # prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.20, right=0.96, top=0.95, bottom=0.12)
    ax = fig.gca()
    xtitle = r'Time [$T_0$]'
    ytitle = 'm = '+str(azi_wavenb)+' Fourier component of gas density'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')

    # set x-range
    if par.mytmin != '#':
        mytmin = par.mytmin
    else:
        mytmin = mytime[0]
    if par.mytmax != '#':
        mytmax = par.mytmax
    else:
        mytmax = mytime[-1]
    ax.set_xlim(mytmin,mytmax)

    # handle labels
    if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
        mylabel = str(par.use_legend)
    else:
        mylabel = str(par.directory)

    # ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.scatter(mytime, cn, color=par.c20[0], s=10, label=mylabel)

    # And save file
    outfile = 'fourier_m'+str(azi_wavenb)+'_'+str(par.directory)+'_'
    if np.isscalar(par.on) == False:
        outfile += str(par.on[0])+'_'+str(par.on[1])
    else:
        outfile += str(par.on)
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)



# Function that plots time evolution of either the maximum value through the disc 
# of the non-axisymmetric gas density, or the minimum of the Rossby number
def plotmaxnaodens_orminrossby():

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

    dens = Field(field='dens', fluid='gas', on=0, directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units)
    # get time
    if dens.fargo3d == 'Yes':
        f1, xpla, ypla, f4, f5, f6, f7, f8, date, omega = np.loadtxt(par.directory+"/planet0.dat",unpack=True)
    else:
        f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(par.directory+"/planet0.dat",unpack=True)

    mytime = np.zeros(len(on))
    quantity = np.zeros(len(on))

    fig = plt.figure(figsize=(8.,8.))

    # ========================
    # loop over output numbers
    # ========================
    for k in range(len(on)):

        print('output number =',str(k),'out of', str(len(on)),end='\r')
        # get disc midplane density: array of size (nrad, nsec)

        # get maximum value throught the disc of the quantity { Sigma - <Sigma> } / <Sigma>
        # where <Sigma> is the azimuthal-averaged radial profile of the gas density
        # maxnaodens[k] = dens.max()
        if par.plot_fourier == 'naodens':
            field = Field(field='dens', fluid='gas', on=on[k], directory=par.directory, physical_units=par.physical_units, nodiff='normnao', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units).data
            quantity[k] = field.max()
            ytitle = r'Max{$\Sigma / \langle\Sigma\rangle - 1$}'
            outfile = 'maxnaodens'+'_'+str(par.directory)+'_'

        # get minimum value throught the disc of the Rossby bymber
        if par.plot_fourier == 'rossby':
            field = Field(field='rossby', fluid='gas', on=on[k], directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units).data
            quantity[k] = field.min()
            ytitle = r'Min(Rossby number)'
            outfile = 'minrossby'+'_'+str(par.directory)+'_'

        # get time
        mytime[k] = date[take_one_point_every*k]/2.0/np.pi  # orbital periods at apla


    # prepare figure
    plt.subplots_adjust(left=0.18, right=0.96, top=0.95, bottom=0.12)
    ax = fig.gca()
    xtitle = r'Time [$T_0$]'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')

    # set x-range
    if par.mytmin != '#':
        mytmin = par.mytmin
    else:
        mytmin = mytime[0]
    if par.mytmax != '#':
        mytmax = par.mytmax
    else:
        mytmax = mytime[-1]
    ax.set_xlim(mytmin,mytmax)

    # handle labels
    if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
        mylabel = str(par.use_legend)
    else:
        mylabel = str(par.directory)

    # plot
    ax.scatter(mytime, quantity, color=par.c20[0], s=10, label=mylabel)

    # And save file
    if np.isscalar(par.on) == False:
        outfile += str(par.on[0])+'_'+str(par.on[1])
    else:
        outfile += str(par.on)
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)