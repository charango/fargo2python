import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import fnmatch
import re
import os

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

    # Range of azimuthal wavenumbers    
    m_min = 0 # 1
    m_max = 10 # int(dens.nsec/8)
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

        print('k = ', k, ' / ', len(on)-1 )
        # get disc midplane density: array of size (nrad, nsec)
        dens = Field(field='dens', fluid='gas', on=on[k], directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units).data

        # total mass
        mass = np.sum(dens*surface)

        # -------------------------------
        # loop over azimuthal wavenumbers
        # -------------------------------
        for m in range(len(azi_wavenb)):

            # real part of Fourier decomposition
            an[m] = np.sum(dens*surface*np.cos(azi_wavenb[m]*pmed2d)) / mass
            # an[m] = np.sum(dens*np.cos(azi_wavenb[m]*pmed2d)) / np.sum(dens)

            # imaginary part of Fourier decomposition
            bn[m] = np.sum(dens*surface*np.sin(azi_wavenb[m]*pmed2d)) / mass
            # bn[m] = np.sum(dens*np.sin(azi_wavenb[m]*pmed2d)) / np.sum(dens)

            # amplitude (the += arises when averaging over mutliple outputs)
            cn[m] += np.sqrt( an[m]*an[m] + bn[m]*bn[m] )

    # final amplitude - divide by len(on) in case of average over multiple output numbers
    for m in range(len(azi_wavenb)):
        cn[m] /= len(on)
        # print(azi_wavenb[m],cn[m])

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
    # ax.set_xscale('log')
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


def plotautocorrelationtimescale():

    # first import global variables
    import par

    # read tqwk0.dat file -> torque on a massless planet
    f1, it, ot, f4, f5, ip, op, f8, f9, time = np.loadtxt(par.directory+"/tqwk0.dat",unpack=True)
    tq = it+ot

    # time in orbital periods at R=1
    time /= (2.0*np.pi)

    tmax = 20 # time.max()
    nbtaustep = int(20.0*tmax)
    tau = np.zeros(nbtaustep)
    acf = np.zeros(nbtaustep)

    nbtimestep = len(tq)

    for k in range(1,nbtaustep): # tau goes from Torb/20 to tmax every Torb/20
        tau[k] = k/20.0
        num = 0.0
        den = 0.0

        for i in range(k,nbtimestep,1): # t goes from tau to TMAX every Torb/20
            num += tq[i]*tq[i-k]
            den += tq[i]*tq[i]

        acf[k] = num/den

    # prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.17, right=0.96, top=0.95, bottom=0.12)
    ax = fig.gca()
    xtitle = r'Lag [T$_{\rm orb}$]'
    ytitle = 'Auto-correlation function'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')

    # handle labels
    if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
        mylabel = str(par.use_legend)
    else:
        mylabel = str(par.directory)

    ax.set_xscale('log')
    ax.set_xlim(tau[1],tau[-1])

    ax.scatter(tau[1:], acf[1:], color=par.c20[0], s=10, label=mylabel)
    ax.plot(tau[1:], acf[1:], color=par.c20[0],linestyle='-')
    ax.plot(tau,0*tau,color=par.c20[0],linestyle='dotted')

    # And save file
    outfile = 'acf_'+str(par.directory)
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)


def plot_alphas():

    # first import global variables
    import par

    # get time-averaged density
    dim = len(fnmatch.filter(os.listdir(par.directory), 'summary*.dat'))
    on = range(0,dim,1)
    for i in range(len(on)):
        dens = Field(field='dens', fluid='gas', on=on[i], directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average='Yes', override_units=par.override_units)
        if i==0:
            nr = dens.nrad
            axidens = np.zeros(nr)
        axidens += np.sum(dens.data,axis=1)
    axidens /= len(on)
    axidens /= dens.nsec

    # get isothermal sound speed then pressure
    command = par.awk_command+' " /^ASPECTRATIO/ " '+par.directory+'/*.par'
    buf = subprocess.getoutput(command)
    aspectratio = float(buf.split()[1])
    command = par.awk_command+' " /^FLARINGINDEX/ " '+par.directory+'/*.par'
    buf = subprocess.getoutput(command)
    flaringindex = float(buf.split()[1])
    cs = aspectratio * dens.rmed**(flaringindex-0.5)  # isothermal sound speed (nrad)
    axipres = axidens*cs*cs

    # get 2D (r,time) binary file with Reynolds stress 
    f = par.directory+'/monitor/gas/reynolds_1d_Y_raw.dat'
    alpha_rey_file_data = np.fromfile(f, dtype='float64')

    # number of outputs in 
    nboutputs = int(len(alpha_rey_file_data)/nr)

    # reshape output as a 2D array
    buffer = alpha_rey_file_data.reshape(nboutputs,nr)  # 2D nrad, nb_outputs

    # time-average 1D radial profile of alpha_reynolds
    alpha_rey = np.sum(buffer, axis=0)/nboutputs/axipres

    # get 2D (r,time) binary file with Maxwell stress 
    f = par.directory+'/monitor/gas/maxwell_1d_Y_raw.dat'
    alpha_max_file_data = np.fromfile(f, dtype='float64')

    # number of outputs in 
    nboutputs = int(len(alpha_max_file_data)/nr)

    # reshape output as a 2D array
    buffer = alpha_max_file_data.reshape(nboutputs,nr)  # 2D nrad, nb_outputs

    # time-average 1D radial profile of alpha_reynolds
    alpha_max = -np.sum(buffer, axis=0)/nboutputs/axipres

    # prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.20, right=0.96, top=0.95, bottom=0.12)
    ax = fig.gca()
    xtitle = 'Radius'
    ytitle = 'Time-averaged alpha coefficients'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
    ax.set_xlim(dens.rmed.min(),dens.rmed.max())

    ax.plot(dens.rmed, alpha_rey, color=par.c20[0], label=r'$\alpha_{\rm Rey}$')
    ax.plot(dens.rmed, alpha_max, color=par.c20[1], label=r'$\alpha_{\rm Max}$')

    ax.legend(frameon=False,fontsize=15)

    # And save file
    outfile = 'alphas_'+str(par.directory)
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)    



# function that plots histogram of quantity (X - <X>) / <X> at different times
# where X is read as par.field
def plot_histofield():

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

    # prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.18, right=0.96, top=0.95, bottom=0.12)
    ax = fig.gca()
    
    ytitle = 'Histogram'
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')

    # handle labels
    if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
        mylabel = str(par.use_legend)
    else:
        mylabel = str(par.directory)

    # ax.set_yscale('log')
    #ax.set_xscale('log')

    if par.fieldmin != '#':
        min_bin = par.fieldmin
    else:
        min_bin = -0.3
    if par.fieldmax != '#':
        max_bin = par.fieldmax
    else:
        max_bin = 0.3
    nb_bins = 30
    mybins = min_bin + (max_bin-min_bin)*np.arange(nb_bins)/(nb_bins-1.0)

    myfield = Field(field='dens', fluid='gas', on=0, directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units)

    nrad = myfield.nrad
    nsec = myfield.nsec
    rmed = myfield.rmed
    myrmin = 1.1*rmed.min()
    imin = np.argmin(np.abs(rmed-myrmin))
    myrmax = 0.9*rmed.max()
    imax = np.argmin(np.abs(rmed-myrmax))

    # ========================
    # loop over output numbers
    # ========================
    for k in range(len(on)):

        print('k = ', k, ' / ', len(on)-1 )

        # get disc midplane X field: array of size (nrad, nsec)
        buf = Field(field=par.whatfield, fluid=par.fluid, on=on[k], directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units)
        axifield = (np.sum(buf.data ,axis=1)/nsec).repeat(nsec).reshape(nrad,nsec)
        myfield = buf.data-axifield

        if par.whatfield == 'dens' or par.whatfield == 'vtheta':
            myfield = (buf.data-axifield)/axifield

        myfieldnew = myfield[imin:imax,:]
        myfieldoned = myfieldnew.reshape((imax-imin)*nsec)

        # # add histogram here
        cmap = matplotlib.cm.get_cmap('Spectral_r')
        if len(on) > 1:
            c20 = cmap(k/(len(on)-1.0))
        else:
            c20 = 'tab:blue'

        counts, bins, patches = plt.hist(myfieldoned, bins=mybins, color=c20, alpha=0.3, rwidth=0.9, density=True)
        if k==0:
            ax.set_ylim(0,2.5*counts.max())
            ax.set_xlim(bins.min(),bins.max())
            if par.whatfield == 'dens':
                xtitle = r'$(\Sigma - \langle \Sigma\rangle_\varphi) / \langle \Sigma\rangle_\varphi$'
                outfile = 'histodens_'
            if par.whatfield == 'vrad':
                xtitle = r'$v_r - \langle v_r \rangle_\varphi$'
                outfile = 'histovrad_'
            if par.whatfield == 'vtheta':
                xtitle = r'$(v_{\phi} - \langle v_{\phi} \rangle_\varphi) / \langle v_{\phi} \rangle_\varphi$'
                outfile = 'histovtheta_'
            ax.set_xlabel(xtitle)

    # And save file
    outfile = outfile+str(par.directory)+'_'
    if np.isscalar(par.on) == False:
        outfile += str(par.on[0])+'_'+str(par.on[1])
    else:
        outfile += str(par.on)
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)


# function that plots histogram of the specific torque
def plot_histotorque():

    # first import global variables
    import par

    # prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.16, right=0.96, top=0.95, bottom=0.12)
    ax = fig.gca()
    xtitle = 'Specific turbulent torque'
    ytitle = 'Histogram'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')

    # now, read tqwk0.dat file
    f1, it, ot, f4, f5, ip, op, f8, f9, time = np.loadtxt(par.directory+"/tqwk0.dat",unpack=True)
    tq = it+ot

    if par.myymin != '#':
        min_bin = par.myymin
    else:
        min_bin = tq.min()
    if par.myymax != '#':
        max_bin = par.myymax
    else:
        max_bin = tq.max()
    nb_bins = 30
    mybins = min_bin + (max_bin-min_bin)*np.arange(nb_bins)/(nb_bins-1.0)

    # plot histogram
    n, bins, patches = plt.hist(x=tq, bins=mybins, color=par.c20[0], alpha=1.0, rwidth=0.9)

    # And save file
    fileout = 'histotorque_'+str(par.directory)+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)


# function that plots time evolution of the disc-averaged Reynolds alpha parameter
def plot_time_alphaRey():

    # first import global variables
    import par

    # Define range of output numbers to consider
    if par.take_one_point_every == '#':
        take_one_point_every = 1
    else:
        take_one_point_every = par.take_one_point_every

    if par.on == 'all':
        import fnmatch
        if isinstance(par.directory, list) == True:
            dir = par.directory[0]
        else:
            dir = par.directory
        if par.fargo3d == 'No':
            nboutputs = len(fnmatch.filter(os.listdir(dir), 'gasdens*.dat'))-len(fnmatch.filter(os.listdir(dir), 'gasdens.ascii*.dat'))
            if par.fargo2d1d == 'Yes':
                nboutputs = len(fnmatch.filter(os.listdir(dir), 'gasdens1D*.dat'))
        else:
            nboutputs = len(fnmatch.filter(os.listdir(dir), 'summary*.dat'))
        # on = [0,nboutputs-1]
        on = range(0,nboutputs-1,take_one_point_every)
    else:
        if np.isscalar(par.on) == False:
            on = range(par.on[0],par.on[1]+1,take_one_point_every)
        else:
            on = [par.on]
        
    # 2D arrays with radius and azimuth
    dens = Field(field='dens', fluid='gas', on=0, directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units)

    # get time
    if dens.fargo3d == 'Yes':
        f1, xpla, ypla, f4, f5, f6, f7, f8, date, omega = np.loadtxt(par.directory+"/planet0.dat",unpack=True)
    else:
        f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(par.directory+"/planet0.dat",unpack=True)

    alpharey = np.zeros(len(on))
    mytime = np.zeros(len(on))

    # ========================
    # loop over output numbers
    # ========================
    for k in range(len(on)):

        print('output number =',str(k+1),'out of', str(len(on)),end='\r')

        # get time
        mytime[k] = date[take_one_point_every*k]/2.0/np.pi  # orbital periods at R=1

        # get 2D field of Reynolds alpha parameter
        buf = Field(field='alpha_reynolds', fluid='gas', on=on[k], directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units).data

        # Option 1: take the mean of alpha_Rey over the disc
        alpharey[k] = np.mean(buf)
        # Option 1: take the median of alpha_Rey over the disc
        # alpharey[k] = np.median(buf)


    # option to write result in 1D ascii file
    if ( ('write_ascii' in open('paramsf2p.dat').read()) and (par.write_ascii == 'Yes') ):
        ascii = open('timealpharey_'+par.directory+'.dat','w')
        for v in range(len(mytime)):
            ascii.write(str(mytime[v])+'\t'+str(alpharey[v])+'\n')


    # prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.16, right=0.96, top=0.95, bottom=0.12)
    ax = fig.gca()
    xtitle = r'Time [$T_0$]'
    ytitle = 'Disc-averaged Reynolds alpha parameter'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')

    # set x-range
    if par.mytmin != '#':
        mytmin = par.mytmin
    else:
        mytmin = 0.0
    if par.mytmax != '#':
        mytmax = par.mytmax
    else:
        mytmax = mytime.max()
    ax.set_xlim(mytmin,mytmax)

    # set x-range
    if par.myymin != '#':
        myymin = par.myymin
    else:
        myymin = alpharey.min()
    if par.myymax != '#':
        myymax = par.myymax
    else:
        myymax = alpharey.max()
    ax.set_ylim(myymin,myymax)

    # plot
    ax.plot(mytime, alpharey, color=par.c20[1])

    # And save file
    outfile = 'timealpharey_'+str(par.directory)
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)    



# function that plots the time evolution of the disc-averaged gravitational stress tensor (alpha form) 
# obtained from the radial and azimuthal components of the self-gravitating acceleration (which need 
# to be output by your FARGO-ADSG run)
def plot_time_alphaGrav():

    # first import global variables
    import par

    # Define range of output numbers to consider
    if par.take_one_point_every == '#':
        take_one_point_every = 1
    else:
        take_one_point_every = par.take_one_point_every

    if par.on == 'all':
        import fnmatch
        if isinstance(par.directory, list) == True:
            dir = par.directory[0]
        else:
            dir = par.directory
        nboutputs = len(fnmatch.filter(os.listdir(dir), 'gasdens*.dat'))-len(fnmatch.filter(os.listdir(dir), 'gasdens.ascii*.dat'))
        # on = [0,nboutputs-1]
        on = range(0,nboutputs-1,take_one_point_every)
    else:
        if np.isscalar(par.on) == False:
            on = range(par.on[0],par.on[1]+1,take_one_point_every)
        else:
            on = [par.on]

    # get isothermal sound speed
    command = par.awk_command+' " /^AspectRatio/ " '+par.directory+'/*.par'
    buf = subprocess.getoutput(command)
    aspectratio = float(buf.split()[1])
    command = par.awk_command+' " /^FlaringIndex/ " '+par.directory+'/*.par'
    buf = subprocess.getoutput(command)
    flaringindex = float(buf.split()[1])
    # check if energy equation was used and then get adiabatic index
    command = par.awk_command+' " /^EnergyEquation/ " '+par.directory+'/*.par'
    buf = subprocess.getoutput(command)
    energyequation = str(buf.split()[1])
    command = par.awk_command+' " /^AdiabaticIndex/ " '+par.directory+'/*.par'
    buf = subprocess.getoutput(command)
    gamma = float(buf.split()[1])
        
    # Read initial density to inherit nrad, nsec, rmed...
    dens = Field(field='dens', fluid='gas', on=0, directory=par.directory, physical_units=par.physical_units, nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units)
    nrad = dens.nrad
    nsec = dens.nsec
    rmed = dens.rmed

    # get time
    f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(par.directory+"/planet0.dat",unpack=True)

    alphagrav = np.zeros(len(on))
    mytime = np.zeros(len(on))

    # ========================
    # loop over output numbers
    # ========================
    for k in range(len(on)):

        print('output number =',str(k+1),'out of', str(len(on)),end='\r')

        # get time
        mytime[k] = date[take_one_point_every*k]/2.0/np.pi  # orbital periods at R=1

        # get radial and azimuthal components of self-gravitating acceleration
        gr = Field(field='sgaccr', fluid='gas', on=on[k], directory=par.directory, physical_units='No', nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units).data
        gphi = Field(field='sgacctheta', fluid='gas', on=on[k], directory=par.directory, physical_units='No', nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units).data

        # azimuthally-averaged pressure
        dens = Field(field='dens', fluid='gas', on=on[k], directory=par.directory, physical_units='No', nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units).data
        if energyequation == 'Yes':
            temp = Field(field='temp', fluid='gas', on=on[k], directory=par.directory, physical_units='No', nodiff='Yes', fieldofview=par.fieldofview, onedprofile='No', slice='midplane', z_average=par.z_average, override_units=par.override_units).data
            pressure = gamma*dens*temp
            axipres = np.sum(pressure,axis=1)/nsec  # azimuthally-averaged pressure (nrad)
        else:
            cs = aspectratio * rmed**(flaringindex-0.5)  # adiabatic sound speed at t=0 (nrad)!
            pressure = dens*((cs*cs).repeat(nsec).reshape(nrad,nsec))  # 2D thermal pressure
            axipres = np.sum(pressure,axis=1)/nsec  # azimuthally-averaged pressure (nrad)

        # Gravitational stress
        H = aspectratio*(rmed**(1+flaringindex))
        H2D = H.repeat(nsec).reshape(nrad,nsec)
        stress = gr*gphi*2.0*H2D/4.0/np.pi   # (nrad,nsec)
        stress *= 2.0 # CUIDADIN (test)

        # radial profile of alphagrav = (2/3) x <TRphi> / <pressure>
        # with <.> = azimuthal average
        alphagrav_R = (2.0/3.0) * np.sum(stress,axis=1) / axipres / nsec

        # final alpha = radial mean of alpha radial profile
        alphagrav[k] = np.mean(alphagrav_R)


    # option to write result in 1D ascii file
    if ( ('write_ascii' in open('paramsf2p.dat').read()) and (par.write_ascii == 'Yes') ):
        ascii = open('timealphagrav_'+par.directory+'.dat','w')
        for v in range(len(mytime)):
            ascii.write(str(mytime[v])+'\t'+str(alphagrav[v])+'\n')


    # prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.16, right=0.96, top=0.95, bottom=0.12)
    ax = fig.gca()
    xtitle = r'Time [$T_0$]'
    ytitle = 'Disc-averaged gravitational alpha parameter'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')

    # set x-range
    if par.mytmin != '#':
        mytmin = par.mytmin
    else:
        mytmin = 0.0
    if par.mytmax != '#':
        mytmax = par.mytmax
    else:
        mytmax = mytime.max()
    ax.set_xlim(mytmin,mytmax)

    # set x-range
    if par.myymin != '#':
        myymin = par.myymin
    else:
        myymin = alphagrav.min()
    if par.myymax != '#':
        myymax = par.myymax
    else:
        myymax = alphagrav.max()
    ax.set_ylim(myymin,myymax)

    # plot
    ax.plot(mytime, alphagrav, color=par.c20[1])

    # And save file
    outfile = 'timealphagrav_'+str(par.directory)
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)  


# function that plots the time evolution of the disc-averaged gravitational stress tensor 
# calculated by external C program alphasg.c (that used in Baruteau+ 2011, see appendix)
def plot_time_alphaGrav_external():

    # first import global variables
    import par

    # get all files alphasg_gl_X.dat where X = output number
    alphafiles = fnmatch.filter(os.listdir('.'), 'alphasg*.dat')
    nboutputs = len(alphafiles)

    # get output number
    on = np.zeros(nboutputs)
    for i in range(nboutputs):
        on[i] = alphafiles[i][11:-4]

    # integer array with output numbers
    inton = [int(x) for x in on]
    sortedon = [0] * nboutputs
    sortedalphafiles = ["" for x in range(nboutputs)]

    for i in range(nboutputs):
        index = np.argsort(inton)[i]
        sortedalphafiles[i] = alphafiles[index]
        sortedon[i] = int(inton[index])

    # get time
    f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(par.directory+"/planet0.dat",unpack=True)

    alphagrav = np.zeros(nboutputs)
    mytime = np.zeros(nboutputs)

    # ========================
    # loop over output numbers
    # ========================
    for k in range(nboutputs):

        # get time
        mytime[k] = date[sortedon[k]]/2.0/np.pi  # orbital periods at R=1

        # read alphasg_glX.dat file
        r, a = np.loadtxt(sortedalphafiles[k],unpack=True)

        alphagrav[k] = np.mean(a)


    # option to write result in 1D ascii file
    if ( ('write_ascii' in open('paramsf2p.dat').read()) and (par.write_ascii == 'Yes') ):
        ascii = open('timealphagravext_'+par.directory+'.dat','w')
        for v in range(len(mytime)):
            ascii.write(str(mytime[v])+'\t'+str(alphagrav[v])+'\n')


    # prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.16, right=0.96, top=0.95, bottom=0.12)
    ax = fig.gca()
    xtitle = r'Time [$T_0$]'
    ytitle = 'Disc-averaged gravitational alpha parameter'
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')

    # set x-range
    if par.mytmin != '#':
        mytmin = par.mytmin
    else:
        mytmin = 0.0
    if par.mytmax != '#':
        mytmax = par.mytmax
    else:
        mytmax = mytime.max()
    ax.set_xlim(mytmin,mytmax)

    # set x-range
    if par.myymin != '#':
        myymin = par.myymin
    else:
        myymin = alphagrav.min()
    if par.myymax != '#':
        myymax = par.myymax
    else:
        myymax = alphagrav.max()
    ax.set_ylim(myymin,myymax)

    # plot
    ax.plot(mytime, alphagrav, color=par.c20[1])

    # And save file
    outfile = 'timealphagravext_'+str(par.directory)
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)    