import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import fnmatch
import os
import re

from mesh import *
from field import *

def plotdisccom():

    # first import global variables
    import par
    
    # first prepare figure
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.18, right=0.94, top=0.94, bottom=0.12)
    ax = fig.gca()
    
    if par.plot_disccom == 'xy':
        xtitle = 'centre-of-mass x'
        ytitle = 'centre-of-mass y'
    if par.plot_disccom == 'tr' or par.plot_disccom == 'tlogr' or par.plot_disccom == 'logtlogr':
        xtitle = 'time [orbits]'
        ytitle = 'centre-of-mass radius'
        if par.physical_units == 'Yes':
            ytitle += ' [au]'
    if par.plot_disccom == 'tpa':
        xtitle = 'time [orbits]'
        ytitle = 'Centre-of-mass position angle [rad]'

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

        # if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
        #     if len(directory) == 1:
        #         mylabel = str(par.use_legend)
        #     else:
        #         mylabel = str(par.use_legend[j])
        # else:
        #     mylabel = str(directory[j])

        if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
            use_legend = par.use_legend
            if isinstance(par.use_legend, str) == True:
                use_legend = [par.use_legend]
            mylabel = str(use_legend[j])
            mylabel = mylabel.replace("_", " ")
        else:
            mylabel = str(directory[j])

        # Test if file gasdens1D0.dat exists to know if simulation in current directory 
        # has been carried out with FARGO2D1D or not
        gasdens1D0_file  = directory[j]+'/gasdens1D0.dat'
        if os.path.isfile(gasdens1D0_file) == True:
            fargo2d1d = 'Yes'
        else:
            fargo2d1d = 'No'

        # Test if simulation has been done with FARGO3D in 3D
        runwas3d = 'No'
        summary_file  = directory[j]+'/summary0.dat'
        if os.path.isfile(summary_file) == True:
            command = par.awk_command+' " /^ZMAX/ " '+directory[j]+'/variables.par'
            if sys.version_info[0] < 3:   # python 2.X
                buf = subprocess.check_output(command, shell=True)
            else:                         # python 3.X
                buf = subprocess.getoutput(command)
            zmax = float(buf.split()[1])
            if zmax != 1:
                runwas3d = 'Yes'


        # DEFAULT CASE (= 2D simulations): we obtain the position of the disc's center-of-mass 
        # by inspecting at the gas density fields. When the code is not FARGO-2D1D, the disc's 
        # center-of-mass is plotted in the reference frame centred on the star. Otherwise, we 
        # plot the distance between the star and the disc's center-of-mass (well, it is actually
        # the same quantity regardless of the version of FARGO code used!)
        if runwas3d == 'No':

            # find how many output numbers were produced for each directory
            if par.fargo3d == 'No':
                if fargo2d1d == 'No':
                    nboutputs = len(fnmatch.filter(os.listdir(directory[j]), 'gasdens*.dat'))
                else:
                    nboutputs = len(fnmatch.filter(os.listdir(directory[j]), 'gasdens1D*.dat'))
            else:
                nboutputs = len(fnmatch.filter(os.listdir(directory[j]), 'summary*.dat'))
            print('number of outputs for directory ',directory[j],': ',nboutputs)

            on = range(0,nboutputs-1,take_one_point_every)

            x_com = np.zeros(len(on))
            y_com = np.zeros(len(on))
            r_com = np.zeros(len(on))
            t_com = np.zeros(len(on))
            pa_com = np.zeros(len(on))
            
            first_time = 0

            # loop over output numbers
            for k in range(len(on)):     

                print('disc com: output number '+str(k)+' / '+str(len(on)-1),end='\r')

                # get 2D gas surface density field (not compatible with 3D yet...)
                dens = Field(field='dens', fluid='gas', on=int(on[k]), directory=directory[j], physical_units=par.physical_units, nodiff='Yes', fieldofview='polar', slice='midplane', onedprofile='No', override_units=par.override_units)

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

                    # get radius and azimuth arrays, infer X and Y for cell centres
                    R = dens.rmed
                    if par.physical_units == 'Yes':
                        R *= (dens.culength / 1.5e11) # in au
                    
                    T = dens.pmed
                    radius_matrix, theta_matrix = np.meshgrid(R,T)
                    X = radius_matrix * np.cos(theta_matrix)       # nsec, nrad
                    Y = radius_matrix * np.sin(theta_matrix)       # nsec, nrad

                    # get time
                    index_orbit_file = 0
                    if dens.fargo3d == 'Yes':
                        f1, xpla, ypla, f4, f5, f6, f7, mpla, date, f10 = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
                    else:
                        if fargo2d1d == 'Yes':
                            f1, xpla, ypla, f4, f5, mpla, f7, date, f9 = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
                            index_orbit_file = 1
                        else:
                            f1, xpla, ypla, f4, f5, mpla, f7, date, f9, f10, f11 = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
                    with open(directory[j]+"/orbit"+str(index_orbit_file)+".dat") as f_in:
                        firstline_orbitfile = np.genfromtxt(itertools.islice(f_in, 0, 1, None), dtype=float)
                    apla = firstline_orbitfile[2]
                

                # mass of each grid cell (2D array)
                mass = dens.data*surface   
                mass = np.transpose(mass)  # (nsec,nrad)

                # Fargo2D1D case: get star's mass and position
                if fargo2d1d == 'Yes':
                    mstar = mpla[int(on[k])]
                    xstar = xpla[int(on[k])]
                    ystar = ypla[int(on[k])]
                else:
                    mstar = 1.0
                    xstar = 0.0
                    ystar = 0.0

                if fargo2d1d == 'No':
                    # x_com and y_com are the coordinates of the centre-of-mass 
                    # of the {disc+star} system in the stellocentric frame
                    x_com[k] = np.sum(mass*X) / (mstar + np.sum(mass))
                    y_com[k] = np.sum(mass*Y) / (mstar + np.sum(mass))
                else:
                    # FARGO2D1D: x_com and y_com are ALSO the coordinates of the centre-of-mass 
                    # of the {disc+star} system in the stellocentric frame. In practice here 
                    # x_com ~ -xstar and y_com ~ -ystar (checked)
                    x_com[k] = ((mstar*xstar + np.sum(mass*X)) / (mstar + np.sum(mass))) - xstar
                    y_com[k] = ((mstar*ystar + np.sum(mass*Y)) / (mstar + np.sum(mass))) - ystar

                # r_com is the distance between the star and the centre-of-mass of the {star+disc} system
                r_com[k] = np.sqrt( x_com[k]*x_com[k] + y_com[k]*y_com[k] )
                # time
                t_com[k] = round(date[k*take_one_point_every]/2./np.pi/apla/np.sqrt(apla),1)

                # For simulations of circumbinary discs, return X and Y of disc's centre of mass
                # to then get its position angle
                if par.central_binary == 'Yes':
                    x_com[k] = np.sum(mass*X) / np.sum(mass)
                    y_com[k] = np.sum(mass*Y) / np.sum(mass)
                    pa_com[k] = math.atan2(y_com[k],x_com[k])+np.pi # + pi to be consistent with simulations outputs
                    t_com[k] = round(date[k*take_one_point_every]/2./np.pi,1)

                #print(mp*xp, np.sum(mass*X), mp*xp+np.sum(mass*X))
                #print(mp*yp, np.sum(mass*Y), mp*yp+np.sum(mass*Y))
                #print(x_com[k],y_com[k],r_com[k])


        # CASE OF 3D SIMULATIONS WITH FARGO3D: we obtain again the position of the center-of-mass 
        # by inspecting at the 3D gas density fields obtained in simulations run in a fixed reference 
        # frame centred on the star
        if runwas3d == 'Yes':
            print('3D FARGO3D simulation detected!')

            # find how many output numbers were produced for each directory
            nboutputs = len(fnmatch.filter(os.listdir(directory[j]), 'summary*.dat'))
            print('number of outputs for directory ',directory[j],': ',nboutputs)

            on = range(0,nboutputs-1,take_one_point_every)

            x_com = np.zeros(len(on))
            y_com = np.zeros(len(on))
            z_com = np.zeros(len(on))
            r_com = np.zeros(len(on))
            t_com = np.zeros(len(on))
            pa_com = np.zeros(len(on))

            # get 2D gas surface density field just to inherit from mesh properties
            buf = Field(field='dens', fluid='gas', on=int(on[0]), directory=directory[j], physical_units=par.physical_units, nodiff='Yes', fieldofview='polar', slice='midplane', onedprofile='No', override_units=par.override_units)

            first_time = 0

            # loop over output numbers
            for k in range(len(on)):     

                print('disc com: output number '+str(k)+' / '+str(len(on)-1),end='\r')

                # get 3D gas volume density field
                f = directory[j]+'/gasdens'+str(int(on[k]))+'.dat'
                #print('file = ', f)
                dens = np.fromfile(f, dtype='float64')
                dens = dens.reshape(buf.nz,buf.nrad,buf.nsec)  # 3D nz, nrad, nsec

                # things we do only when entering for loop
                if first_time == 0:
                
                    first_time = 1
                
                    # get volume of every cell
                    Redge,Cedge,Aedge = np.meshgrid(buf.redge, buf.tedge, buf.pedge)   # nz+1, nrad+1, Nsec+1
                    if par.physical_units == 'Yes':
                        Redge *= (buf.culength / 1.5e11) # in au
                    r2 = Redge*Redge

                    jacob  = r2[:-1,:-1,:-1] * np.cos(Cedge[:-1,:-1,:-1])
                    dphi   = Aedge[:-1,:-1,1:] - Aedge[:-1,:-1,:-1]     # same as 2pi/nsec
                    dr     = Redge[:-1,1:,:-1] - Redge[:-1,:-1,:-1]     # same as Rsup-Rinf
                    dtheta = Cedge[1:,:-1,:-1] - Cedge[:-1,:-1,:-1]
                    volume = jacob * dr * dphi * dtheta     # nz, nrad, nsec

                    # get radius and azimuth arrays, infer X and Y for cell centres
                    R = buf.rmed
                    if par.physical_units == 'Yes':
                        R *= (buf.culength / 1.5e11) # in au
                    P = buf.pmed
                    T = buf.tmed

                    radius_matrix, phi_matrix, theta_matrix = np.meshgrid(R,P,T, indexing='ij')   # nrad nsec nz
                    X = radius_matrix * np.cos(theta_matrix) * np.cos(phi_matrix)  
                    Y = radius_matrix * np.cos(theta_matrix) * np.sin(phi_matrix)
                    Z = radius_matrix * np.sin(theta_matrix)

                    # get time
                    f1, xpla, ypla, zpla, f5, f6, f7, mpla, date, f10 = np.loadtxt(directory[j]+"/planet0.dat",unpack=True)
                    with open(directory[j]+"/orbit0.dat") as f_in:
                        firstline_orbitfile = np.genfromtxt(itertools.islice(f_in, 0, 1, None), dtype=float)
                    apla = firstline_orbitfile[2]
                

                # mass of each grid cell (3D array)
                mass = dens*volume            # nz, nrad, nsec
                mass = np.swapaxes(mass,1,2)  # nz, nsec, nrad
                mass = np.transpose(mass)     # nrad, nsec, nz
                #print('np.sum(mass) = ', np.sum(mass))

                '''
                # is there a planet?
                mp = mpla[int(on[k])]
                xp = xpla[int(on[k])]
                yp = ypla[int(on[k])]
                zp = zpla[int(on[k])]
                '''
                mstar = 1.0
                xstar = 0.0
                ystar = 0.0
                zstar = 0.0

                # get x-, y- and z-coordinates of centre-of-mass
                x_com[k] = np.sum(mass*X) / (mstar + np.sum(mass))
                y_com[k] = np.sum(mass*Y) / (mstar + np.sum(mass))
                z_com[k] = np.sum(mass*Z) / (mstar + np.sum(mass))
                r_com[k] = np.sqrt( x_com[k]*x_com[k] + y_com[k]*y_com[k] + z_com[k]*z_com[k] )
                # test:
                #r_com[k] = np.sqrt( x_com[k]*x_com[k] + y_com[k]*y_com[k] )
                #print('xcom, ycom, zcom, rcom = ', x_com[k], y_com[k], z_com[k], r_com[k])
                # time
                t_com[k] = round(date[k*take_one_point_every]/2./np.pi/apla/np.sqrt(apla),1)

                # For simulations of circumbinary discs, return X and Y of disc's centre of mass
                if par.central_binary == 'Yes':
                    x_com[k] = np.sum(mass*X) / np.sum(mass)
                    y_com[k] = np.sum(mass*Y) / np.sum(mass)
                    t_com[k] = round(date[k*take_one_point_every]/2./np.pi,1)
                    pa_com[k] = math.atan2(y_com[k],x_com[k])+np.pi # + pi to be consistent with simulations outputs


        # find minimum anx maximum time over directories
        if par.plot_disccom == 'tr':
            mytmin = np.minimum(mytmin,t_com[0])
        else:
            mytmin = np.minimum(mytmin,t_com[1])
        mytmax = np.maximum(mytmax,t_com.max())

        # display data as scatter plot for each directory
        if par.plot_disccom == 'xy':
            ax.scatter(x_com, y_com, s=30, marker='+', c=t_com, cmap='nipy_spectral', alpha=1.0)
        if par.plot_disccom == 'tr':
            ax.scatter(t_com, r_com, s=30, marker='+', alpha=1.0, color=par.c20[j],label=mylabel)
        if par.plot_disccom == 'tlogr':
            ax.set_yscale('log')
            ax.scatter(t_com[1:len(t_com)-1], r_com[1:len(t_com)-1], s=20, marker='+', alpha=1.0, color=par.c20[j],label=mylabel)
        if par.plot_disccom == 'logtlogr':
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.scatter(t_com[1:len(t_com)-1], r_com[1:len(t_com)-1], s=20, marker='+', alpha=1.0, color=par.c20[j],label=mylabel)
        if par.central_binary == 'Yes' and par.plot_disccom == 'tpa':
            ax.scatter(t_com, pa_com, s=30, marker='+', alpha=1.0, color=par.c20[j],label=mylabel)
            ax.set_ylim(0.0,2.0*np.pi)

        # save data in ascii file
        if par.central_binary == 'Yes' and par.plot_disccom == 'tpa':
            fileout = open('pa_com_'+str(directory[j])+'_every'+str(take_one_point_every)+'.dat','w')
            for i in range(1,len(t_com)):
                fileout.write(str(t_com[i])+'\t'+str(pa_com[i])+'\n')
            fileout.close()
        else:
            fileout = open('log10rcom_'+str(directory[j])+'.dat','w')
            fileout.write('# time[orb]\t log10(rcom)\n')
            for i in range(1,len(t_com)):
                fileout.write(str(t_com[i])+'\t'+str(np.log10(r_com[i]))+'\n')
            fileout.close()
            
    # set x-range
    if par.plot_disccom != 'xy':
        if par.mytmin != '#':
            mytmin = par.mytmin
        if par.mytmax != '#':
            mytmax = par.mytmax
        ax.set_xlim(mytmin,mytmax)
    else:
        if ( ('myymin' in open('paramsf2p.dat').read()) and (par.myymin != '#') and (('myymax' in open('paramsf2p.dat').read()) and (par.myymax != '#')) ):
            ymin = par.myymin
            ymax = par.myymax
            ax.set_ylim(ymin,ymax)
            ax.set_xlim(ymin,ymax)

    # set y-range
    if ( ('myymin' in open('paramsf2p.dat').read()) and (par.myymin != '#') and (('myymax' in open('paramsf2p.dat').read()) and (par.myymax != '#')) ):
        ymin = par.myymin
        ymax = par.myymax
        ax.set_ylim(ymin,ymax)
    
    # finally add legend
    ax.set_axisbelow(False)
    ax.grid(axis='both', which='major', ls='-', alpha=0.8)

    if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != 'None'):
        legend = plt.legend(loc='lower right',fontsize=15,facecolor='white',edgecolor='white',framealpha=0.85,numpoints=1,bbox_transform=plt.gcf().transFigure)
        for line, text in zip(legend.get_lines(), legend.get_texts()):
            text.set_color(line.get_color())

    # And save file
    if len(directory) == 1:           
       outfile = 'com_'+par.plot_disccom+'_'+str(directory[0])
    else:
       outfile = 'com_'+par.plot_disccom
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=160)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)
