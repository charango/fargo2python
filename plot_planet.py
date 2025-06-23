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
    if par.plot_planet[1] == 'adot':
        ytitle = 'da/dt'
    if par.plot_planet[1] == 'r':
        ytitle = 'Orbital radius'
        if par.physical_units == 'Yes':
            ytitle += ' [au]'
    if par.plot_planet[1] == 'rdot':
        ytitle = 'dr/dt'
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

    if ( ('myymax' in open('paramsf2p.dat').read()) and ('myymin' in open('paramsf2p.dat').read()) ):
        if ( (par.myymax != '#') and (par.myymin != '#') ):
            ax.set_ylim(par.myymin,par.myymax)

    # loop over directories
    for j in range(len(directory)):

        # Locally check with which version of FARGO each simulation was run
        fargo3d = 'No'
        fargo_orig = 'Yes'
        fargo2d1d = 'No'
        if os.path.isfile(directory[j]+'/summary0.dat') == True:
            # Simulation was carried out with Fargo3D
            fargo3d = 'Yes'
            fargo_orig = 'No'
        else:
            if os.path.isfile(directory[j]+'/used_azi.dat') == True:
            # Simulation was carried out with Dusty FARGO-ADSG
                fargo_orig = 'No'
            else:
                # Simulations were carried out with the original FARGO code or with FARGO-2D1D
                if os.path.isfile(directory[j]+'/gasdens1D0.dat') == True:
                    fargo2d1d = 'Yes'
                    fargo_orig = 'Yes'
                

        # work out physical units
        if par.physical_units == 'Yes':
            if fargo3d == 'No' and fargo_orig == 'No':
            # units.dat contains physical units of mass [kg], length [m], time [s], and temperature [k] 
                cumass, culength, cutime, cutemp = np.loadtxt(directory[j]+"/units.dat",unpack=True)
            if fargo3d == 'Yes':
            # get units via variable.par file
                command = par.awk_command+' " /^UNITOFLENGTHAU/ " '+directory[j]+'/variables.par'
                # check which version of python we're using
                if sys.version_info[0] < 3:   # python 2.X
                    buf = subprocess.check_output(command, shell=True)
                else:                         # python 3.X
                    buf = subprocess.getoutput(command)
                if buf:
                    culength = float(buf.split()[1])*1.5e11  #from au to meters
                else:
                    print('UNITOFLENGTHAU is absent in variables.par, I will assume it is equal to unity, meaning that your code unit of length was 1 au!')
                    culength = 1.5e11 # 1 au in meters
                command = par.awk_command+' " /^UNITOFMASSMSUN/ " '+directory[j]+'/variables.par'
                # check which version of python we're using
                if sys.version_info[0] < 3:   # python 2.X
                    buf = subprocess.check_output(command, shell=True)
                else:                         # python 3.X
                    buf = subprocess.getoutput(command)
                if buf:
                    cumass = float(buf.split()[1])*2e30  #from Msol to kg
                else:
                    print('UNITOFMASSMSUN is absent in variables.par, I will assume it is equal to unity, meaning that your code unit of mass was 1 Solar mass!')
                    cumass = 2.0e30 # 1 Solar mass in kg
                # unit of time = sqrt( pow(L,3.) / 6.673e-11 / M );
                cutime = np.sqrt( culength**3.0 / 6.673e-11 / cumass)
            # case we overrride code units
            if par.override_units == 'Yes':
                if par.new_unit_length == 0.0:
                    sys.exit('override_units set to yes but new_unit_length is not defined in params.dat, I must exit!')
                else:
                    if par.verbose == 'Yes':
                        print('new unit of length in meters : ', par.new_unit_length)
                if par.new_unit_mass == 0.0:
                    sys.exit('override_units set to yes but new_unit_mass is not defined in params.dat, I must exit!')
                else:
                    if par.verbose == 'Yes':
                        print('new unit of mass in kg : ', par.new_unit_mass)
                cumass = par.new_unit_mass
                culength = par.new_unit_length
                # Deduce new units of time and temperature:
                # T = sqrt( pow(L,3.) / 6.673e-11 / M )
                # U = mmw * 8.0841643e-15 * M / L;
                cutime = np.sqrt( culength**3 / 6.673e-11 / cumass)
                cutemp = 2.35 * 8.0841643e-15 * cumass / culength
                if par.verbose == 'Yes':
                    print('### NEW UNITS SPECIFIED: ###')
                    print('new unit of length [m] = ',culength)
                    print('new unit of mass [kg]  = ',cumass)
                    print('new unit of time [s] = ',cutime)
                    print('new unit of temperature [K] = ',cutemp)
                
        if fargo3d == 'No':
            if fargo_orig == 'No':
                time, e, a, M, V, PA, mylambda, varpi = np.loadtxt(directory[j]+"/orbit0.dat",unpack=True)
            else:
                if fargo2d1d == 'Yes':
                    time, e, a, M, V, PA = np.loadtxt(directory[j]+"/orbit1.dat",unpack=True)
                else:
                    time, e, a, M, V, PA = np.loadtxt(directory[j]+"/orbit0.dat",unpack=True)
        else:
            command = par.awk_command+' "{print NF; exit}" '+directory[j]+'/orbit0.dat'
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

        # CUIDADIN!!
        if fargo2d1d == 'Yes':
            nbplanets = 1

        if (nbplanets <= 1 and par.plot_planet[1] == 'p'):
            sys.exit('ERROR: you requested to plot orbital period ratio, but there is only one planet simulated in directory: ', directory[j])
           
        if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
            use_legend = par.use_legend
            #print('use_legend = ', use_legend)
            if isinstance(par.use_legend, str) == True:
                use_legend = [par.use_legend]
            mylabel = str(use_legend[j])
            mylabel = mylabel.replace("_", " ")
        else:
            mylabel = str(directory[j])
        
        # now, read orbitN.dat or bigplanetN.dat files
        for k in range(nbplanets):
            
            if fargo3d == 'No':
                if fargo_orig == 'No':
                    time, e, a, M, V, PA, mylambda, varpi = np.loadtxt(directory[j]+"/orbit"+str(k)+".dat",unpack=True)
                else:
                    time, e, a, M, V, PA = np.loadtxt(directory[j]+"/orbit"+str(k)+".dat",unpack=True)
            else:
                if nbcol == 11:
                    time, e, a, M, V, argPA, phiangle, incl, longAN, PA, mylambda = np.loadtxt(directory[j]+"/orbit"+str(k)+".dat",unpack=True)
                else:
                    time, e, a, M, V, argPA, phiangle, incl, longAN, PA = np.loadtxt(directory[j]+"/orbit"+str(k)+".dat",unpack=True)
            
            if (par.mytmin == '#'):
                xmin = time.min() / (2.0*np.pi*rpla_0*np.sqrt(rpla_0))  # time in orbital periods at inner planet's initial location
                umin = 0
            else:
                xmin = par.mytmin
                umin = np.argmin(np.abs(time/(2.0*np.pi*rpla_0*np.sqrt(rpla_0))-xmin))
            if (par.mytmax == '#'):
                xmax = time.max() / (2.0*np.pi*rpla_0*np.sqrt(rpla_0))  # time in orbital periods at inner planet's initial location
                umax = len(time)-1
            else:
                xmax = par.mytmax
                umax = np.argmin(np.abs(time/(2.0*np.pi*rpla_0*np.sqrt(rpla_0))-xmax))
                
            if par.plot_planet[0] == 't':                              
                x = time[umin:umax+1] / (2.0*np.pi*rpla_0*np.sqrt(rpla_0))  # time in orbital periods at inner planet's initial location
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
            
            if par.plot_planet[1] == 'adot':
                y = (a[umin+1:umax+1]-a[umin:umax])/(time[umin+1:umax+1]-time[umin:umax])
                x = time[umin:umax] / (2.0*np.pi*rpla_0*np.sqrt(rpla_0))  # time in orbital periods at inner planet's initial location

            if par.plot_planet[1] == 'r' or par.plot_planet[1] == 'rdot' or par.plot_planet[1] == 'm':
                if fargo3d == 'No':
                    if fargo2d1d == 'Yes':
                        f1, xpla, ypla, f4, f5, mpla, f7, date, f9 = np.loadtxt(directory[j]+"/bigplanet"+str(k)+".dat",unpack=True)
                    else:
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
                    if fargo3d == 'No':
                        y = np.sqrt( xpla*xpla + ypla*ypla )
                    else:
                        y = np.sqrt( xpla*xpla + ypla*ypla + zpla*zpla )
                    y = y[umin:umax+1]
                    if par.physical_units == 'Yes':
                        y *= (culength / 1.5e11) # in au
                if par.plot_planet[1] == 'rdot':
                    if fargo3d == 'No':
                        r = np.sqrt( xpla*xpla + ypla*ypla )
                    else:
                        r = np.sqrt( xpla*xpla + ypla*ypla + zpla*zpla )
                    y = (r[umin+1:umax+1]-r[umin:umax])/(date[umin+1:umax+1]-date[umin:umax])
                    x = mytime[umin:umax]
                if par.plot_planet[1] == 'm':
                    y = mpla[umin:umax+1]
                    if par.physical_units == 'Yes':
                        y *= (cumass / 2e30 / 3e-6) # in Earth masses
                # new (Nov. 2023): display in y-axis log scale (indirect term
                # project)
                if par.log_xyplots_y == 'Yes':
                    y = np.abs(y)
                    ax.set_yscale('log')
                    ytitle = str('|')+ytitle+str('|')
                        
            if par.plot_planet[1] == 'p':
                if k == 0:
                    p0 = (a[umin:umax+1])**(1.5)
                if k == 1:
                    p1 = (a[umin:umax+1])**(1.5)

            if par.plot_planet[1] == 'mmr':
                if k == 0:
                    lambda0 = mylambda[umin:umax+1]
                    varpi0  = PA[umin:umax+1]
                if k == 1:
                    lambda1 = mylambda[umin:umax+1]
                    varpi1  = PA[umin:umax+1]
                    
            if par.plot_planet[1] != 'p' and par.plot_planet[1] != 'mmr':
                if len(directory) == 1:
                    ax.legend_ = None
                else:
                    ax.legend(frameon=False,fontsize=15)
                ax.plot(x[::par.take_one_point_every], y[::par.take_one_point_every], color=par.c20[k*len(directory)+j], lw=2., linestyle = 'solid', label=mylabel)
                fig.add_subplot(ax)
                
        if par.plot_planet[1] == 'p':
            y = p1/p0
            if len(directory) == 1:
                ax.legend_ = None
            ax.plot(x[::par.take_one_point_every], y[::par.take_one_point_every], color=par.c20[j], lw=2., linestyle = 'solid', label=mylabel)
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

    if len(directory) != 1:
        legend = plt.legend(loc='upper right',fontsize=15,facecolor='white',edgecolor='white',framealpha=0.85,numpoints=1,bbox_transform=plt.gcf().transFigure)
        for line, text in zip(legend.get_lines(), legend.get_texts()):
            text.set_color(line.get_color())
        #ax.legend(frameon=False,fontsize=15)
            
  
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
