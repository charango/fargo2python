import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import fnmatch
import os
import subprocess
import sys
import re

import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, LogFormatter)
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mesh import *
from field import *


def plotspacetimediagram():
    
    # first import global variables
    import par
    
    # several output numbers and several directories
    on = par.on
    on = range(on[0],on[1]+1,par.take_one_point_every)
    directory = par.directory
    if isinstance(par.directory, str) == True:
        directory = [par.directory]       
 

    # =====================
    # loop over directories
    # =====================
    for j in range(len(directory)):  
        
        # read field at time=0
        myfield0  = Field(field=par.whatfield, fluid=par.fluid, on=0, directory=directory[j], physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile=par.onedprofile, z_average=par.z_average, override_units=par.override_units)

        strfield = myfield0.strname

        # In case myfieldmin or myfieldmax aren't set: we simply take the min or max of the field at time = 0
        if par.fieldmin == '#':
            ymin = myfield0.data.min()
        else:
            ymin = par.fieldmin
        if par.fieldmax == '#':
            ymax = myfield0.data.max()
        else:
            ymax = par.fieldmax

        # get r, xmin, xmax for displaying puroposes
        R = myfield0.redge
        if par.physical_units == 'Yes':
            R *= (myfield0.culength / 1.5e11) # in au
        if (par.myrmin == '#'):
            xmin = R.min()
        else:
            xmin = par.myrmin
        if (par.myrmax == '#'):
            xmax = R.max()
        else:
            xmax = par.myrmax

        # define and allocate array for space-time diagram
        spacetime_array = np.zeros((len(myfield0.rmed),len(on)-1))

        # ------------------------
        # loop over output numbers
        # ------------------------
        for k in range(len(on)):     

            print('directory number '+str(j)+' / '+str(len(directory)-1),' and output number '+str(k)+' / '+str(len(on)-1),end='\r')

            # read field
            myfield  = Field(field=par.whatfield, fluid=par.fluid, on=on[k], directory=directory[j], physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile=par.onedprofile, z_average=par.z_average, override_units=par.override_units)

            if par.nodiff == 'No':
                array = (myfield.data-myfield0.data)/myfield0.data
            else:
                array = myfield.data
                # conversion in physical units
                if par.physical_units == 'Yes':
                    array = myfield.data * myfield.unit
                if par.log_xyplots_y == 'Yes' and (par.whatfield == 'vrad' or par.whatfield == 'vy'):
                    #print('1D vrad displayed with log y-scale')
                    array = np.abs(array)
                    
            axiarray = np.sum(array,axis=1)/myfield.nsec

            if par.onedprofile == 'Cut':
                axiarray = array[:,0]    # azimuthal cut at zero azimuth (j=0)

            if par.onedprofile == 'Median':
                axiarray = np.median(array,axis=1)   # median over azimuth of the density profile

            # save into spacetime_array array
            if k != len(on)-1:
                spacetime_array[:,k] = axiarray

        # ------------
        # Final figure (for each directory)
        # ------------
        fig = plt.figure(figsize=(8.,8.))
        plt.subplots_adjust(left=0.16, right=0.95, top=0.88, bottom=0.12)
        ax = fig.gca()

        if par.physical_units == 'Yes':
            xtitle = 'radius [au]'
        else:
            xtitle = 'radius [code units]'
        ax.set_xlabel(xtitle)
        if par.log_xyplots_x == 'Yes':
            ax.set_xscale('log')

        ax.set_ylabel(r'Time [$T_0$]')

        mycolormap=par.mycolormap
        if par.log_xyplots_y == 'Yes':
            mynorm = matplotlib.colors.LogNorm(vmin=ymin,vmax=ymax)
        else:
            mynorm = matplotlib.colors.Normalize(vmin=ymin,vmax=ymax)
        
        CF = ax.pcolormesh(R,np.arange(len(on)),np.transpose(spacetime_array),cmap=mycolormap,norm=mynorm,rasterized=True)

        # plot color-bars
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="2.5%", pad=0.12)
        cb = plt.colorbar(CF, cax=cax, orientation='horizontal')
        cax.xaxis.tick_top()
        cax.xaxis.set_tick_params(direction='out')
        #cax.xaxis.set_major_locator(plt.MaxNLocator(3))
        cax.xaxis.set_major_locator(plt.MaxNLocator(4))
        # title on top
        cax.xaxis.set_label_position('top')
        cax.set_xlabel(strfield)
        cax.xaxis.labelpad = 8
        if par.log_colorscale == 'Yes' and par.fieldmin != 'auto' and par.fieldmax != 'auto':
            cax.xaxis.set_major_locator(ticker.LogLocator(base=10.0,numticks=8))
        if par.log_colorscale == 'Yes' and (par.fieldmin == 'auto' or par.fieldmax == 'auto'):
            cax.xaxis.set_tick_params(direction='out')
        
        # save figure
        prefix = 'SpaceTime_axi'
        if par.onedprofile == 'Cut':
            prefix = 'SpaceTime_cut'
        if par.onedprofile == 'Median':
            prefix = 'SpaceTime_median'
        outfile = prefix+par.fluid+'_'+par.whatfield+'_'+str(directory[j])
        fileout = outfile+'.pdf'
        if par.saveaspdf == 'Yes':
            plt.savefig('./'+fileout, dpi=160)
        if par.saveaspng == 'Yes':
            plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)
        plt.close(fig)  # close figure as we reopen figure at every output number

        # save 2D array in binary file
        fileout = outfile+'.dat'
        FILEOUT = open(fileout,'wb')        # binary format
        np.transpose(spacetime_array).tofile(FILEOUT)
        FILEOUT.close()