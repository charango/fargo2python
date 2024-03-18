import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import fnmatch
import os
import subprocess
import sys
import re

from mpl_toolkits.axes_grid1 import make_axes_locatable

from mesh import *
from field import *

def plotdust():

    # first import global variables
    import par
    
    # several directories and several output numbers are possible
    directory = par.directory
    if isinstance(par.directory, str) == True:
        directory = [par.directory]
    on = par.on
    if isinstance(par.on, int) == True:
        on = [par.on]
    if par.movie == 'Yes':
        on = range(par.on[0],par.on[1]+1,par.take_one_point_every)
        
    # get time and length units
    myfield = Field(field=par.whatfield, fluid=par.fluid, on=on[0], directory=directory[0], physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile='No', override_units=par.override_units)
    cutime = myfield.cutime
    culength = myfield.culength
        
    # first prepare figure
    if par.plot_dust[0] == 'r':
        xtitle = 'Radius '
        if par.physical_units == 'Yes':
            xtitle += ' [au]'
    if par.plot_dust[0] == 'phi':
        xtitle = 'Azimuth [rad]'
    if par.plot_dust[1] == 'phi':
        ytitle = 'Azimuth [rad]'
    if par.plot_dust[1] == 'size':
        ytitle = 'Dust size [meter]'
    if par.plot_dust[1] == 'Stokes':
        ytitle = 'Stokes number'
    if par.plot_dust[1] == 'vr':
        ytitle = 'Dust radial velocity'
        if par.physical_units == 'Yes':
            ytitle += r' [km s$^{-1}$]'
    if par.plot_dust[1] == 'vphi':
        ytitle = 'Dust azimuthal velocity'
        if par.physical_units == 'Yes':
            ytitle += r' [km s$^{-1}$]'     

    ymin = 1e8
    ymax = -1e8
    xmin = 1e8
    xmax = -1e8

    if par.movie == 'No':
        # first prepare figure
        fig = plt.figure(figsize=(8.,8.))
        plt.subplots_adjust(left=0.16, right=0.96, top=0.95, bottom=0.12)
        ax = fig.gca()
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)
        ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
        
    # now plot: one figure for each output number, several directories can be dealt with
    for k in range(len(on)):   # loop over output numbers
        if par.movie == 'Yes':
            print('animation: output number '+str(k)+' / '+str(len(on)-1),end='\r')
            fig = plt.figure(figsize=(8.,8.))
            plt.subplots_adjust(left=0.16, right=0.96, top=0.95, bottom=0.12)
            ax = fig.gca()
            ax.set_xlabel(xtitle)
            ax.set_ylabel(ytitle)
            ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')

        # get time
        myfield  = Field(field=par.whatfield, fluid=par.fluid, on=on[k], directory=directory[0], physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile=par.onedprofile, override_units=par.override_units)
        mylabel = myfield.strtime
            
        for j in range(len(directory)):   # loop over directories

            # read dustsysatX.dat files
            (rd, td, vrd, vtd, Stokes, sizedust) = np.loadtxt(directory[j]+'/dustsystat'+str(on[k])+'.dat', unpack=True)
            if par.physical_units == 'Yes':
                rd *= (culength / 1.5e11)    # in au
                vrd *= 1e-3*culength/cutime  # in km / s
                vtd *= 1e-3*culength/cutime  # in km / s
                # sizedust is already in meters
            if par.plot_dust[0] == 'r':
                x = rd
                if (par.myrmin == '#'):
                    xmin = x.min()
                else:
                    xmin = par.myrmin
                if (par.myrmax == '#'):
                    xmax = x.max()
                else:
                    xmax = par.myrmax
            if par.plot_dust[0] == 'phi':
                x = td
                if (par.myphimin == '#'):
                    xmin = x.min()
                else:
                    xmin = par.myphimin
                if (par.myphimax == '#'):
                    xmax = x.max()
                else:
                    xmax = par.myphimax
            if par.plot_dust[1] == 'phi':
                y = td
                if (par.myphimin == '#'):
                    ymin = y.min()
                else:
                    ymin = par.myphimin
                if (par.myphimax == '#'):
                    ymax = y.max()
                else:
                    ymax = par.myphimax
            if par.plot_dust[1] == 'size':
                y = sizedust
                if (par.sizemin == '#'):
                    ymin = y.min()
                else:
                    ymin = par.sizemin
                if (par.sizemax == '#'):
                    ymax = y.max()
                else:
                    ymax = par.sizemax
                ax.set_yscale('log')
            if par.plot_dust[1] == 'Stokes':
                y = Stokes
                ax.set_yscale('log')
            if par.plot_dust[1] == 'vr':
                y = vrd
            if par.plot_dust[1] == 'vphi':
                y = vtd
            if par.plot_dust[1] == 'vphi' or par.plot_dust[1] == 'vr' or par.plot_dust[1] == 'Stokes':           
                if (par.fieldmin == '#'):
                    if y.min() < ymin:
                        ymin = y.min()
                else:
                    ymin = par.fieldmin
                if (par.fieldmax == '#'):
                    if y.max() > ymax:
                        ymax = y.max()
                else:
                    ymax = par.fieldmax
            #
            ax.set_xlim(xmin,xmax)
            ax.set_ylim(ymin,ymax)
            if par.movie == 'Yes':
                mycolor = par.c20[j]
            else:
                mycolor = par.c20[k*len(directory)+j]
            ax.scatter(x, y, s=1, c=mycolor, alpha=0.3)
            ax.legend(frameon=False,fontsize=15)
            fig.add_subplot(ax)

            # show time in plot's top-right corner
            xstr = 0.99*xmax
            ystr = 0.97*ymax
            colorstr = 'black'
            ax.text(xstr,ystr,mylabel,fontsize=18,color=colorstr,weight='bold',horizontalalignment='right',verticalalignment='top')

        # save file
        if len(directory) == 1:           
            outfile = par.plot_dust[0]+'_'+par.plot_dust[1]+'_'+str(directory[0])+'_'+str(on[k]).zfill(4)
        else:
            outfile = par.plot_dust[0]+'_'+par.plot_dust[1]+'_'+str(on[k]).zfill(4)
        fileout = outfile+'.pdf'
        if par.saveaspdf == 'Yes':
            plt.savefig('./'+fileout, dpi=160)
        if par.saveaspng == 'Yes':
            plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)

    # finally concatenate png if movie requested
    if par.movie == 'Yes':
        if len(directory) == 1:
            # png files that have been created above
            allpngfiles = [par.plot_dust[0]+'_'+par.plot_dust[1]+'_'+str(directory[0])+'_'+str(on[x]).zfill(4)+'.png' for x in range(len(on))]
            # input files for ffpmeg
            input_files = par.plot_dust[0]+'_'+par.plot_dust[1]+'_'+str(directory[0])+'_%04d.png'
            # output file for ffmpeg
            filempg = par.plot_dust[0]+'_'+par.plot_dust[1]+'_'+str(directory[0])+'_'+str(on[0])+'_'+str(on[len(on)-1])+'.mpg'
        else:
            # png files that have been created above
            allpngfiles = [par.plot_dust[0]+'_'+par.plot_dust[1]+'_'+str(on[x]).zfill(4)+'.png' for x in range(len(on))]
            # input files for ffpmeg
            input_files = par.plot_dust[0]+'_'+par.plot_dust[1]+'_%04d.png'
            # output file for ffmpeg
            filempg = par.plot_dust[0]+'_'+par.plot_dust[1]+'_'+str(on[0])+'_'+str(on[len(on)-1])+'.mpg'
        # call to python-ffmpeg
        import ffmpeg
        (
            ffmpeg            
            .input(input_files, framerate=10, start_number=str(on[0]))
            # framerate=10 means the video will play at 10 of the original images per second
            .output(filempg, r=30, pix_fmt='yuv420p', **{'qscale:v': 3})
            # r=30 means the video will play at 30 frames per second
            .overwrite_output()
            .run()
        )
        # erase png files        
        allfiles = ' '.join(allpngfiles)
        os.system('rm -f '+allfiles)
