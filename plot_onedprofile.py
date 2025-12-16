import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import fnmatch
import os
import subprocess
import sys
import re

from mesh import *
from field import *


def plotonedprofile():
    
    # first import global variables
    import par
    
    # several output numbers and several directories
    on = par.on
    if isinstance(par.on, int) == True:
        on = [on]
    if par.movie == 'Yes':
        on = range(on[0],on[1]+1,par.take_one_point_every)
    directory = par.directory
    if isinstance(par.directory, str) == True:
        directory = [par.directory]       
 
    if par.physical_units == 'Yes':
        xtitle = 'radius [au]'
    else:
        xtitle = 'radius [code units]'

    # if fieldmin and fieldmax are undefined, find out min and max of
    # y-values to be displayed on plot by going though all directories
    # and all output numbers:
    if (par.fieldmin == '#') and (par.fieldmax == '#'):
        if par.verbose == 'Yes':
            print('fieldmin and fieldmax are unspecified, I will set them automatically...')
        ymin = 1e8
        ymax = -1e8
        
        for j in range(len(directory)):  # loop over directories
            if par.nodiff == 'No':
                myfield0  = Field(field=par.whatfield, fluid=par.fluid, on=0, directory=directory[j], physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile=par.onedprofile, z_average=par.z_average, override_units=par.override_units)

            for k in range(len(on)):     # loop over output numbers
                myfield  = Field(field=par.whatfield, fluid=par.fluid, on=on[k], directory=directory[j], physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile=par.onedprofile, z_average=par.z_average, override_units=par.override_units)

                if j==0 and k==0:
                    R = myfield.rmed
                    if par.physical_units == 'Yes':
                        R *= (myfield.culength / 1.5e11) # in au
                    myrmin = par.myrmin
                    if (par.myrmin == '#'):
                        myrmin = R.min()
                    imin = np.argmin(np.abs(R-myrmin))
                    myrmax = par.myrmax
                    if (par.myrmax == '#'):
                        myrmax = R.max()
                    imax = np.argmin(np.abs(R-myrmax))   

                if par.nodiff == 'No':
                    array = (myfield.data-myfield0.data)/myfield0.data
                else:
                    array = myfield.data
                    # conversion in physical units
                    if par.physical_units == 'Yes':
                        array = myfield.data * myfield.unit
                    if par.log_xyplots_y == 'Yes' and (par.whatfield == 'vrad' or par.whatfield == 'vy'):
                        array = np.abs(array)
                        
                axiarray = np.sum(array[imin:imax,:],axis=1)/myfield.nsec
                if par.onedprofile == 'Cut':
                    axiarray = array[:,0]    # azimuthal cut at zero azimuth (j=0)
                
                if axiarray.min() < ymin:
                    ymin = axiarray.min()
                if axiarray.max() > ymax:
                    ymax = axiarray.max()

        if par.verbose == 'Yes':              
            print('fieldmin = ', ymin)
            print('fieldmin = ', ymax)
                
    else:
        if (par.fieldmin != '#'):
            ymin = par.fieldmin
        else:
            ymin = 0.0
        if (par.fieldmax != '#'):    
            ymax = par.fieldmax
        else:
            ymax = 2.0*ymin  # CUIDADIN
                    
    if par.movie == 'No':
        # first prepare figure
        fig = plt.figure(figsize=(8.,8.))
        plt.subplots_adjust(left=0.20, right=0.96, top=0.94, bottom=0.12)
        ax = fig.gca()
        ax.set_xlabel(xtitle)
        if par.log_xyplots_y == 'Yes':
            ax.set_yscale('log')
        if par.log_xyplots_x == 'Yes':
            ax.set_xscale('log')

    # -----------------------------
    # then plot via double for loop
    # -----------------------------
    for k in range(len(on)):     # loop over output numbers
        
        if par.movie == 'Yes':
            print('animation: output number '+str(k)+' / '+str(len(on)-1),end='\r')
            # first prepare figure
            fig = plt.figure(figsize=(8.,8.))
            plt.subplots_adjust(left=0.20, right=0.96, top=0.94, bottom=0.12)
            ax = fig.gca()
            ax.set_xlabel(xtitle)
            if par.log_xyplots_y == 'Yes':
                ax.set_yscale('log')
            if par.log_xyplots_x == 'Yes':
                ax.set_xscale('log')
                
        for j in range(len(directory)):  # loop over directories
            # it does not take much time to read all fields again...
            myfield  = Field(field=par.whatfield, fluid=par.fluid, on=on[k], directory=directory[j], physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile=par.onedprofile, z_average=par.z_average, override_units=par.override_units)
            
            # stuff we do only once: set xmin and xmax
            if k == 0 and j == 0:
                R = myfield.redge
                if par.physical_units == 'Yes':
                    R *= (myfield.culength / 1.5e11) # in au
                if (par.myrmin == '#'):
                    xmin = R.min()
                else:
                    xmin = par.myrmin
                if (par.myrmax == '#'):
                    xmax = R.max()
                else:
                    xmax = par.myrmax
                strfield = myfield.strname
            
            if par.nodiff == 'No':
                myfield0  = Field(field=par.whatfield, fluid=par.fluid, on=0, directory=directory[j], physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, onedprofile=par.onedprofile, z_average=par.z_average, override_units=par.override_units)
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

            #axiarray = axiarray[2:-2]  # CUIDADIN!!
            R = myfield.rmed
            #R = R[2:-2] # CUIDADIN!!
            
            if par.physical_units == 'Yes':
                R *= (myfield.culength / 1.5e11) # in au
            mylabel = myfield.strtime

            if len(directory) > 1:
                if ('use_legend' in open('paramsf2p.dat').read()) and (par.use_legend != '#'):
                    if len(directory) == 1:
                        mylabel = str(par.use_legend) + ', '+ mylabel
                    else:
                        mylabel = str(par.use_legend[j]) + ', '+ mylabel
                else:
                    mylabel = str(directory[j]) + ', '+ mylabel
            if par.movie == 'Yes':
                mycolor = par.c20[j]
            else:
                mycolor = par.c20[k*len(directory)+j]

            ax.plot(R, axiarray, color=mycolor, lw=2., linestyle = 'solid', label=mylabel)
            ax.set_ylabel(strfield)

            if ( ('dynamical_colorscale' in open('paramsf2p.dat').read()) and (par.dynamical_colorscale == 'Yes') ):
                ymin = axiarray.min()
                ymax = axiarray.max()

            ax.set_ylim(ymin,ymax)
            ax.set_xlim(xmin,xmax)
            ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
            #plt.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            if par.movie == 'Yes':
                ax.legend(frameon=False,fontsize=15,loc='upper left')
            else:
                ax.legend(frameon=False,fontsize=15)
            fig.add_subplot(ax)

            # option to write result in 1D ascii file
            if ( ('write_ascii' in open('paramsf2p.dat').read()) and (par.write_ascii == 'Yes') ):
                 ascii = open('1D'+directory[j]+par.whatfield+str(on[k])+'.dat','w')
                 for v in range(len(R)):
                    ascii.write(str(R[v])+'\t'+str(axiarray[v])+'\n')
                    

        # save file
        if len(directory) == 1:           
            outfile = 'axi'+par.fluid+'_'+par.whatfield+'_'+str(directory[0])+'_'+str(on[k]).zfill(4)
            if par.movie == 'Yes' and par.take_one_point_every != 1:
                outfile = 'axi'+par.fluid+'_'+par.whatfield+'_'+str(directory[0])+'_'+str(k).zfill(4)
        else:
            outfile = 'axi'+par.fluid+'_'+par.whatfield+'_'+str(on[k]).zfill(4)
            if par.movie == 'Yes' and par.take_one_point_every != 1:
                outfile = 'axi'+par.fluid+'_'+par.whatfield+'_'+str(k).zfill(4)
        fileout = outfile+'.pdf'
        if par.saveaspdf == 'Yes':
            plt.savefig('./'+fileout, dpi=160)
        if par.saveaspng == 'Yes':
            plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=120)
        if par.movie == 'Yes':
            plt.close(fig)  # close figure as we reopen figure at every output number

            
    # finally concatenate png if movie requested
    if par.movie == 'Yes':
        if len(directory) == 1:
            # png files that have been created above
            allpngfiles = ['axi'+par.fluid+'_'+par.whatfield+'_'+str(directory[0])+'_'+str(on[x]).zfill(4)+'.png' for x in range(len(on))]
            if par.take_one_point_every != 1:
                allpngfiles = ['axi'+par.fluid+'_'+par.whatfield+'_'+str(directory[0])+'_'+str(x).zfill(4)+'.png' for x in range(len(on))]
            # input files for ffpmeg
            input_files = 'axi'+par.fluid+'_'+par.whatfield+'_'+str(directory[0])+'_%04d.png'
            # output file for ffmpeg
            filempg = 'axi'+par.fluid+'_'+par.whatfield+'_'+str(directory[0])+'_'+str(on[0])+'_'+str(on[len(on)-1])+'.mpg'    
        else:
            # png files that have been created above
            allpngfiles = ['axi'+par.fluid+'_'+par.whatfield+'_'+str(on[x]).zfill(4)+'.png' for x in range(len(on))]
            if par.take_one_point_every != 1:
                allpngfiles = ['axi'+par.fluid+'_'+par.whatfield+'_'+str(x).zfill(4)+'.png' for x in range(len(on))]
            # input files for ffpmeg
            input_files = 'axi'+par.fluid+'_'+par.whatfield+'_%04d.png'
            # output file for ffmpeg
            filempg = 'axi'+par.fluid+'_'+par.whatfield+'_'+str(on[0])+'_'+str(on[len(on)-1])+'.mpg'
        # options
        if par.take_one_point_every != 1:
            str_on_start_number = str(0)
        else:
            str_on_start_number = str(on[0])
        if par.nodiff == 'Yes':
            filempg = re.sub('.mpg', '_nodiff.mpg', filempg)
        if par.z_average == 'Yes':
            filempg = re.sub('.mpg', '_zave.mpg', filempg)
        # call to python-ffmpeg
        import ffmpeg
        (
            ffmpeg            
            .input(input_files, framerate=10, start_number=str_on_start_number)
            # framerate=10 means the video will play at 10 of the original images per second
            .output(filempg, r=30, pix_fmt='yuv420p', **{'qscale:v': 3})
            # r=30 means the video will play at 30 frames per second
            .overwrite_output()
            .run()
        )
        # erase png files
        allfiles = ' '.join(allpngfiles)
        os.system('rm -f '+allfiles)
