import matplotlib
import matplotlib.pyplot as plt
import sys
import os
import re
import subprocess

from shutil import which

# special color scale for 1D plots:
c20 = [(31, 119, 180), (255, 128, 128), (0, 153, 76), (255, 187, 120),
       (214, 39, 40), (174, 199, 232), (152, 223, 138), (255, 152, 150),    
       (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
       (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
       (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
'''
c20 = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (255, 187, 120),
       (214, 39, 40), (174, 199, 232), (152, 223, 138), (255, 152, 150),    
       (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
       (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
       (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
'''
'''
c20 = [(31, 119, 180), (255, 127, 14), (255, 187, 120), (174, 199, 232),     
       (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
       (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
       (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
       (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
'''   
for i in range(len(c20)):
    r, g, b = c20[i]    
    c20[i] = (r / 255., g / 255., b / 255.)


params = open("paramsf2p.dat",'r')
lines_params = params.readlines()
params.close()                 

# Global booleans set to No by default
plot_disccom = 'No'
plot_discmass = 'No'
plot_libcross = 'No'
plot_turb = 'No'

par = []                       # allocating a dictionary
var = []                       # allocating a dictionary
regex = re.compile(',')

for line in lines_params:      
    try:
        line.split()[0][0]     # check if blank line (GWF)
    except:
        continue
    if (line.split()[0][0]=='#'): # check if line starts with a # (comment)
        continue
    else:
        name, value = line.split()[0:2]
    try:
        float(value)           # first trying with float
    except ValueError:         # if it is not float
        try:
            int(value)         # we try with integer
        except ValueError:     # if it is not integer nor float
            if(regex.search(value) != None):  # we try with array with several values separated by a comma (,)
                if name == 'on' or name == 'mmr_integers':
                    value = [int(x) for x in value.split(',')]
                if name == 'directory' or name == 'plot_planet' or name == 'plot_dust' or name == 'use_legend' or name == 'filename':
                    value = [str(x) for x in value.split(',')]
            else:
                value = '"' + value + '"'   # if none of the above tests work, we know value is a string
    par.append(value)
    var.append(name)

for i in range(len(par)):
    exec(var[i]+"="+str(par[i]))

# Check whether awk or gawk are installed!
if which('gawk') is not None:
    awk_command = 'gawk'
else:
    if which('awk') is not None:
        awk_command = 'awk'
    else:
        print('neither gawk not awk are installed on your system! I cannot use them to extract relevant parameters from your .par parameter file in the simulation directory. Please install either awk or gawk.')
    
if not('verbose' in open('paramsf2p.dat').read()):
    verbose = 'No'

# Below we work out specific parameters and/or error messages
# was simulation carried out with Fargo3D?
fargo3d = 'No'
fargo_orig = 'No'
if isinstance(directory, str) == False:
    summary0_file = directory[0]+'/summary0.dat'
    usedazi_file  = directory[0]+'/used_azi.dat'
else:
    summary0_file = directory+'/summary0.dat'
    usedazi_file  = directory+'/used_azi.dat'
if os.path.isfile(summary0_file) == True:
    # Simulations were carried out with Fargo3D
    fargo3d = 'Yes'
else:
    # Simulations were carried out with Fargo2D
    fargo3d = 'No'
    if os.path.isfile(usedazi_file) == True:
    # Simulations were carried out with Dusty FARGO-ADSG
        fargo_orig = 'No'
    else:
    # Simulations were carried out with the original FARGO code
        fargo_orig = 'Yes'

# global boolean: if True, then plot 1D or 2D fields
plot_field = True

if ( (plot_tqwk != 'No') or (plot_planet != 'No') or (plot_disccom != 'No') or (plot_discmass != 'No') or (plot_libcross != 'No') or (plot_turb != 'No') ):
    plot_field = False
    movie = 'No'
    if plot_planet[1] == 'mmr':
        import fnmatch
        nbplanets = len(fnmatch.filter(os.listdir(directory), 'orbit*.dat'))
        if nbplanets <= 1:
            sys.exit('ERROR: mean-motion resonant angles requested, but only one planet was included in the simulation!')
        else:
            print('integers used to compute mean-motion resonant angles : ', mmr_integers)
if plot_dust != 'No':
    plot_field = False
            
# case an animation is requested
if movie == 'Yes':
    saveaspdf = 'No'
    saveaspng = 'Yes'
    if isinstance(on, int) == True:
        sys.exit("ERROR: you requested an animation but specified a single output number: on needs to be specified as min,max -- eg: 1,10 for output numbers from 1 to 10")

# units
if whatfield == 'stokes' or whatfield == 'betacooling':
    physical_units = 'Yes'
if override_units == 'Yes':
    if new_unit_length == 0.0:
        sys.exit('override_units set to yes but new_unit_length is not defined in params.dat, I must exit!')
    if new_unit_mass == 0.0:
        sys.exit('override_units set to yes but new_unit_mass is not defined in params.dat, I must exit!')

if (take_one_point_every == '#'):
    take_one_point_every = 1

if (nodiff == 'No'):
    log_colorscale = 'No'

# case magnetic field is displayed in FARGO3D runs: convert fluid to '' as 
# fieldnames are simply, eg, bx123.dat and not gasbx123.dat
if whatfield == 'bx' or whatfield == 'by' or whatfield == 'bz':
    fluid = ''

# case polargrid Test has been used in simulation
if whatfield == 'Test':
    fluid = ''

# case direct or indirect torque of disc on planet is computed
if whatfield == 'direct_torque' or whatfield == 'indirect_torque':
    nodiff = 'Yes'

# case where all fluids should be displayed!
if fluid == 'all':
    import numpy as np
    allfluids = 'Yes'
    if isinstance(directory, str) == False:
        input_file = directory[0]+'/dustsizes.dat'
        if os.path.isfile(input_file):
            dust_id, dust_size, dust_gas_ratio = np.loadtxt(input_file,unpack=True)
        else:
            input_file = directory[0]+'/duststokesnb.dat'
            if os.path.isfile(input_file):
                dust_id, dust_stokes = np.loadtxt(input_file,unpack=True)
                dust_size = np.zeros(len(dust_id))
    else:
        input_file = directory+'/dustsizes.dat'
        if os.path.isfile(input_file):
            dust_id, dust_size, dust_gas_ratio = np.loadtxt(input_file,unpack=True)
        else:
            input_file = directory+'/duststokesnb.dat'
            if os.path.isfile(input_file):
                dust_id, dust_stokes = np.loadtxt(input_file,unpack=True)
                dust_size = np.zeros(len(dust_id))
        
    nbfluids = 1+len(dust_id)  # +1 = gas
    # redefine fluid
    fluids = [str(x) for x in range(nbfluids)]
    fluids[0] = 'gas'
    for i in range(1,nbfluids):
        fluids[i] = 'dust'+str(i)
    fluid = 'allfluids'
else:
    allfluids = 'No'
    fluids = [fluid]
    
# Color map
if not('mycolormap' in open('paramsf2p.dat').read()):
    mycolormap = 'nipy_spectral'
