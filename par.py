import matplotlib
import matplotlib.pyplot as plt
import sys
import os
import re
import subprocess

# special color scale for 1D plots:
c20 = [(31, 119, 180), (255, 127, 14), (174, 199, 232), (255, 187, 120),    
       (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
       (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
       (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
       (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
for i in range(len(c20)):
    r, g, b = c20[i]    
    c20[i] = (r / 255., g / 255., b / 255.)


params = open("paramsf2p.dat",'r')
lines_params = params.readlines()
params.close()                 

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
                if name == 'directory' or name == 'plot_planet' or name == 'plot_dust':
                    value = [str(x) for x in value.split(',')]
            else:
                value = '"' + value + '"'   # if none of the above tests work, we know value is a string
    par.append(value)
    var.append(name)

for i in range(len(par)):
    exec(var[i]+"="+str(par[i]))

# Below we work out specific parameters and/or error messages
# was simulation carried out with Fargo3D?
fargo3d = 'No'
if isinstance(directory, str) == False:
    input_file = directory[0]+'/summary0.dat'
else:
    input_file = directory+'/summary0.dat'
if os.path.isfile(input_file) == True:
    fargo3d = 'Yes'    
    
# global boolean: if True, then plot 1D or 2D fields
plot_field = True
if plot_tqwk != 'No' or plot_planet != 'No':
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
if whatfield == 'stokes':
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

# case where all fluids should be displayed!
if fluid == 'all':
    import numpy as np
    allfluids = 'Yes'
    if isinstance(directory, str) == False:
        input_file = directory[0]+'/dustsizes.dat'
    else:
        input_file = directory+'/dustsizes.dat'
    dust_id, dust_size, dust_gas_ratio = np.loadtxt(input_file,unpack=True)
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
