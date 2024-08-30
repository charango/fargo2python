# =================================================================== 
#                        FARGO to PYTHON
# code written by Clement Baruteau (CB), with contributions from 
# Pablo Benitez-Llambay and Gaylor Wafflard-Fernandez
# =================================================================== 
# 
# program can run with either Python 2.X or Python 3.X.
#
# NB: logarithmic color scale requires matplotblib version > 3.0 to
# display color bar properly
#
# =================================================================== 

# =========================================
#            TO DO LIST
# =========================================
# - think of a better handling of 3D case? Right now I can only plot the primitive fields in a
# vertical plane, not the fields built upon primitive ones (like the vorticity etc.).
# - show RH around planet
# - axifield (minor): add 1D visc evolving profiles, add planet position via dashed line
# =========================================

# Plotting parameters
import matplotlib
matplotlib.rcParams.update({'font.size': 24})
matplotlib.rc('font', family='Arial')

# =====================
# 1. IMPORT GLOBAL VARIABLES
# =====================
import par
if par.allfluids == 'Yes':
    matplotlib.rcParams.update({'font.size': 20})

# =====================
# 2. DISPLAY 2D FIELDS 
# =====================
if par.onedprofile == 'No' and par.plot_field == True:
    from plot_twodfield import *
    plottwodfield()

# =====================
# 3. DISPLAY AZIMUTHALLY-AVERAGED RADIAL 1D PROFILES
# =====================
if par.onedprofile == 'Yes' and par.plot_field == True:
    from plot_onedprofile import *
    plotonedprofile()

# =====================
# 4. DISPLAY TIME EVOLUTION OF DISC TORQUE OR POWER ON PLANET(s)
# =====================
if par.plot_tqwk != 'No':
    from plot_tqwk import *
    plottqwk()

# =====================
# 5. DISPLAY TIME EVOLUTION OF PLANET ORBITAL ELEMENTS
# =====================
if par.plot_planet != 'No':
    from plot_planet import *
    plotplanet()

# =====================
# 6. DISPLAY STUFF FOR DUST PARTICLES (only for Dusty-FARGO ADSG)
# =====================
if par.plot_dust != 'No' and par.fargo3d == 'No':
    from plot_dust import *
    plotdust()

# =====================
# 7. DISC'S CENTER OF MASS
# =====================
if ( ('plot_disccom' in open('paramsf2p.dat').read()) and (par.plot_disccom != 'No') ):
    from plot_disccom import *
    plotdisccom()

# =====================
# 8. DISC MASS (time)
# =====================
if ( ('plot_discmass' in open('paramsf2p.dat').read()) and (par.plot_discmass != 'No') ):
    from plot_discmass import *
    plotdiscmass()

# =====================
# 9. ratio of orbit-crossing and librating inverse vortensities (time)
# =====================
if ( ('plot_libcross' in open('paramsf2p.dat').read()) and (par.plot_libcross != 'No') ):
    from plot_libcross import *
    plotlibcross()