'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

SPICE data filepaths
'''

import os
import numpy as np
import spiceypy as spice

base_dir = os.path.join(
	os.path.dirname( os.path.realpath( __file__ ) ),
	os.path.join( '..', '..', 'data', 'spice' ) 
	)

# can be changed if you want to load different SPICE kernels.
leapseconds_kernel = os.path.join( base_dir, 'lsk/naif0012.tls' )
de432              = os.path.join( base_dir, 'spk/de432s.bsp'   )
pck00010           = os.path.join( base_dir, 'pck/pck00010.tpc' )



# This function loads the ephemeris from the SPICE kernel in a way consistent with the code from orbit_calculations.py
# This function was written by Jimmy.
def load_ephemeris(target, times, frame, observer):
	if type( target ) == str:
		return np.array( spice.spkezr( target, times, frame, 'NONE', observer )[ 0 ] )
	else:
		return np.array(spice.spkez(target, times, frame, 'NONE', observer)[0])