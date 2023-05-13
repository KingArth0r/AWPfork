
import numpy as np
import spiceypy as spice

# This function loads the ephemeris from the spice kernel in a way consistent with the code from orbit_calculations.py
# This function was written by Jimmy.
def load_ephemeris(target, times, frame, observer):
	if type( target ) == str:
		return np.array( spice.spkezr( target, times, frame, 'NONE', observer )[ 0 ] )
	else:
		return np.array(spice.spkez(target, times, frame, 'NONE', observer)[0])