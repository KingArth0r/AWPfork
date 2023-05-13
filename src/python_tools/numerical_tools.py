'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Numerical Tools Library
'''

# Python standard libraries
import math

# 3rd party libraries
import numpy    as np
import spiceypy as spice


r2d     = 180.0 / np.pi
d2r     = 1.0  / r2d
sec2day = 1.0 / 3600.0 / 24.0
fps2kms = 0.0003048
mi2km   = 1.60934

frame_transform_dict = {
	3: spice.pxform,
	6: spice.sxform
}

def norm( v ):
	return np.linalg.norm( v )

def normed( v ):
	return v / np.linalg.norm( v )
