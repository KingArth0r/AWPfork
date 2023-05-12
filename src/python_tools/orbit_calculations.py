'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Orbit Calculations Library
'''

# 3rd party libraries
import numpy as np
import math
import spiceypy as spice
import matplotlib.pyplot as plt

# AWP library
import numerical_tools as nt
import lamberts_tools as lt
import planetary_data as pd
import spice_data
import spice_tools as st
import spice_data as sd

def interplanetary_porkchop(config):
    _config = {
        'planet0': 'Earth',
        'planet1': 'MARS BARYCENTER',
        'departure0': '2020-07-01',
        'departure1': '2020-09-01',
        'arrival0': '2020-11-01',
        'arrival1': '2022-01-24',
        'mu': pd.sun['mu'],
        'step': 86400, # Seconds to days (default step is 1 day)
        'frame': 'ECLIPJ2000',
        'observer': 'SOLAR SYSTEM BARYCENTER',
        'cutoff_v': 20.0,
        'c3_levels': None,
        'vinf_levels': None,
        'tof_levels': None,
        'dv_levels': None,
        'dv_cmap': 'RdPu_r',
        'figsize': (10,20),
        'lw': 1.5,
        'title': 'Porkchop Plot',
        'fontsize': 15,
        'show': True,
        'filename': 'C:/Users/space/OneDrive/Documents/UW/514/Final Project Code/porkchop.png',
        'filename_dv': None,
        'dpi': 300,
        'load': False
    }

    spice.furnsh(sd.leapseconds_kernel) # keeps track of leapseconds
    spice.furnsh(sd.pck00010)
    spice.furnsh(sd.de432)

    for key in config.keys():
        _config[key] = config[key]
    cutoff_c3 = _config['cutoff_v']**2 # cutoff_c3 should be escape velocity ** 2

    # arrays for arrivals and departures
    et_departures = np.arange(
        spice.utc2et(_config['departure0']),
        spice.utc2et(_config['departure1']) + _config['step'],
        _config['step'])
    et_arrivals = np.arange(
        spice.utc2et(_config['arrival0']),
        spice.utc2et(_config['arrival1']) + _config['step'],
        _config['step'])
    # calculate the number of days
    ds = len(et_departures) # departure days
    as_ = len(et_arrivals) # arrival days
    total = ds*as_ # total days

    print(f'Departure days: {ds}.')
    print(f'Arrival days: {as_}.')
    print(f'Total combinations: {total}')

    #create arrays for C3, v_inf, and tof
    C3_shorts = np.zeros((as_, ds))
    C3_longs = np.zeros((as_, ds))
    v_inf_shorts = np.zeros((as_, ds))
    v_inf_longs = np.zeros((as_, ds))
    tofs = np.zeros((as_, ds))

    # create coords
    X = np.arange(ds)
    Y = np.arange(as_)

    # loop through every combination
    for na in Y:
        for nd in X:
            # state of planet0 (Earth) at departure
            state_depart = st.load_ephemeris(
                _config['planet0'],
                [et_departures[nd]],
                _config['frame'], _config['observer'])[0]

            # state of planet1 (Mars) at arrival
            state_arrive = st.load_ephemeris(
                _config['planet1'],
                [et_arrivals[na]],
                _config['frame'], _config['observer'])[0]

            # calculate time of flight
            tof = et_arrivals[na] - et_departures[nd]

            tm = 1
            mu = _config['mu']

            # try lambert's solver. If lambert's solver fails, exception is thrown
            try:
                v_sc_depart_short, v_sc_arrive_short = lt.lamberts_universal_variables(
                    state_depart[ :3 ], state_arrive[ :3 ],
                    tof, {'tm' : 1, 'mu' : _config[ 'mu' ]})
            except Exception as e:
                print(e)
                v_sc_depart_short = np.array([1000, 1000, 1000])
                v_sc_arrive_short = np.array([1000, 1000, 1000])
            try:
                v_sc_depart_long, v_sc_arrive_long = lt.lamberts_universal_variables(
                    state_depart[ :3 ],state_arrive[ :3 ],
                    tof, {'tm' : 1, 'mu' : _config[ 'mu' ]})
            except Exception as e:
                print(e)
                v_sc_depart_long = np.array([1000, 1000, 1000])
                v_sc_arrive_long = np.array([1000, 1000, 1000])

            C3_short = nt.norm(v_sc_depart_short - state_depart[3:])**2
            C3_long = nt.norm(v_sc_depart_long - state_depart[3:])**2
            if C3_short > cutoff_c3: C3_short = cutoff_c3
            if C3_long > cutoff_c3: C3_long = cutoff_c3

            v_inf_short = nt.norm(v_sc_arrive_short - state_arrive[3:])
            v_inf_long = nt.norm(v_sc_arrive_long - state_arrive[3:])

            if v_inf_short > _config['cutoff_v']: v_inf_short = _config['cutoff_v']
            if v_inf_long > _config['cutoff_v']: v_inf_long = _config['cutoff_v']

            C3_shorts[na, nd] = C3_short
            C3_longs[na, nd] = C3_long
            v_inf_shorts[na, nd] = v_inf_short
            v_inf_longs[na,nd] = v_inf_long
            tofs[na, nd] = tof

        print("Made it") #print statement for debugging

    tofs /= (3600.0*24.0) #convert from seconds to days

    # total delta V
    dv_shorts = v_inf_shorts + np.sqrt(C3_shorts)
    dv_longs = v_inf_longs + np.sqrt(C3_longs)

    # levels - creating margins of porkchop plot, feel free to change the numbers
    if _config['c3_levels'] is None:
        _config['c3_levels'] = np.arange(10, 50, 2)
    if _config['vinf_levels'] is None:
        _config['vinf_levels'] = np.arange(0, 15, 1)
    if _config['tof_levels'] is None:
        _config['tof_levels'] = np.arange(100, 500, 20)
    if _config['dv_levels'] is None:
        _config['dv_levels'] = np.arange(3, 20, 0.5)

    lw = _config['lw']


    c3levels = _config['c3_levels']
    vinflevels = _config['vinf_levels']
    toflevels = _config['tof_levels']
    color1 = 'm'
    color2 = 'deepskyblue'
    color3 = 'white'


    fig, ax = plt.subplots(figsize = _config['figsize'])

    #contours
    c0 = ax.contour(C3_shorts, levels = c3levels, colors = color1, linewidths = lw)
    c1 = ax.contour(C3_longs, levels = c3levels, colors = color1, linewidths = lw)
    c2 = ax.contour(v_inf_shorts, levels = vinflevels, colors =color2, linewidths = lw)
    c3 = ax.contour(v_inf_longs, levels = vinflevels, colors = color2, linewidths = lw)
    c4 = ax.contour(tofs, levels = toflevels, colors = color3, linewidths = lw*0.6)

    plt.clabel(c0, fmt = '%i')
    plt.clabel(c1, fmt = '%i')
    plt.clabel(c2, fmt = '%i')
    plt.clabel(c3, fmt = '%i')
    plt.clabel(c4, fmt = '%i')
    plt.plot([0], [0], color1)
    plt.plot([0], [0], color2)
    plt.plot([0], [0], color3)

    plt.legend(
        [ r'C3 ( $\dfrac{km^2}{s^2}$)',
            r'$V_{\infty}\; (\dfrac{km}{s})$',
            r'Time of Flight (days)'],
        bbox_to_anchor = (1.005, 1.01),
        fontsize = 10)

    ax.set_title(_config['title'], fontsize = _config['fontsize'])
    plt.show()
    '''
    if _config['show']:
        
    if _config['filename'] is not None:
        plt.savefig(_config['filename'], dpi = _config['dpi'])
        print(f"Saved {_config['filename']}")
    '''
    plt.close()

    # Delta V plot

    fig, ax = plt.subplots(figsize = _config['figsize'])
    c0 = ax.contour(
        dv_shorts, levels = _config['dv_levels'],
        cmap = _config['dv_cmap'], linewidths = lw)
    c1 = ax.contour(
        dv_longs, levels=_config['dv_levels'],
        cmap=_config['dv_cmap'], linewidths=lw)
    c2 = ax.contour(tofs, levels = _config['tof_levels'], colors = 'c', linewidths = lw * 0.6)
    plt.clabel(c0, fmt='%.1f')
    plt.clabel(c1, fmt='%.1f')
    plt.clabel(c2, fmt='%i')


    ax.set_title(_config['title'], fontsize = _config['fontsize'])
    ax.set_ylabel('Arrival (Days past %s)' % _config['arrival0'], fontsize = _config['fontsize'])
    ax.set_xlabel('Departure (Days past %s)' % _config['departure0'], fontsize = _config['fontsize'])

    if _config['show']:
        plt.show()

    plt.close()


'''
ECLIPSE_MAP = {
	'umbra'   : ( ( 1,  3 ), ( -1, -3 ) ),
	'penumbra': ( ( 2, -1 ), (  1, -2 ) ),
	'either'  : ( ( 3,  2 ), ( -2, -3 ) )
}
'''

def esc_v( r, mu = pd.earth[ 'mu' ] ):
	'''
	Calculate escape velocity at given radial distance from body
	'''
	return math.sqrt( 2 * mu / r )
 
def state2coes( state, args = {} ):
	_args = {
		'et'        : 0,
		'mu'        : pd.earth[ 'mu' ],
		'deg'       : True,
		'print_coes': False
	}
	for key in args.keys():
		_args[ key ] = args[ key ]

	rp,e,i,raan,aop,ma,t0,mu,ta,a,T = spice.oscltx( 
		state, _args[ 'et' ], _args[ 'mu' ] )

	if _args[ 'deg' ]:
		i    *= nt.r2d
		ta   *= nt.r2d
		aop  *= nt.r2d
		raan *= nt.r2d

	if _args[ 'print_coes' ]:
		print( 'a'   , a    )
		print( 'e'   , e    )
		print( 'i'   , i    )
		print( 'RAAN', raan )
		print( 'AOP' , aop  )
		print( 'TA'  , ta   )
		print()

	return [ a, e, i, ta, aop, raan ]

def state2period( state, mu = pd.earth['mu'] ):

	# specific mechanical energy
	epsilon = nt.norm( state[ 3:6 ] ) ** 2 / 2.0 - mu / nt.norm( state[ :3 ] )

	# semi major axis
	a = -mu / ( 2.0 * epsilon )
 	
	# period
	return 2 * math.pi * math.sqrt( a ** 3 / mu )

def coes2state( coes, mu = pd.earth[ 'mu' ], deg = True ):
	a, e, i, ta, aop, raan = coes
	if deg:
		i    *= nt.d2r
		ta   *= nt.d2r
		aop  *= nt.d2r
		raan *= nt.d2r

	rp = a * ( 1 - e )

	return spice.conics( [ rp, e, i, raan, aop, ta, 0, mu], 0 )

def state2ap( state, mu = pd.earth[ 'mu' ] ):
	h       = nt.norm( np.cross( state[ :3 ], state[ 3: ] ) )
	epsilon = nt.norm( state[ 3: ] ) ** 2 / 2.0 - mu / nt.norm( state[ :3 ] )
	e       = math.sqrt( 2 * epsilon * h ** 2 / mu ** 2 + 1 )
	a       = h ** 2 / mu / ( 1 - e ** 2 )
	ra      = a * ( 1 + e )
	rp      = a * ( 1 - e )
	return  ra, rp

def two_body_ode( t, state, mu = pd.earth[ 'mu' ] ):
	# state = [ rx, ry, rz, vx, vy, vz ]

	r = state[ :3 ]
	a = -mu * r / np.linalg.norm( r ) ** 3

	return np.array( [
		state[ 3 ], state[ 4 ], state[ 5 ],
		    a[ 0 ],     a[ 1 ],     a[ 2 ] ] )

def calc_close_approach( turn_angle, v_inf, mu = pd.sun[ 'mu' ] ):
	'''
	Calculate periapsis distance in flyby trajectory
	'''
	return mu * ( 1 / math.sin( turn_angle ) - 1 ) / v_inf ** 2

def calc_vinfinity( tof, args ):

	r1_planet1 = spice.spkgps( args[ 'planet1_ID' ],
		args[ 'et0' ] + tof, args[ 'frame' ], args[ 'center_ID' ] )[ 0 ]

	v0_sc_depart, v1_sc_arrive = lt.lamberts_universal_variables(
		args[ 'state0_planet0' ][ :3 ], r1_planet1, tof,
		{ 'mu': args[ 'mu' ], 'tm': args[ 'tm' ] } )

	vinf = nt.norm( v0_sc_depart - args[ 'state0_planet0' ][ 3: ] )
	return args[ 'vinf' ] - vinf

def vinfinity_match( planet0, planet1, v0_sc, et0, tof0, args = {} ):
	'''
	Given an incoming v-infinity vector to planet0, calculate the
	outgoing v-infinity vector that will arrive at planet1 after
	time of flight (tof) where the incoming and outgoing v-infinity
	vectors at planet0 have equal magnitude
	'''
	_args = {
		'et0'       : et0,
		'planet1_ID': planet1,
		'frame'     : 'ECLIPJ2000',
		'center_ID' : 0,
		'mu'        : pd.sun[ 'mu' ],
		'tm'        : 1,
		'diff_step' : 1e-3,
		'tol'       : 1e-4
	}
	for key in args.keys():
		_args[ key ] = args[ key ]

	_args[ 'state0_planet0' ] = spice.spkgeo( planet0, et0,
		_args[ 'frame' ], _args[ 'center_ID' ] )[ 0 ]

	_args[ 'vinf' ] = nt.norm( v0_sc - _args[ 'state0_planet0' ][ 3: ] )

	tof, steps = nt.newton_root_single_fd(
		calc_vinfinity, tof0, _args )

	r1_planet1 = spice.spkgps( planet1, et0 + tof,
		_args[ 'frame' ], _args[ 'center_ID' ] )[ 0 ]

	v0_sc_depart, v1_sc_arrive = lt.lamberts_universal_variables(
		_args[ 'state0_planet0' ][ :3 ], r1_planet1, tof,
		{ 'mu': _args[ 'mu' ], 'tm': _args[ 'tm' ] } )

	return tof, v0_sc_depart, v1_sc_arrive

'''
def check_eclipse( et, r, body, frame = 'J2000', r_body = 0 ):
	r_sun2body  = spice.spkpos(
		str( body[ 'SPICE_ID' ] ), et, frame, 'LT', 'SUN' )[ 0 ]
	delta_ps    = nt.norm( r_sun2body )
	s_hat       = r_sun2body / delta_ps
	proj_scalar = np.dot( r, s_hat )

	if proj_scalar <= 0.0:
		return -1

	proj     = proj_scalar * s_hat
	rej_norm = nt.norm( r - proj )

	if r_body == 0:
		if check_umbra( delta_ps, body[ 'diameter' ], proj_scalar, rej_norm, r_body ):
			return 2
		elif check_penumbra( delta_ps, body[ 'diameter' ], proj_scalar, rej_norm, r_body ):
			return 1
		else:
			return -1

def check_solar_eclipse_latlons( et, body0, body1, frame = 'J2000' ):
	r_sun2body = spice.spkpos(
		str( body0[ 'SPICE_ID' ] ), et, frame, 'LT', 'SUN' )[ 0 ]
	r = spice.spkpos(
		str( body1[ 'SPICE_ID' ] ), et, frame, 'LT',
		str( body0[ 'SPICE_ID' ] ) )[ 0 ]

	delta_ps    = nt.norm( r_sun2body )
	s_hat       = r_sun2body / delta_ps
	proj_scalar = np.dot( r, s_hat )

	if proj_scalar <= 0.0:
		return -1, None

	proj     = proj_scalar * s_hat
	rej_norm = nt.norm( r - proj )
	umbra    = check_umbra( delta_ps, body0[ 'diameter' ],
		proj_scalar, rej_norm, body1[ 'radius' ] )

	if umbra:
		args = { 'r': r, 's_hat': s_hat, 'radius': body1[ 'radius' ] }
		try:
			sigma = nt.newton_root_single_fd( eclipse_root_func,
				proj_scalar - body1[ 'radius' ], args )[ 0 ]
		except RuntimeError:
			return -1, None

		r_eclipse = sigma * s_hat - r
		r_bf      = np.dot(
			spice.pxform( frame, body1[ 'body_fixed_frame' ], et ),
			r_eclipse )
		latlon        = np.array( spice.reclat( r_bf ) )
		latlon[ 1: ] *= nt.r2d

		return 2, latlon
	else:
		return -1, None

def calc_solar_eclipse_latlons( ets, body0, body1, frame = 'J2000' ):
	eclipses = np.zeros( len( ets ) )
	latlons  = []

	for n in range( len( ets ) ):
		eclipse = check_solar_eclipse(
			ets[ n ], body0, body1, frame )
		if eclipse[ 0 ] == 2:
			latlons.append( eclipse[ 1 ] )
	return np.array( latlons )

def eclipse_root_func( sigma, args ):
	return nt.norm( sigma * args[ 's_hat' ] - args[ 'r' ] ) - args[ 'radius' ]

def check_umbra( delta_ps, Dp, proj_scalar, rej_norm, r_body = 0 ):
	Xu     = ( Dp * delta_ps ) / ( pd.sun[ 'diameter' ] - Dp )
	alphau = math.asin( Dp / ( 2 * Xu ) )
	zeta   = ( Xu - proj_scalar ) * math.tan( alphau )
	return rej_norm - r_body <= zeta

def check_penumbra( delta_ps, Dp, proj_scalar, rej_norm, r_body = 0 ):
	Xp     = ( Dp * delta_ps ) / ( pd.sun[ 'diameter' ] + Dp )
	alphap = math.asin( Dp / ( 2 * Xp ) )
	kappa  = ( Xp + proj_scalar ) * math.tan( alphap )
	return rej_norm - r_body <= kappa

def calc_eclipse_array( ets, rs, body, frame = 'J2000', r_body = 0 ):
	eclipses = np.zeros( rs.shape[ 0 ] )

	for n in range( len( ets ) ):
		eclipses[ n ] = check_eclipse(
			ets[ n ], rs[ n ], body, frame, r_body )

	return eclipses

def find_eclipses( ets, a, method = 'either', v = False, vv = False ):
	diff          = np.diff( a )
	idxs          = ECLIPSE_MAP[ method ]
	ecl_entrances = np.where( np.isin( diff, idxs[ 0 ] ) )[ 0 ]
	ecl_exits     = np.where( np.isin( diff, idxs[ 1 ] ) )[ 0 ]

	if len( ecl_entrances ) == 0:
		return {}

	if ecl_entrances[ 0 ] > ecl_exits[ 0 ]:
		ecl_entrances = np.insert( ecl_entrances, 0, 0 )

	if len( ecl_entrances ) > len( ecl_exits ):
		ecl_exits = np.append( ecl_exits, len( ets ) - 1 )

	ecls                = {}
	ecls[ 'idxs' ]      = []
	ecls[ 'ets'  ]      = []
	ecls[ 'durations' ] = []
	for pair in zip( ecl_entrances, ecl_exits ):
		_ets = [ ets[ pair[ 0 ] ], ets[ pair[ 1 ] ] ]
		ecls[ 'idxs'      ].append( pair )
		ecls[ 'ets'       ].append( _ets )
		ecls[ 'durations' ].append( _ets[ 1 ] - _ets[ 0 ] )

	ecls[ 'total_time' ] = sum( ecls[ 'durations' ] )
	ecls[ 'max_time'   ] = max( ecls[ 'durations' ] )
	ecls[ 'ratio'      ] = ecls[ 'total_time' ] / ( ets[ -1 ] - ets[ 0 ] )

	if v or vv:
		print( '\n******** ECLIPSE SUMMARY START ********' )
		print( f'Number of eclipses: {len(ecls["idxs"])}' )
		print( 'Eclipse durations (seconds): ', end = '' )
		print( [ float(f'{a:.2f}') for a in ecls[ "durations" ] ] )
		print( f'Max eclipse duration: {ecls["max_time"]:.2f} seconds' )
		print( f'Eclipse time ratio: {ecls["ratio"]:.3f}' )
		if vv:
			print( 'Eclipse entrances and exits:' )
			for n in range( len( ecls[ 'ets' ] ) ):
				print(
					spice.et2utc( ecls[ 'ets' ][ n ][ 0 ], 'C', 1 ),
					'-->',
					spice.et2utc( ecls[ 'ets' ][ n ][ 1 ], 'C', 1 )
				)
		print( '******** ECLIPSE SUMMARY END ********\n' )

	return ecls 
'''