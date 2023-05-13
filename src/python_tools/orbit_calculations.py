# 3rd party libraries
import numpy as np
import spiceypy as spice
import matplotlib.pyplot as plt
import lamberthub as lh

# AWP library
import numerical_tools as nt
import planetary_data as pd
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
        'step': 86400, # default is equal to 1 day (in seconds)
        'frame': 'ECLIPJ2000',
        'observer': 'SOLAR SYSTEM BARYCENTER',
        'cutoff_v': 20.0,
        'c3_levels': None,
        'vinf_levels': None,
        'tof_levels': None,
        'dv_levels': None,
        'dv_cmap': 'RdPu_r',
        'figsize': (10,20),
        'lw': 1.5, # line width
        'title': 'Porkchop Plot',
        'fontsize': 15,
        'show': True,
        'filename': None,    # file name for base plot, if you want it saved
        'filename_dv': None, # file name for delta-v specific plot
        'dpi': 300,
        'load': False, # load the plot while the function is being run?
        'integrator': lh.izzo2015
    }

    # SPICE kernels contain neccacary planetary positional data
    spice.furnsh(sd.leapseconds_kernel)
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
            state_depart = sd.load_ephemeris(
                _config['planet0'],
                [et_departures[nd]],
                _config['frame'], _config['observer'])[0]

            # state of planet1 (Mars) at arrival
            state_arrive = sd.load_ephemeris(
                _config['planet1'],
                [et_arrivals[na]],
                _config['frame'], _config['observer'])[0]

            # calculate time of flight
            tof = et_arrivals[na] - et_departures[nd]

            mu = _config['mu']

            # Try lambert's solver. If lambert's solver fails, exception is thrown
            try:
                v_sc_depart_short, v_sc_arrive_short = _config['integrator'](_config['mu'], state_depart[:3],
                                                        state_arrive[ :3 ], tof, prograde= True, full_output= False)
            except Exception as e:
                v_sc_depart_short = np.array([1000, 1000, 1000])
                v_sc_arrive_short = np.array([1000, 1000, 1000])
            try:
                v_sc_depart_long, v_sc_arrive_long = _config['integrator'](_config['mu'], state_depart[:3],
                                                        state_arrive[:3], tof, prograde=False, full_output=False)
            except Exception as e:
                v_sc_depart_long = np.array([1000, 1000, 1000])
                v_sc_arrive_long = np.array([1000, 1000, 1000])


            C3_short = np.linalg.norm(v_sc_depart_short - state_depart[3:])**2
            C3_long = np.linalg.norm(v_sc_depart_long - state_depart[3:])**2
            if C3_short > cutoff_c3: C3_short = cutoff_c3
            if C3_long > cutoff_c3: C3_long = cutoff_c3

            v_inf_short = np.linalg.norm(v_sc_arrive_short - state_arrive[3:])
            v_inf_long = np.linalg.norm(v_sc_arrive_long - state_arrive[3:])

            if v_inf_short > _config['cutoff_v']: v_inf_short = _config['cutoff_v']
            if v_inf_long > _config['cutoff_v']: v_inf_long = _config['cutoff_v']

            C3_shorts[na, nd] = C3_short
            C3_longs[na, nd] = C3_long
            v_inf_shorts[na, nd] = v_inf_short
            v_inf_longs[na,nd] = v_inf_long
            tofs[na, nd] = tof


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

    # labels and plots
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


    if _config['show']:
        plt.show()

    if _config['filename'] is not None:
        plt.savefig(_config['filename'], dpi = _config['dpi'])
        print(f"Saved {_config['filename']}")
    plt.close()


    # Delta V plot

    fig, ax = plt.subplots(figsize = _config['figsize'])

    # contours
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
    if _config['filename_dv'] is not None:
        plt.savefig(_config['filename_dv'], dpi = _config['dpi'])
        print(f"Saved {_config['filename_dv']}")
    plt.close()

