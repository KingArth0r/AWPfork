'''
main.py --- Primary tester file.

This is the tester that generates the time data for our plots.
This is computationally expensive - it took Jimmy's machine about 3 hours to run.

Attribution - Jimmy
'''
import orbit_calculations as ob
import lamberthub as lh
import numpy as np
import time

X = np.linspace(.25, 1, 20)


# test gooding 1990
Y = np.zeros(20)

for k in range(len(X)):
    print('Gooding iteration:', k)
    tick = time.perf_counter()

    ob.interplanetary_porkchop({'step': 86400 * X[k], 'show': False, "integrator": lh.gooding1990})

    tock = time.perf_counter()
    Y[k] = tock - tick # time passed

np.savetxt("gooding1990.csv", [X,Y], delimiter=",")

# test Avanzini 2008

for k in range(len(X)):
    print('Avanzini iteration:', k)
    tick = time.perf_counter()

    ob.interplanetary_porkchop({'step': 86400 * X[k], 'show': False, "integrator": lh.avanzini2008})

    tock = time.perf_counter()
    Y[k] = tock - tick # time passed

np.savetxt("avanzini2008.csv", [X,Y], delimiter=",")

Y = np.zeros(20)

# test Arora 2013

for k in range(len(X)):
    print('Arora iteration:', k)
    tick = time.perf_counter()

    ob.interplanetary_porkchop({'step': 86400 * X[k], 'show': False, "integrator": lh.arora2013})

    tock = time.perf_counter()
    Y[k] = tock - tick # time passed

np.savetxt("arora2013.csv", [X,Y], delimiter=",")

Y = np.zeros(20)

# test Izzo 2015

for k in range(len(X)):
    print('Izzo iteration:', k)
    tick = time.perf_counter()

    ob.interplanetary_porkchop({'step': 86400 * X[k], 'show': False, "integrator": lh.izzo2015})

    tock = time.perf_counter()
    Y[k] = tock - tick # time passed

np.savetxt("izzo2015.csv", [X,Y], delimiter=",")