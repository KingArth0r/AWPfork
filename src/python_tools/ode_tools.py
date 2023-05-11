'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Ordinary Differential Equations (ODEs) Tools Library
'''

# Python standard libraries

# 3rd party libraries

# AWP library

def euler_step(f ,t ,y ,h):
	return y + h*f(t, y)

# Define the acceleration function for a two-body problem
def acceleration(q, mu):
    r_norm = np.linalg.norm(q)
    return -mu * q / r_norm**3

# Symplectic Verlet integrator function
def verlet_step(q0, v0, mu, total_time, time_step):
    num_steps = int(total_time / time_step)

    q = np.zeros((num_steps + 1, len(q0)))
    v = np.zeros((num_steps + 1, len(q0)))

    q[0] = q0
    v[0] = v0

    for n in range(num_steps):
        # Update the position using the current velocity
        q[n + 1] = q[n] + time_step * v[n]

        # Compute the acceleration at the updated position
        a = acceleration(q[n + 1], mu)

        # Update the velocity using the acceleration at the updated position
        v[n + 1] = v[n] + time_step * a

    return q, v

def rk4_step( f, t, y, h ):
	'''
	Calculate one RK4 step
	'''
	k1 = f( t, y )
	k2 = f( t + 0.5 * h, y + 0.5 * k1 * h )
	k3 = f( t + 0.5 * h, y + 0.5 * k2 * h )
	k4 = f( t +       h, y +       k3 * h )

	return y + h / 6.0 * ( k1 + 2 * k2 + 2 * k3 + k4 )


methods = {
	'rk4': rk4_step,
	'euler': euler_step,
	'verlet': verlet_step
}



def explicit_midpoint(f, t, y, h):
    k1 = f(t,y)
    k2 = f(t + 1/2 * h, y + 1/2 * h * k1)

    return y + h(k2)

