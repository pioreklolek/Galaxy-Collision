from vpython import vector, color, sqrt, sphere, rate, scene
from math import fsum
import numpy as np
from numpy import sum as numpy_sum

RECORD_MODE = True
CSV_RECORD_MODE = True
SIM_ID = "sim_7"
SIM_DESCRIPTION = "top-down collision, andromeda falls along Y axis"
CSV_TRACK_INTERVAL = 50

INTEGRATOR = "leapfrog"

MAX_DT_STEP = 20000

# CONSTANTS
RANDOM_SEED = 9999

# Universal gravitational constant
G = 6.673e-11


scene.width = 1300
scene.height = 650

# Solar mass in kg (assume average stellar mass)
SOLAR_MASS = 2.000e30

# Precalculated bounds to solar mass
MIN_SOLAR_MASS = SOLAR_MASS * 0.5
MAX_SOLAR_MASS = SOLAR_MASS * 250
AVG_SOLAR_MASS = SOLAR_MASS * 3.0

# Scale distances for galactic scales
DIST_SCALE = 1e20  # 1e20

CSV_DENSITY_RADIUS = DIST_SCALE * 1.5

# Galactic parameters
MAX_ORBITAL_RADIUS = DIST_SCALE * 15
MIN_ORBITAL_RADIUS = DIST_SCALE * 0.15

MILKY_WAY_GALAXY_THICKNESS = DIST_SCALE * 1.5
ANDROMEDA_GALAXY_THICKNESS = DIST_SCALE * 0.5

COLLISION_THRESHOLD = DIST_SCALE * 5.0

# Milky Way contains about 300 billion stars
NUM_STARS_MILKY_WAY = 700 #700
# Andromeda Galaxy contains about 1 trillion (10^12) stars
NUM_STARS_ANDROMEDA = 1400 #1400

#Earth-Like tracking 
EARTH_ORBITAL_RADIUS = DIST_SCALE * 2.46
EARTH_ORBITAL_BAND = DIST_SCALE * 0.15 # +- tolerance 
MAX_EARTH_ANALOGS = 6 # max num , of tracked stars

TRACK_EARTH_IN_MILKY_WAY = True
TRACK_EARTH_IN_ANDROMEDA = False

EARTH_COLOR = vector(0,1,0)
TURN_ON_LOGS = True

TRACK_INTERVAL = 50
LERP_SPEED = 0.05


# Graphical constants
STAR_RADIUS = 0.025
dt = 1e16 # 1e16 or 1e17


# FUNCTIONS

# Limit x between lower and upper
def clamp(x, lower, upper):
    return max(min(x, upper), lower)


# Return the force due to gravity on an object
def gravity(mass1, mass2, radius):
    return G * mass1 * mass2 / radius**2


 # Return the acceleration due to gravity on an object.
def g_accel(mass, radius):
    eps = MIN_ORBITAL_RADIUS
    return G * mass / (radius**2 + eps**2)


# # old g-accel
# def g_accel(mass, radius):
#     radius = max(radius, MIN_ORBITAL_RADIUS)
#     return G * mass / radius / radius

# Calculate acceleration on an object caused by galaxy
def accel(obj, galaxy):
    r_galaxy = galaxy.pos - obj.pos
    # We have a = F / m = G * m_center / r ^2
    return r_galaxy.norm() * g_accel(galaxy.mass, r_galaxy.mag)
