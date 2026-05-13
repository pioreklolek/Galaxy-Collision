from app.globals import *
from app.star import Star


class Galaxy(object):
    def __init__(self, num_stars, pos, vel, radius, thickness, color,rng, label=""):
        self.pos = pos
        self.vel = vel
        self.radius = radius
        self.label = label

        # Gaussian distributions
        sigma_mass = AVG_SOLAR_MASS / 3.0
        masses = [clamp(rng.gauss(mu=AVG_SOLAR_MASS, sigma=sigma_mass), MIN_SOLAR_MASS, MAX_SOLAR_MASS)
                  for i in range(num_stars)]

        # Galaxy mass is sum of all stars
        self.mass_without_dark_matter = fsum(masses)

        #dark matter 
        self.mass = self.mass_without_dark_matter * (1.0 + DARK_MATTER_FACTOR)

        # Gaussian distribution of positions
        sigma_x = radius * 0.1
        sigma_y = thickness * 0.10
        sigma_z = radius * 0.1

        # Generate list of all positions
        positions = []
        for i in range(num_stars):
            pos = vector(
                clamp(rng.gauss(mu=0, sigma=sigma_x), -radius, radius),
                clamp(rng.gauss(mu=0, sigma=sigma_y), -thickness, thickness),
                clamp(rng.gauss(mu=0, sigma=sigma_z), -radius, radius)
            )

            # Limit radius to avoid particles shooting to nowhere
            if pos.mag < MIN_ORBITAL_RADIUS:
                pos.mag = MIN_ORBITAL_RADIUS

            positions.append(pos)

        def calc_orbital_velocity(center_mass, radius):
            return sqrt(G * center_mass / radius)

        # Generate list of all stars
        stars = []
        up = vector(0.0, 1.0, 0.0)
        for i in range(num_stars):
            # Find normalized vector along direction of travel
            absolute_pos = positions[i] + self.pos
            relative_pos = positions[i]
            vec = relative_pos.cross(up).norm()
            relative_vel = vec * \
                calc_orbital_velocity(self.mass, relative_pos.mag)
            absolute_vel = relative_vel + vel

            star = Star(
                mass=masses[i],
                radius=STAR_RADIUS,
                pos=absolute_pos,
                vel=absolute_vel,
                color=color
            )
            star.home_galaxy_label = label
            star.initial_orbital_radius = positions[i].mag  
            stars.append(star)

        self.stars = np.array(stars)
