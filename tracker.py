from globals import *
import random

class StarTracker:
    #select stars whose orbits at step 0 is close to Earth's , marks them green

    def __init__(self, galaxies, earth_radius, earth_band):
        #galaxies - list of Galaxy obj to search
        # earth radius - oribtal radius of Earth-Like (meters), import from Globals
        #earth_band - tolerance +- to orbital radius (meters) , import from globals 
        self.earth_analogs = []
        self._select(galaxies,earth_radius,earth_band)
        print(f"StarTracker: {len(self.earth_analogs)} Earth-Like stars selected!")

    def _select(self,galaxies,earth_radius,earth_band):
        candidates = []
        for galaxy in galaxies:
            for star in galaxy.stars:
                orbital_radius = (star.pos - galaxy.pos).mag
                if abs(orbital_radius - earth_radius) < earth_band:
                    candidates.append(star)
            
        if MAX_EARTH_ANALOGS is not None and len(candidates) > MAX_EARTH_ANALOGS:
            candidates = random.sample(candidates,MAX_EARTH_ANALOGS)

        for star in candidates:
            star.obj.color = EARTH_COLOR
            star.is_earth_analog = True
            self.earth_analogs.append(star)

    def record_step(self,step,galaxy_center):
        pass

    def save_to_file(self,filepath,sim_params):
        pass