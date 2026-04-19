from globals import *
import numpy as np
import random

class StarTracker:
    #select stars whose orbits at step 0 is close to Earth's , marks them green

    def __init__(self, galaxies, earth_radius, earth_band):
        #galaxies - list of Galaxy obj to search
        # earth radius - oribtal radius of Earth-Like (meters), import from Globals
        #earth_band - tolerance +- to orbital radius (meters) , import from globals 
        self._cluster_radius = 5.0
        self._galaxy_center = None
        self._target_center = None
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

    def _find_two_centroids(self, stars, sample_radius_scene):
        if len(stars) == 0:
            return None, 0, None, 0, sample_radius_scene 

        positions = np.array([[s.pos.x / DIST_SCALE,
                            s.pos.y / DIST_SCALE,
                            s.pos.z / DIST_SCALE] for s in stars])

        r2 = sample_radius_scene ** 2

        def find_cluster(pos_array):
            counts = np.zeros(len(pos_array), dtype=int)
            for i in range(len(pos_array)):
                d2 = ((pos_array - pos_array[i]) ** 2).sum(axis=1)
                counts[i] = np.sum(d2 < r2)
            seed = pos_array[np.argmax(counts)]
            mask = ((pos_array - seed) ** 2).sum(axis=1) < r2
            cluster = pos_array[mask]
            centroid = cluster.mean(axis=0)
            return centroid, mask, cluster
        
        centroid1, mask1, cluster1 = find_cluster(positions)

        remaining = positions[~mask1]
        if len(remaining) > 0:
            centroid2, _, cluster2 = find_cluster(remaining)
            c2 = vector(centroid2[0], centroid2[1],centroid2[2])
            count2 = len(cluster2)
        else:
            c2, count2 = None, 0

        spreads = np.sqrt(((cluster1 - centroid1) ** 2).sum(axis=1))
        next_radius = float(np.clip(spreads.mean() * 2.5, 3.0, 30.0))

        c1 = vector(centroid1[0], centroid1[1], centroid1[2])
        
        return c1, len(cluster1), c2, count2, next_radius
    
    def _calc_ejection_rate(self, total_stars, count1, count2):
    #return % of stars ejected in collision
        if total_stars == 0:
            return 0
        return (total_stars - count1 - count2) / total_stars * 100.0



    def record_step(self, step, merged_stars):
        if step % TRACK_INTERVAL != 0:
            # ... existing lerp code ...
            if self._target_center is not None:
                if self._galaxy_center is None:
                    self._galaxy_center = self._target_center
                else:
                    self._galaxy_center += (self._target_center - self._galaxy_center) * LERP_SPEED
            return self._galaxy_center, None   # <-- added None

        c1, count1, c2, count2, self._cluster_radius = self._find_two_centroids(
            merged_stars, self._cluster_radius
        )
        if c1 is not None:
            largest = c1 if (c2 is None or count1 >= count2) else c2
            self._target_center = largest

        total = len(merged_stars)
        ejected = total - count1 - count2
        ejection_rate = self._calc_ejection_rate(total, count1, count2)

        if TURN_ON_LOGS:
            print(f"[step {step:>6}] "
                  f"cluster1: {count1:>4} | cluster2: {count2:>4} | "
                  f"ejected: {ejected:>4} | ejection rate: {ejection_rate:>5.1f}%")

        if self._target_center is not None:
            if self._galaxy_center is None:
                self._galaxy_center = self._target_center
            else:
                self._galaxy_center += (self._target_center - self._galaxy_center) * LERP_SPEED

        cluster_stats = {                      # <-- new
            "count1": count1,
            "count2": count2,
            "ejected": ejected,
            "ejection_rate": ejection_rate,
            "total": total,
        }
        return self._galaxy_center, cluster_stats   # <-- added stats

    def save_to_file(self,filepath,sim_params):
        pass