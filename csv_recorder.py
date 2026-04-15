import csv
import os
import math
from globals import *

META_FILE = "sim_meta.csv"
DATA_FILE = "sim_data.csv"

META_FIELDS = [
    "sim_id", "description",
    "random_seed", "num_stars_mw", "num_stars_andromeda",
    "dt", "earth_orbital_radius", "earth_orbital_band",
    "max_earth_analogs", "track_mw", "track_andromeda",
    "mw_initial_pos", "mw_initial_vel",
    "and_initial_pos", "and_initial_vel",
    "collision_threshold", "max_orbital_radius",
    "mw_thickness", "and_thickness"
]

DATA_FIELDS = [
    "sim_id", "step", "star_id", "home_galaxy",
    # Position
    "pos_x", "pos_y", "pos_z",
    # Velocity
    "vel_x", "vel_y", "vel_z", "speed",
    # Distance
    "dist_from_mw_center", "dist_from_and_center", "dist_from_nearest_center",
    # Orbit health
    "orbital_plane_tilt_deg",   # angle between angular momentum and Y axis
                                # ~0 = flat orbit, ~90 = vertical orbit
    "eccentricity_proxy",       # current radius / initial radius
                                # ~1 = stable, >> 1 or << 1 = distorted
    # Environment
    "local_star_density",       # num stars within CSV_DENSITY_RADIUS
    "radial_velocity",          # velocity component pointing away from nearest center
                                # positive = moving outward, negative = falling in
    # Binding
    "is_bound",                 # 1 = gravitationally bound to nearest galaxy, 0 = ejected
    # Phase
    "collision_phase"           # "pre", "collision", "post"
]

class CsvRecorder:
    def __init__(self, sim_id, description, milky_way, andromeda):
        self._sim_id = sim_id
        self._mw = milky_way
        self._and = andromeda
        self._phase = "pre"

        self._write_meta(description, milky_way, andromeda)

        data_file_is_new = not os.path.exists(DATA_FILE) or os.path.getsize(DATA_FILE) == 0
        self._data_file = open(DATA_FILE, "a", newline="")
        self._writer = csv.DictWriter(self._data_file, fieldnames=DATA_FIELDS)

        if data_file_is_new:
            self._writer.writeheader()

        print(f"CsvRecorder: ready — sim_id='{sim_id}'")

    def set_phase(self, phase):
        """ Call with 'pre', 'collision', or 'post' from main.py """
        self._phase = phase

    def record_step(self, step, earth_analogs, all_stars):
        """
        earth_analogs : list of Star — the green Earth-like stars from StarTracker
        all_stars     : list of Star — every star, used for density calculation
        """
        if step % CSV_TRACK_INTERVAL != 0:
            return

        # Build a fast numpy array of all star positions for density lookup
        all_pos = self._build_pos_array(all_stars)

        for idx, star in enumerate(earth_analogs):
            row = self._compute_row(step, idx, star, all_pos)
            self._writer.writerow(row)

    def close(self):
        self._data_file.flush()
        self._data_file.close()
        print(f"CsvRecorder: closed — data written to '{DATA_FILE}'")
    

    def _compute_row(self, step, star_id, star, all_pos):
        mw_pos  = self._mw.pos
        and_pos = self._and.pos

        # --- distances ---
        dist_mw  = (star.pos - mw_pos).mag
        dist_and = (star.pos - and_pos).mag
        nearest_center_pos = mw_pos if dist_mw <= dist_and else and_pos
        dist_nearest = min(dist_mw, dist_and)

        # --- speed ---
        vel = star.vel
        speed = vel.mag

        # --- orbital plane tilt ---
        # Angular momentum vector L = r x v (relative to nearest center)
        r_vec = star.pos - nearest_center_pos
        # vpython vectors: use cross product
        L = r_vec.cross(vel)
        if L.mag > 0:
            y_axis = vector(0, 1, 0)
            cos_angle = L.dot(y_axis) / L.mag
            cos_angle = max(-1.0, min(1.0, cos_angle))   # clamp for safety
            tilt_deg = math.degrees(math.acos(abs(cos_angle)))
            # 0° = orbit flat in XZ plane (healthy), 90° = vertical orbit
        else:
            tilt_deg = 0.0

        # --- eccentricity proxy ---
        current_radius = (star.pos - nearest_center_pos).mag
        init_r = star.initial_orbital_radius if star.initial_orbital_radius > 0 else 1.0
        eccentricity_proxy = current_radius / init_r

        # --- local star density ---
        star_scene = star.pos / DIST_SCALE
        star_np = [star_scene.x, star_scene.y, star_scene.z]
        density_radius_scene = CSV_DENSITY_RADIUS / DIST_SCALE
        d2 = ((all_pos - star_np) ** 2).sum(axis=1)
        # subtract 1 to exclude the star itself
        local_density = int(numpy_sum(d2 < density_radius_scene ** 2)) - 1

        # --- radial velocity ---
        r_norm = r_vec.norm()
        radial_vel = vel.dot(r_norm)   # positive = moving away from center

        # --- is bound ---
        # Bound if total specific energy < 0: KE + PE < 0
        nearest_mass = self._mw.mass if dist_mw <= dist_and else self._and.mass
        ke = 0.5 * speed ** 2
        pe = -G * nearest_mass / (dist_nearest + 1e-10)
        is_bound = 1 if (ke + pe) < 0 else 0

        return {
            "sim_id":                   self._sim_id,
            "step":                     step,
            "star_id":                  f"{star.home_galaxy_label}_star_{star_id}",
            "home_galaxy":              star.home_galaxy_label,
            "pos_x":                    f"{star.pos.x:.4e}",
            "pos_y":                    f"{star.pos.y:.4e}",
            "pos_z":                    f"{star.pos.z:.4e}",
            "vel_x":                    f"{vel.x:.4e}",
            "vel_y":                    f"{vel.y:.4e}",
            "vel_z":                    f"{vel.z:.4e}",
            "speed":                    f"{speed:.4e}",
            "dist_from_mw_center":      f"{dist_mw:.4e}",
            "dist_from_and_center":     f"{dist_and:.4e}",
            "dist_from_nearest_center": f"{dist_nearest:.4e}",
            "orbital_plane_tilt_deg":   f"{tilt_deg:.2f}",
            "eccentricity_proxy":       f"{eccentricity_proxy:.4f}",
            "local_star_density":       local_density,
            "radial_velocity":          f"{radial_vel:.4e}",
            "is_bound":                 is_bound,
            "collision_phase":          self._phase,
        }
    
    def _build_pos_array(self, all_stars):
        import numpy as np
        return np.array([
            [s.obj.pos.x, s.obj.pos.y, s.obj.pos.z]
            for s in all_stars
        ])
    
    def _write_meta(self, description, milky_way, andromeda):
        file_exists = os.path.exists(META_FILE) and os.path.getsize(META_FILE) > 0
        with open(META_FILE, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=META_FIELDS)
            if not file_exists:
                writer.writeheader()
            writer.writerow({
                "sim_id":               self._sim_id,
                "description":          description,
                "random_seed":          RANDOM_SEED,
                "num_stars_mw":         NUM_STARS_MILKY_WAY,
                "num_stars_andromeda":  NUM_STARS_ANDROMEDA,
                "dt":                   dt,
                "earth_orbital_radius": EARTH_ORBITAL_RADIUS,
                "earth_orbital_band":   EARTH_ORBITAL_BAND,
                "max_earth_analogs":    MAX_EARTH_ANALOGS,
                "track_mw":             TRACK_EARTH_IN_MILKY_WAY,
                "track_andromeda":      TRACK_EARTH_IN_ANDROMEDA,
                "mw_initial_pos":       f"{milky_way.pos}",
                "mw_initial_vel":       f"{milky_way.vel}",
                "and_initial_pos":      f"{andromeda.pos}",
                "and_initial_vel":      f"{andromeda.vel}",
                "collision_threshold":  COLLISION_THRESHOLD,
                "max_orbital_radius":   MAX_ORBITAL_RADIUS,
                "mw_thickness":         MILKY_WAY_GALAXY_THICKNESS,
                "and_thickness":        ANDROMEDA_GALAXY_THICKNESS,
            })