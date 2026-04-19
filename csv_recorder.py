import csv
import os
import math
import numpy as np
from globals import *

META_FILE = "sim_meta.csv"
DATA_FILE = "sim_data.csv"
CLUSTER_FILE = "sim_cluster.csv"

META_FIELDS = [
    "sim_id", "description",
    "integrator",
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
    # Distance
    "dist_from_mw_center", "dist_from_and_center", "dist_from_nearest_center",
    # Velocity
    "speed", "radial_velocity", "tangential_velocity",
    # Orbit health
    "orbital_plane_tilt_deg",
    "eccentricity_proxy",
    # Environment
    "local_star_density",
    # Binding
    "specific_energy",          # ujemna = związana, dodatnia = ucieka (zastępuje is_bound)
    "nearest_galaxy",           # "MW" lub "AND" — do której galaktyki liczymy energię
    # Phase
    "collision_phase"
]

CLUSTER_FIELDS = [
    "sim_id", "step",
    "cluster1", "cluster2",
    "ejected", "total_merged",
    "ejection_rate_pct",
    "collision_phase"
]

class CsvRecorder:
    def __init__(self, sim_id, description,milky_way, andromeda,integrator=INTEGRATOR):
        self._sim_id = sim_id
        self._integrator = integrator
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

        cluster_file_is_new = not os.path.exists(CLUSTER_FILE) or os.path.getsize(CLUSTER_FILE) == 0
        self._cluster_file = open(CLUSTER_FILE, "a", newline="")
        self._cluster_writer = csv.DictWriter(self._cluster_file, fieldnames=CLUSTER_FIELDS)
        if cluster_file_is_new:
            self._cluster_writer.writeheader()

    def set_phase(self, phase):
        self._phase = phase

    def record_step(self, step, earth_analogs, all_stars):
        if step % CSV_TRACK_INTERVAL != 0:
            return

        all_pos = self._build_pos_array(all_stars)

        for idx, star in enumerate(earth_analogs):
            row = self._compute_row(step, idx, star, all_pos)
            self._writer.writerow(row)

    def _compute_row(self, step, star_id, star, all_pos):
        mw_pos  = self._mw.pos
        and_pos = self._and.pos

        # --- distances ---
        dist_mw  = (star.pos - mw_pos).mag
        dist_and = (star.pos - and_pos).mag
        nearest_is_mw = dist_mw <= dist_and
        nearest_center_pos = mw_pos if nearest_is_mw else and_pos
        nearest_mass       = self._mw.mass if nearest_is_mw else self._and.mass
        nearest_label      = "MW" if nearest_is_mw else "AND"
        dist_nearest       = min(dist_mw, dist_and)

        # --- velocity decomposition ---
        vel   = star.vel
        speed = vel.mag
        r_vec  = star.pos - nearest_center_pos
        r_norm = r_vec.norm()

        radial_vel      = vel.dot(r_norm)                       # + = ucieka, - = opada
        tangential_vel  = (vel - r_norm * radial_vel).mag       # prędkość prostopadła do r

        # --- orbital plane tilt ---
        L = r_vec.cross(vel)
        if L.mag > 0:
            y_axis    = vector(0, 1, 0)
            cos_angle = max(-1.0, min(1.0, L.dot(y_axis) / L.mag))
            tilt_deg  = math.degrees(math.acos(abs(cos_angle)))
        else:
            tilt_deg = 0.0

        # --- eccentricity proxy ---
        current_radius = r_vec.mag
        init_r         = star.initial_orbital_radius if star.initial_orbital_radius > 0 else 1.0
        eccentricity_proxy = current_radius / init_r

        # --- local star density ---
        star_np = [star.pos.x / DIST_SCALE,
                   star.pos.y / DIST_SCALE,
                   star.pos.z / DIST_SCALE]
        density_radius_scene = CSV_DENSITY_RADIUS / DIST_SCALE
        d2 = ((all_pos - star_np) ** 2).sum(axis=1)
        local_density = int(np.sum(d2 < density_radius_scene ** 2)) - 1

        # --- specific energy (zastępuje is_bound) ---
        ke             = 0.5 * speed ** 2
        pe             = -G * nearest_mass / (dist_nearest + 1e-10)
        specific_energy = ke + pe   # < 0 = związana, > 0 = ucieka

        return {
            "sim_id":                   self._sim_id,
            "step":                     step,
            "star_id":                  f"{star.home_galaxy_label}_star_{star_id}",
            "home_galaxy":              star.home_galaxy_label,
            "dist_from_mw_center":      f"{dist_mw:.4e}",
            "dist_from_and_center":     f"{dist_and:.4e}",
            "dist_from_nearest_center": f"{dist_nearest:.4e}",
            "speed":                    f"{speed:.4e}",
            "radial_velocity":          f"{radial_vel:.4e}",
            "tangential_velocity":      f"{tangential_vel:.4e}",
            "orbital_plane_tilt_deg":   f"{tilt_deg:.2f}",
            "eccentricity_proxy":       f"{eccentricity_proxy:.4f}",
            "local_star_density":       local_density,
            "specific_energy":          f"{specific_energy:.4e}",
            "nearest_galaxy":           nearest_label,
            "collision_phase":          self._phase,
        }

    def _build_pos_array(self, all_stars):
        return np.array([
            [s.pos.x / DIST_SCALE,
             s.pos.y / DIST_SCALE,
             s.pos.z / DIST_SCALE]
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
                "integrator":           self._integrator,
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
    
    def record_cluster_step(self, step, stats):
        """stats dict: count1, count2, ejected, total, ejection_rate"""
        if stats is None:
            return
        self._cluster_writer.writerow({
            "sim_id":            self._sim_id,
            "step":              step,
            "cluster1":          stats["count1"],
            "cluster2":          stats["count2"],
            "ejected":           stats["ejected"],
            "total_merged":      stats["total"],
            "ejection_rate_pct": f"{stats['ejection_rate']:.2f}",
            "collision_phase":   self._phase,
        })

    def close(self):
        self._data_file.flush()
        self._data_file.close()
        self._cluster_file.flush()       # <-- new
        self._cluster_file.close()       # <-- new
        print(f"CsvRecorder: closed — data written to '{DATA_FILE}' and '{CLUSTER_FILE}'")
