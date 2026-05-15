import csv
import os
import numpy as np
from app.globals import *

HABITABILITY_FILE = "sim_habitability.csv"

HABITABILITY_FIELDS = [
    "sim_id",
    "scan_phase",           # "initial" - przed zderzeniem  lub "final" - po zderzeniu
    "step",
    "galaxy",               
    "total_stars",

    #kryterium 1 - tylko odlegosc od centrum , niezaleznie od osi
    "habitable_by_radius",
    "habitable_by_radius_pct",

    #kryterium 2 - odleglosc od centrum oraz ograniczenie plaszczyzny oś Y
    "habitable_by_radius_and_plane",
    "habitable_by_radius_and_plane_pct",

    "earth_orbital_radius",
    "earth_orbital_band",
    "plane_band_y",         # maksymalne odchylenie osi Y 
]


class HabitabilityRecorder:
    #Skanuje galaktyki pod kątem % gwiazd zdolnych do zamieszkania.

    PLANE_BAND_Y: float = EARTH_ORBITAL_BAND

    def __init__(self, sim_id: str):
        self._sim_id = sim_id
        self._file_is_new = (
            not os.path.exists(HABITABILITY_FILE)
            or os.path.getsize(HABITABILITY_FILE) == 0
        )
        self._file = open(HABITABILITY_FILE, "w", newline="")
        self._writer = csv.DictWriter(self._file, fieldnames=HABITABILITY_FIELDS)
        if self._file_is_new:
            self._writer.writeheader()
        print(f"HabitabilityRecorder: ready — sim_id='{sim_id}'")

    # ------------------------------------------------------------------
    def scan(self, galaxies, step: int, phase: str):
        #phase    – "initial" lub "final"
        
        r_min = EARTH_ORBITAL_RADIUS - EARTH_ORBITAL_BAND
        r_max = EARTH_ORBITAL_RADIUS + EARTH_ORBITAL_BAND
        y_max = self.PLANE_BAND_Y

        for galaxy, label in galaxies:
            stars = galaxy.stars
            total = len(stars)
            if total == 0:
                continue

            cx, cy, cz = galaxy.pos.x, galaxy.pos.y, galaxy.pos.z

            positions = np.array([
                [s.pos.x - cx, s.pos.y - cy, s.pos.z - cz]
                for s in stars
            ])

            # Odleglosc 3D od centrum galaktyki
            dist3d = np.sqrt((positions ** 2).sum(axis=1))

            #kryt 1 
            mask_radius = (dist3d >= r_min) & (dist3d <= r_max)

            #kryt 2 
            mask_plane = mask_radius & (np.abs(positions[:, 1]) <= y_max)

            count_radius = int(mask_radius.sum())
            count_plane  = int(mask_plane.sum())

            pct_radius = count_radius / total * 100.0
            pct_plane  = count_plane  / total * 100.0

            row = {
                "sim_id":                           self._sim_id,
                "scan_phase":                       phase,
                "step":                             step,
                "galaxy":                           label,
                "total_stars":                      total,
                "habitable_by_radius":              count_radius,
                "habitable_by_radius_pct":          f"{pct_radius:.4f}",
                "habitable_by_radius_and_plane":    count_plane,
                "habitable_by_radius_and_plane_pct": f"{pct_plane:.4f}",
                "earth_orbital_radius":             f"{EARTH_ORBITAL_RADIUS:.4e}",
                "earth_orbital_band":               f"{EARTH_ORBITAL_BAND:.4e}",
                "plane_band_y":                     f"{y_max:.4e}",
            }
            self._writer.writerow(row)

            if TURN_ON_LOGS:
                print(
                    f"[Habitability | {phase} | step {step}] "
                    f"{label}: total={total} | "
                    f"by_radius={count_radius} ({pct_radius:.2f}%) | "
                    f"by_radius+plane={count_plane} ({pct_plane:.2f}%)"
                )

        self._file.flush()

    def close(self):
        self._file.flush()
        self._file.close()
        print(f"HabitabilityRecorder: closed — data written to '{HABITABILITY_FILE}'")