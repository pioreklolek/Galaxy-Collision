"""
liczenie zmian gestosci w galaktykach 

Grid geometry
-------------
  - GRID_N      : cells per axis (default 3  → 3×3×3 = 27 cells)
  - GRID_HALF   : half-size of the whole grid in simulation units
                  default = 2 × MAX_ORBITAL_RADIUS
                  → each cell is (4 × MAX_ORBITAL_RADIUS / GRID_N) on a side
  - Centre-cell (1,1,1) is always aligned with the galaxy/cluster centre.

CSV columns
-----------
  sim_id, step, target, cx, cy, cz,
  grid_half, cell_size,
  ix, iy, iz,          ← cell index 0..GRID_N-1
  cell_cx, cell_cy, cell_cz,   ← world-space centre of that cell
  star_count,
  collision_phase
"""

import csv
import os
import numpy as np

DENSITY_FILE     = "sim_density.csv"
DENSITY_INTERVAL = 10          # liczymy co 10 stepow
GRID_N           = 3           #liczba komorek , np 3 = 3x3x3 ze srodkiem w srodku galaktki/klastra
# grid half-size as a multiple of MAX_ORBITAL_RADIUS (set after import)
GRID_HALF_FACTOR = 1.5         # badamy obszar 1.5x wiekszy od  wielosci poczatkowej galaktyki 

DENSITY_FIELDS = [
    "sim_id", "step", "target",
    "cx", "cy", "cz",           # centre of this grid (galaxy / cluster)
    "grid_half", "cell_size",
    "ix", "iy", "iz",
    "cell_cx", "cell_cy", "cell_cz",
    "star_count",
    "collision_phase",
]


class DensityGrid:
    def __init__(self, sim_id: str, max_orbital_radius: float):
        self._sim_id    = sim_id
        self._grid_half = max_orbital_radius * GRID_HALF_FACTOR
        self._cell_size = (self._grid_half * 2) / GRID_N
        self._phase     = "pre"

        file_is_new = (
            not os.path.exists(DENSITY_FILE)
            or os.path.getsize(DENSITY_FILE) == 0
        )
        self._file   = open(DENSITY_FILE, "w", newline="")
        self._writer = csv.DictWriter(self._file, fieldnames=DENSITY_FIELDS)
        if file_is_new:
            self._writer.writeheader()

        print(
            f"DensityGrid: ready — grid {GRID_N}³, "
            f"half={self._grid_half:.3e}, cell={self._cell_size:.3e}"
        )

    # ── public API ────────────────────────────────────────────────────────────

    def set_phase(self, phase: str):
        self._phase = phase

    def record_step(self, step: int, centres: dict, all_stars: list, phase: str = None):

        if phase is not None:
            self._phase = phase

        pos_array = np.array(
            [[s.pos.x, s.pos.y, s.pos.z] for s in all_stars],
            dtype=np.float64,
        )

        for label, centre in centres.items():
            self._record_target(step, label, centre, pos_array)

        self._file.flush()

    def close(self):
        self._file.flush()
        self._file.close()
        print(f"DensityGrid: closed — data written to '{DENSITY_FILE}'")


    def _record_target(self, step, label, centre, pos_array):
        """Write 27 rows (one per cell) for a single galaxy/cluster centre."""
        cx, cy, cz      = centre.x, centre.y, centre.z
        half            = self._grid_half
        cell            = self._cell_size
        n               = GRID_N

        #liczenie pozycji gwaizdy od centrum
        rel = pos_array - np.array([cx, cy, cz])

        #sprwadzanie ktore gwaizdy sa w boxie
        inside_mask = (
            (rel[:, 0] >= -half) & (rel[:, 0] < half) &
            (rel[:, 1] >= -half) & (rel[:, 1] < half) &
            (rel[:, 2] >= -half) & (rel[:, 2] < half)
        )
        rel_inside = rel[inside_mask]          # shape (M, 3)

        ix_arr = np.floor((rel_inside[:, 0] + half) / cell).astype(int)
        iy_arr = np.floor((rel_inside[:, 1] + half) / cell).astype(int)
        iz_arr = np.floor((rel_inside[:, 2] + half) / cell).astype(int)

        ix_arr = np.clip(ix_arr, 0, n - 1)
        iy_arr = np.clip(iy_arr, 0, n - 1)
        iz_arr = np.clip(iz_arr, 0, n - 1)

        counts = np.zeros((n, n, n), dtype=int)
        np.add.at(counts, (ix_arr, iy_arr, iz_arr), 1)

        # 1 rzad csv dla komorki
        for ix in range(n):
            for iy in range(n):
                for iz in range(n):
                    cell_cx = cx + (-half + (ix + 0.5) * cell)
                    cell_cy = cy + (-half + (iy + 0.5) * cell)
                    cell_cz = cz + (-half + (iz + 0.5) * cell)

                    self._writer.writerow({
                        "sim_id":          self._sim_id,
                        "step":            step,
                        "target":          label,
                        "cx":              f"{cx:.4e}",
                        "cy":              f"{cy:.4e}",
                        "cz":              f"{cz:.4e}",
                        "grid_half":       f"{half:.4e}",
                        "cell_size":       f"{cell:.4e}",
                        "ix":              ix,
                        "iy":              iy,
                        "iz":              iz,
                        "cell_cx":         f"{cell_cx:.4e}",
                        "cell_cy":         f"{cell_cy:.4e}",
                        "cell_cz":         f"{cell_cz:.4e}",
                        "star_count":      int(counts[ix, iy, iz]),
                        "collision_phase": self._phase,
                    })