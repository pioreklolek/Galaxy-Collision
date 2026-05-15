from __future__ import division
import random
import pickle
from app.globals import *
from app.galaxy import Galaxy
from app.tracker import StarTracker
from app.csv_recorder import CsvRecorder
from app.gravity_calc import Gravity_calc
from app.density_grid import DensityGrid, DENSITY_INTERVAL
from app.habitability_recorder import HabitabilityRecorder
import numpy as np

HABITABILITY_FINAL_STEP = 19_900


def main():
    scene.objects.clear()

    rng = random.Random(RANDOM_SEED)
    print(rng.gauss(0, 1))
    
    milky_way = Galaxy(
    num_stars=NUM_STARS_MILKY_WAY,
    pos=vector(-5, 0, 0) * DIST_SCALE, 
    vel=vector(0, 0, 0),               
    radius=MAX_ORBITAL_RADIUS,
    thickness=MILKY_WAY_GALAXY_THICKNESS,
    color=vector(0.9, 0.9, 1),
    rng=rng,
    label="milky_way"
    )

    andromeda = Galaxy(
        num_stars=NUM_STARS_ANDROMEDA,
        pos=vector(3, 0, 0) * DIST_SCALE,  
        vel=vector(0, 3, 0),               
        radius=MAX_ORBITAL_RADIUS,
        thickness=ANDROMEDA_GALAXY_THICKNESS,
        color=vector(0, 0.5, 1),
        rng=rng,
        label="andromeda"
    )

    galaxies_to_track = []
    if TRACK_EARTH_IN_MILKY_WAY:
        galaxies_to_track.append(milky_way)
    if TRACK_EARTH_IN_ANDROMEDA:
        galaxies_to_track.append(andromeda)

    tracker = StarTracker(
        galaxies=galaxies_to_track,
        earth_radius=EARTH_ORBITAL_RADIUS,
        earth_band=EARTH_ORBITAL_BAND
    )

    recorder = CsvRecorder(
        sim_id=SIM_ID,
        description=SIM_DESCRIPTION,
        milky_way=milky_way,
        andromeda=andromeda
    )

    density_grid = DensityGrid(
    sim_id=SIM_ID,
    max_orbital_radius=MAX_ORBITAL_RADIUS,
)
    
    galaxies_for_hab = [
        (milky_way, "milky_way"),
        (andromeda, "andromeda"),
    ]
 
    hab_recorder = HabitabilityRecorder(sim_id=SIM_ID)
 
    # SKAN POCZĄTKOWY  habitat 
    hab_recorder.scan(galaxies_for_hab, step=0, phase="initial")


    #recording 
    recording = RECORD_MODE
    
    frames = []
    camera_centres = []

    def snapshot(all_stars):
        return [(s.obj.pos.x, s.obj.pos.y, s.obj.pos.z,
                 s.obj.color.x, s.obj.color.y, s.obj.color.z)
                for s in all_stars]

    all_stars = list(milky_way.stars) + list(andromeda.stars)

    collision_happened = False
    track_step = 0
    final_hab_scanned = False

    last_density_centres = {
    "milky_way": milky_way.pos,
    "andromeda": andromeda.pos,
}

    while True:
        rate(100) #100 lub 1000

        galaxy_separation = (andromeda.pos - milky_way.pos).mag

        if not collision_happened and galaxy_separation < COLLISION_THRESHOLD:
            collision_happened = True
            recorder.set_phase("collision")

            for star in milky_way.stars:
                if not star.is_earth_analog:
                    star.obj.color = vector(1, 0.5, 0)
                star.is_merged = True
            for star in andromeda.stars:
                if not star.is_earth_analog:
                    star.obj.color = vector(1, 0.5, 0)
                star.is_merged = True

        
        #Gravity_calc.step_leapfrog(milky_way,andromeda,dt)
        Gravity_calc.step_euler(milky_way,andromeda,dt)

        
        if collision_happened:
            track_step += 1
            
            if track_step > 1:
                recorder.set_phase("post")

            #wykonaj skan koncowy raz , gdy step dojdzie do HABITABILITY_FINAL_STEP
            if not final_hab_scanned and track_step >= HABITABILITY_FINAL_STEP:
                hab_recorder.scan(galaxies_for_hab, step=track_step, phase="final")
                final_hab_scanned = True


            if track_step >= MAX_DT_STEP:
                print(f"Simulation ended after {track_step} post-collision steps.")
                break
                
            
            merged = [s for s in milky_way.stars if s.is_merged]
            merged += [s for s in andromeda.stars if s.is_merged]

            camera_center, cluster_stats = tracker.record_step(  
                step=track_step, merged_stars=merged
            )
            if camera_center is not None:
                scene.center = camera_center


            #liczenie zmian gestosci:
            if track_step % DENSITY_INTERVAL == 0:

                if (cluster_stats is not None
                        and cluster_stats.get("c1_x") is not None
                        and cluster_stats.get("c2_x") is not None):
                    # post-kolizja  
                    last_density_centres = {
                        "cluster1": vector(
                            cluster_stats["c1_x"] * DIST_SCALE,
                            cluster_stats["c1_y"] * DIST_SCALE,
                            cluster_stats["c1_z"] * DIST_SCALE,
                        ),
                        "cluster2": vector(
                            cluster_stats["c2_x"] * DIST_SCALE,
                            cluster_stats["c2_y"] * DIST_SCALE,
                            cluster_stats["c2_z"] * DIST_SCALE,
                        ),
                    }
                elif not collision_happened:
                    # pre-kolizja galaktyki na swoich pozycjach
                    last_density_centres = {
                        "milky_way": milky_way.pos,
                        "andromeda": andromeda.pos,
                    }

                density_grid.record_step(         
                    step=track_step,
                    centres=last_density_centres, 
                    all_stars=all_stars,
                    phase=recorder._phase,
                )
 
            if CSV_RECORD_MODE:
                recorder.record_step(
                    step=track_step,
                    earth_analogs=tracker.earth_analogs,
                    all_stars=all_stars
                )
                recorder.record_cluster_step(                      # new 
                    step=track_step,
                    stats=cluster_stats
                )
            if recording and track_step % TRACK_INTERVAL == 0:
                frames.append(snapshot(all_stars))
                camera_centres.append((scene.center.x, scene.center.y, scene.center.z))

    if CSV_RECORD_MODE:
        recorder.close()
    density_grid.close()
    hab_recorder.close()

    if recording:
        with open("simulation.pkl", "wb") as f:
            pickle.dump({
                "num_mw": len(milky_way.stars),
                "frames": frames,
                "camera_centers": camera_centres
            }, f)
        print(f"Saved {len(frames)} frames to simulation.pkl")
        
if __name__ == '__main__':
    main()