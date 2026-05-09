from __future__ import division
import random
import pickle
from globals import *
from galaxy import Galaxy
from tracker import StarTracker
from csv_recorder import CsvRecorder
from gravity_calc import Gravity_calc
import numpy as np


def main():
    scene.objects.clear()

    rng = random.Random(RANDOM_SEED)
    print(rng.gauss(0, 1))
    
    milky_way = Galaxy(
    num_stars=NUM_STARS_MILKY_WAY,
    pos=vector(0, -9, 0) * DIST_SCALE, 
    vel=vector(0, 3, 0),               
    radius=MAX_ORBITAL_RADIUS,
    thickness=MILKY_WAY_GALAXY_THICKNESS,
    color=vector(0.9, 0.9, 1),
    rng=rng,
    label="milky_way"
    )
    andromeda = Galaxy(
        num_stars=NUM_STARS_ANDROMEDA,
        pos=vector(0, 6, 0) * DIST_SCALE,  
        vel=vector(0, -3, 0),               
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
    while True:
        rate(100)

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

        
        Gravity_calc.step_leapfrog(milky_way,andromeda,dt)
        #Gravity_calc.step_euler(milky_way,andromeda,dt)

        
        if collision_happened:
            track_step += 1
            
            if track_step > 1:
                recorder.set_phase("post")

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