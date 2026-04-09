from __future__ import division
import random
from globals import *
from galaxy import Galaxy
from tracker import StarTracker
import numpy as np


def main():
    scene.objects.clear()

    rng = random.Random(RANDOM_SEED)
    print(rng.gauss(0, 1))
    
    milky_way = Galaxy(
        num_stars=NUM_STARS_MILKY_WAY,
        pos=vector(-9, 0, 0) * DIST_SCALE,
        vel=vector(0, 0, 0),
        radius=MAX_ORBITAL_RADIUS,
        thickness=MILKY_WAY_GALAXY_THICKNESS,
        color=vector(0.9, 0.9, 1),
        rng=rng
    )
    andromeda = Galaxy(
        num_stars=NUM_STARS_ANDROMEDA,
        pos=vector(6, 0, 0) * DIST_SCALE,
        vel=vector(0, 3, 0),
        radius=MAX_ORBITAL_RADIUS,
        thickness=ANDROMEDA_GALAXY_THICKNESS,
        color=vector(0, 0.5, 1),
        rng=rng
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

    collision_happened = False
    track_step = 0
    
    while True:
        rate(100)

        galaxy_separation = (andromeda.pos - milky_way.pos).mag
        if not collision_happened and galaxy_separation < COLLISION_THRESHOLD:
            collision_happened = True

            for star in milky_way.stars:
                if not star.is_earth_analog:
                    star.obj.color = vector(1, 0.5, 0)
                star.is_merged = True
            for star in andromeda.stars:
                if not star.is_earth_analog:
                    star.obj.color = vector(1, 0.5, 0)
                star.is_merged = True

        #leapfrog gravity algorithm
        # half - kick
        for star in milky_way.stars:
            a = accel(star,milky_way) + accel(star, andromeda)
            star.vel += 0.5 * a * dt
        for star in andromeda.stars:
            a = accel(star,milky_way) + accel(star, andromeda)
            star.vel += 0.5 * a * dt
        a_mw = accel(milky_way,andromeda)
        a_and = accel(andromeda,milky_way)
        milky_way.vel += 0.5 * a_mw * dt
        andromeda.vel += 0.5 * a_and * dt

        #drift 
        for star in milky_way.stars:
            star.pos += star.vel * dt
        for star in andromeda.stars:
            star.pos += star.vel * dt
        milky_way.pos += milky_way.vel * dt
        andromeda.pos += andromeda.vel * dt

        #half kick 
        for star in milky_way.stars:
            a = accel(star,milky_way) + accel(star,andromeda)
            star.vel += 0.5 * a * dt
        for star in andromeda.stars:
            a = accel(star,milky_way) + accel(star,andromeda)
            star.vel += 0.5 * a * dt

        a_mw = accel(milky_way, andromeda)
        a_and = accel(andromeda, milky_way)
        milky_way.vel += 0.5 * a_mw * dt
        andromeda.vel += 0.5 * a_and * dt
        
        if collision_happened:
            track_step += 1
            merged = [s for s in milky_way.stars if s.is_merged]
            merged += [s for s in andromeda.stars if s.is_merged]

            camera_center = tracker.record_step(step=track_step, merged_stars=merged)
            if camera_center is not None:
                scene.center = camera_center
        
if __name__ == '__main__':
    main()