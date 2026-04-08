from __future__ import division
import random
from globals import *
from galaxy import Galaxy
from tracker import StarTracker
import numpy as np


def find_galaxy_centroid(stars, sample_radius_scene=5.0):
    if len(stars) == 0:
        return None, 0, sample_radius_scene

    positions = np.array([[s.pos.x / DIST_SCALE,
                           s.pos.y / DIST_SCALE,
                           s.pos.z / DIST_SCALE] for s in stars])

    r2 = sample_radius_scene ** 2
    neighbor_counts = np.zeros(len(stars), dtype=int)
    for i in range(len(positions)):
        diffs = positions - positions[i]
        dist2 = (diffs ** 2).sum(axis=1)
        neighbor_counts[i] = np.sum(dist2 < r2)

    seed_idx = np.argmax(neighbor_counts)
    seed_pos = positions[seed_idx]

    diffs = positions - seed_pos
    dist2 = (diffs ** 2).sum(axis=1)
    cluster_positions = positions[dist2 < r2]
    centroid = cluster_positions.mean(axis=0)

    spreads = np.sqrt(((cluster_positions - centroid) ** 2).sum(axis=1))
    next_radius = float(np.clip(spreads.mean() * 2.5, 3.0, 30.0))

    return vector(centroid[0], centroid[1], centroid[2]), len(cluster_positions), next_radius

def calc_ejection_rate(total_stars, stars_in_cluster):
    #return % of stars ejected in collision
    if total_stars == 0:
        return 0
    return (total_stars - stars_in_cluster) / total_stars * 100.0

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
    TRACK_INTERVAL = 50
    cluster_radius = 5.0
    galaxy_center = None

    target_center = None
    LERP_SPEED = 0.05

    
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

            if track_step % TRACK_INTERVAL == 0:
                merged = [s for s in milky_way.stars
                          if s.is_merged]
                merged += [s for s in andromeda.stars
                           if s.is_merged]

                centroid, cluster_count, cluster_radius = find_galaxy_centroid(
                    merged, sample_radius_scene=cluster_radius
                )
                if centroid is not None:
                    target_center = centroid

                #logs 
                ejection_percent = calc_ejection_rate(len(merged), cluster_count)

                if TURN_ON_LOGS:
                    ejected = len(merged) - cluster_count
                    print(f"[step {track_step:>6}] "
                          f"cluster: {cluster_count:>4} stars | "
                          f"ejected: {ejected:>4} stars | "
                          f"ejection rate: {ejection_percent:>5.1f}%")



            if target_center is not None:
                if galaxy_center is None:
                    galaxy_center = target_center
                else:
                    galaxy_center += (target_center - galaxy_center) * LERP_SPEED
                scene.center = galaxy_center
        
        tracker.record_step(step=track_step,galaxy_center=galaxy_center)
        
if __name__ == '__main__':
    main()