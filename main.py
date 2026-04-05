from __future__ import division
import random
from globals import *
from galaxy import Galaxy
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
        rng=rng
    )
    andromeda = Galaxy(
        num_stars=NUM_STARS_ANDROMEDA,
        pos=vector(3, 0, 0) * DIST_SCALE,
        vel=vector(0, 3, 0),
        radius=MAX_ORBITAL_RADIUS,
        thickness=ANDROMEDA_GALAXY_THICKNESS,
        color=vector(0, 0.5, 1),
        rng=rng
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

        if not collision_happened and andromeda.pos.mag < 1.1920057081525512e+20:
            collision_happened = True

        for star in milky_way.stars:
            star.vel += accel(star, andromeda) * dt
            star.vel += accel(star, milky_way) * dt
            star.pos += star.vel * dt
            if andromeda.pos.mag < 1.1920057081525512e+20:
                star.obj.color = vector(1, 0.5, 0)

        for star in andromeda.stars:
            star.vel += accel(star, milky_way) * dt
            star.vel += accel(star, andromeda) * dt
            star.pos += star.vel * dt
            if andromeda.pos.mag < 1.1920057081525512e+20:
                star.obj.color = vector(1, 0.5, 0)

        milky_way.vel += accel(milky_way, andromeda) * dt
        milky_way.pos += milky_way.vel * dt
        andromeda.vel += accel(andromeda, milky_way) * dt
        andromeda.pos += andromeda.vel * dt

        if collision_happened:
            track_step += 1

            if track_step % TRACK_INTERVAL == 0:
                merged = [s for s in milky_way.stars
                          if s.obj.color == vector(1, 0.5, 0)]
                merged += [s for s in andromeda.stars
                           if s.obj.color == vector(1, 0.5, 0)]

                centroid, _, cluster_radius = find_galaxy_centroid(
                    merged, sample_radius_scene=cluster_radius
                )
                if centroid is not None:
                    target_center = centroid

            if target_center is not None:
                if galaxy_center is None:
                    galaxy_center = target_center
                else:
                    galaxy_center += (target_center - galaxy_center) * LERP_SPEED
                scene.center = galaxy_center

if __name__ == '__main__':
    main()