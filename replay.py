import pickle
from vpython import sphere, vector, rate, scene

STAR_RADIUS = 0.025
scene.width = 1300
scene.height = 650
REPLAY_SPEED = 5 # lower slower

with open("simulation.pkl", "rb") as f:
    data = pickle.load(f)

frames = data["frames"]
camera_centers = data["camera_centers"]
num_stars = len(frames[0])

spheres = [sphere(pos=vector(0,0,0), radius=STAR_RADIUS) for _ in range(num_stars)]

print(f"Replaying {len(frames)} frames, {num_stars} stars...")


for frame, cam in zip(frames, camera_centers):  
    rate(REPLAY_SPEED)
    cx, cy, cz = cam
    scene.center = vector(cx, cy, cz)       
    for i, (x, y, z, cr, cg, cb) in enumerate(frame):
        spheres[i].pos = vector(x, y, z)
        spheres[i].color = vector(cr, cg, cb)


    