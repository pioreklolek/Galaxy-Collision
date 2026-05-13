from app.globals import * 

class Gravity_calc:

    @staticmethod
    def step_euler(milky_way, andromeda, dt):
        for star in milky_way.stars:
            star.vel += accel(star, andromeda) * dt
            star.vel += accel(star, milky_way) * dt
            star.pos += star.vel * dt
        for star in andromeda.stars:
            star.vel += accel(star, milky_way) * dt
            star.vel += accel(star, andromeda) * dt
            star.pos += star.vel * dt
        milky_way.vel += accel(milky_way, andromeda) * dt
        milky_way.pos += milky_way.vel * dt
        andromeda.vel += accel(andromeda, milky_way) * dt
        andromeda.pos += andromeda.vel * dt

    @staticmethod
    def step_leapfrog(milky_way, andromeda, dt):
        # half-kick
        for star in milky_way.stars:
            star.vel += 0.5 * (accel(star, milky_way) + accel(star, andromeda)) * dt
        for star in andromeda.stars:
            star.vel += 0.5 * (accel(star, milky_way) + accel(star, andromeda)) * dt
        milky_way.vel += 0.5 * accel(milky_way, andromeda) * dt
        andromeda.vel += 0.5 * accel(andromeda, milky_way) * dt

        # drift
        for star in milky_way.stars:
            star.pos += star.vel * dt
        for star in andromeda.stars:
            star.pos += star.vel * dt
        milky_way.pos += milky_way.vel * dt
        andromeda.pos += andromeda.vel * dt

        # half-kick
        for star in milky_way.stars:
            star.vel += 0.5 * (accel(star, milky_way) + accel(star, andromeda)) * dt
        for star in andromeda.stars:
            star.vel += 0.5 * (accel(star, milky_way) + accel(star, andromeda)) * dt
        milky_way.vel += 0.5 * accel(milky_way, andromeda) * dt
        andromeda.vel += 0.5 * accel(andromeda, milky_way) * dt