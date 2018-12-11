#!/usr/bin/python3

from random import random, seed
from trace import *

print("-Demo Animation 2-")

trace = Trace('demo4.png', sky_color=(0, 0, 0))
animation = Animation(
    trace,
    start_frame=0,
    end_frame=59,
    out_file="demoAnimation2.mp4", )
scene = Object(trace, [Location(0, 0, -10)])

seed(1)

for i in range(0, 100):
    r = random()
    g = random()
    b = random()
    length = math.sqrt(r * r + g * g + b * b)
    r /= length
    g /= length
    b /= length
    scale = random()
    rotation = Rotation(random() * math.pi, random() * math.pi, random() * math.pi)
    Cube(trace, [
        Location(
            (random() - random()) * 5,
            (random() - random()) * 5,
            (random() - random()) * 5),
        rotation,
        Scale(scale, scale, scale)
        ],
        material=Material(trace, (r, g, b), reflectiveness=0.25), parent=scene)
    animation.animators.append(
        LinearAnimator(rotation, 'z',
        	[LinearKeyframe(0, 0),
        	LinearKeyframe(60, math.pi)]))

Light(trace)

animation.render()
