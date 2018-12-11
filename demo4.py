#!/usr/bin/python3
from random import random, seed
from trace import *

print("-Demo 4-")

trace = Trace('demo4.png', sky_color=(0, 0, 0))

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
    Cube(trace, [
        Location(
            (random() - random()) * 5,
            (random() - random()) * 5,
            (random() - random()) * 5),
        Rotation(random() * math.pi, random() * math.pi, random() * math.pi),
        Scale(scale, scale, scale)
        ],
        material=Material(trace, (r, g, b), reflectiveness=0.25), parent=scene)

Light(trace)

trace.trace()
