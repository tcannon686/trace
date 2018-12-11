#!/usr/bin/python3
from trace import *

print("-Demo 3-")

trace = Trace(out_file="demo3.png", sky_color=(0.25, 0, 0.25), quality=32)

scene = Object(trace, [Location(0, 0, -8), Rotation(math.pi / 4, math.pi / 8, 0)])
Sphere(trace, material=Material(trace, (1, 0, 0)), parent=scene)
Plane(trace, [Location(0, -1, 0), Scale(4, 4, 4)], material=Material(trace, (0, 0, 1), reflectiveness=0.5), parent=scene)

for i in range(0, 3):
    Light(trace, [Rotation(i * 2 * math.pi / 3, 0, 0), Location(0, 4, -5)])

trace.trace()
