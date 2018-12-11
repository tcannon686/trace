#!/usr/bin/python3
from trace import *

print("-Demo 2-")

trace = Trace(out_file="demo2.png", sky_color=(0, 0, 0), render_iterations=10, render_samples=4)

scene = Object(trace, [Location(0, 0, -10), Rotation(1, 0, 0)])

for i in range(0, 4):
    for j in range(0, 4):
        for k in range(0, 4):
            Sphere(trace,
                [Location(i - 4 / 2 + 1 / 2, -1 + k, j - 4 / 2 + 1 / 2),
                    Scale(0.5, 0.5, 0.5)],
                material=Material(trace, reflectiveness=1.0), parent=scene)
Plane(trace, [Location(0, -1.5, 0), Scale(10, 10, 10)], parent=scene)

Light(trace, [Location(0, 5, 10)], parent=scene)
trace.trace()
