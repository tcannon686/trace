#!/usr/bin/python3
from trace import *

print("-Demo 1-")

trace = Trace(
    sky_color=(0.5, 0, 0.5),
    out_file="demo1.png",
    code_path = "demo1.ray",
    quality=32
)

scene = Object(trace, [Location(0, 0, -8)])
Sphere(
    trace,
    [
        Location(0, 0.5, 0),
        Rotation(1, 0, 0)
    ],
    Material(
        trace,
        (1, 1, 1),
        (1.0, 1.0, 1.0),
        shininess=25,
        reflectiveness=0.5),
    parent=scene
)

Plane(
    trace,
    [
        Location(0, -1, 0),
        Scale(5, 5, 5)
    ],
    Material(
        trace,
        (1, 0, 0),
        (0, 0, 0),
        reflectiveness=0.5
    ),
    parent=scene
)

Light(
    trace,
    [
        Location(0, 5, 5)
    ],
    parent=scene
)

trace.trace()
