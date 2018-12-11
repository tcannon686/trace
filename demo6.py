#!/usr/bin/python3
"""
This is a demo of the so-called Menger Sponge
(https://en.wikipedia.org/wiki/Menger_sponge).
This demo takes a considerable amount of time to 
generate the code for the renderer, but have 
patience, it doesn't take too long to render. The
total process should take around a minute.
"""

from trace import *

print("-Demo 6-")

trace = Trace(
    sky_color=(0.0, 0.0, 0.0),
    out_file="demo6.png",
)

material = Material(
    trace,
    (1.0, 0.0, 1.0)
)

def sponge(trace, maxiters=1, parent=None, iteration=0):
    if parent == None:
        parent = Object(transforms=[Scale(1/3, 1/3, 1/3)])
    
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                if (i == 0 and j != 0 and k != 0)\
                    or (j == 0 and i != 0 and k != 0)\
                    or (k == 0 and i != 0 and j != 0)\
                    or (i != 0 and j != 0 and k != 0):
                    if iteration < maxiters:
                        sponge(
                            trace,
                            maxiters, 
                            Object(
                                trace,
                                transforms=[
                                    Location(i, j, k),
                                    Scale(1/3, 1/3, 1/3)
                                ],
                                parent=parent),
                            iteration + 1)
                    else:
                        Cube(
                            trace,
                            transforms=[
                                Location(i, j, k)
                            ],
                            material=material,
                            parent=parent)
    return parent
scene = Object(trace, [Location(0, 0, -6)])
carpet = Object(trace, [Rotation(0, math.pi / 4, 0), Rotation(math.pi / 4, 0, 0)], scene)

Light(trace, [Location(0, 0, 4)], parent=scene)

sponge(trace, 3, carpet)

trace.trace()

