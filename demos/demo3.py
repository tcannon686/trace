#!/usr/bin/python3
from random import random, seed
from raytrace import *

def run(title, out_file, render_basic=False, code_path=None, render_window=False):
	print(title)

	t = Trace(out_file, sky_color=(0, 0, 0), basic=render_basic, code_path=code_path, render_window=render_window)

	scene = Object(t, [Location(0, 0, -10)])

	seed(1)

	for i in range(0, 100):
		r = random()
		g = random()
		b = random()
		length = math.sqrt(r * r + g * g + b * b)
		r /= length
		g /= length
		b /= length
		scale = random() * 0.5
		Cube(t, [
		    Location(
		        (random() - random()) * 5,
		        (random() - random()) * 5,
		        (random() - random()) * 5),
		    Rotation(random() * math.pi, random() * math.pi, random() * math.pi),
		    Scale(scale, scale, scale)
		    ],
		    material=Material(t, (r, g, b), reflectiveness=0.25), parent=scene)

	Light(t)

	t.trace()
