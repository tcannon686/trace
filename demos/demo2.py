#!/usr/bin/python3
from raytrace import *

def run(title, out_file, render_basic=False, code_path=None, render_window=False):

	print(title)

	t = Trace(out_file=out_file, sky_color=(0.25, 0, 0.25), quality=32,
		basic=render_basic, code_path=code_path, render_window=render_window)

	scene = Object(t, [Location(0, 0, -8), Rotation(math.pi / 4, math.pi / 8, 0)])
	Sphere(t, material=Material(t, (1, 0, 0)), parent=scene)
	Plane(t, [Location(0, -1, 0), Scale(2, 2, 2)], material=Material(t, (0, 0, 1), reflectiveness=0.5), parent=scene)

	for i in range(0, 3):
		Light(t, [Rotation(i * 2 * math.pi / 3, 0, 0), Location(0, 4, -5)])

	t.trace()
