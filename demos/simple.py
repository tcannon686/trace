#!/usr/bin/python3

from raytrace import *

def run(title, out_file, render_basic=False, code_path=None, render_window=False):
	print(title)

	t = Trace(
		sky_color=(0.5, 0, 0.5),
		out_file=out_file,
		quality=32,
		render_window=render_window,
		basic=render_basic,
		code_path=code_path
	)

	scene = Object(t, [Location(0, 0, -8)])

	Plane(
		t,
		[
			Location(0, -1, 0),
			Scale(2.5, 2.5, 2.5)
		],
		Material(
			t,
			(1, 0, 0),
			(0, 0, 0),
			reflectiveness=0.5
		),
		parent=scene
	)

	Light(
		t,
		[
			Location(0, 5, 5)
		],
		parent=scene
	)

	t.trace()