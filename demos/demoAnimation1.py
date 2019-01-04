#!/usr/bin/python3
from raytrace import *

def run(title, out_file, render_basic=False, code_path=None, render_window=False):

	print(title)

	t = Trace(out_file=out_file, sky_color=(1, 1, 1),
		ambient_color=(0.5, 0.5, 0.5),
		render_iterations=5, render_samples=4, shadow_samples=1, ao_samples=4,
		basic=render_basic, code_path=code_path, render_window=render_window)

	rot = Rotation(0, 0, 0)
	scene = Object(t,
		[Location(0, 0, -4),
		Rotation(0, math.pi / 8, 0),
		rot])

	count = 3

	base = Material(t, reflectiveness=.5)
	mat = Material(t, diffuse=(1, 0, 0), shininess=1, specular=(0, 0, 0), reflectiveness=0.0)
	shiny = Material(t, diffuse=(1, 1, 1), shininess=100, ambient=(0, 0, 0), alpha=0.0, ior=1.5)

	for i in range(0, count):
		for j in range(0, count):
		        Sphere(t,
		        [Location(i - count / 2 + 1 / 2, 0.5 if i == count // 2 and j == count // 2 else 0, j - count / 2 + 1 / 2),
		            Scale(0.5, 0.5, 0.5)],
		        material=shiny if i == count // 2 and j == count // 2 else mat, parent=scene)
	Plane(t, [Location(0, -0.5, 0), Scale(10000, 10000, 10000)], material=base, parent=scene)

	Light(t, [Location(0, 10, 0)], parent=scene)

	a = Animation(
		t,
		[
			LinearAnimator(
				rot, 'x',
				[LinearKeyframe(0, 0),
				LinearKeyframe(120, 2 * math.pi)])
		],
		start_frame=0,
		end_frame=119,
		out_file=out_file)
	if render_window:
		t.trace()
	else:
		a.render()


