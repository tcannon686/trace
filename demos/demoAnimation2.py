
from raytrace import *

def run(title, out_file, render_basic=False, code_path=None, render_window=False):
	print(title)
	t = Trace(
		out_file=out_file,
		sky_color=(1, 1, 1),
		render_threads=4, render_samples=4, shadow_samples=10,
		render_section_size=128,
		render_iterations=4,
		quality=24,
		width=800, height=600,
		ambient_color=(0.25, 0.25, 0.25),
		ao_samples=10,
		basic=render_basic,
		code_path=code_path,
		render_window=render_window)

	scene = Object(t, [Location(0, 0, -6), Rotation(0, math.pi / 8, 0), Rotation(math.pi / 16, 0, 0)])


	rot = Rotation()
	rotator = Object(t, [rot], parent=scene)
	colors = [
		    (1, 0, 0),
		    (1, 1, 1),
		    (0, 0, 1)
		]
	cindex = 0
	for i in range(-2, 3):
		Cube(
			t,
			[
				Location(0, i * 4, 0),
				Scale(0.25 / 2, 4 / 2, 0.25 / 2)
			],
			Material(t, diffuse=(0.25, 0.25, 0.25)),
			parent=rotator)


		for j in range(-4, 4):
			Cube(
				t,
				[
					Location(0, j / 2 + i * 4, 0),
					Rotation(j * math.pi / 4, 0, 0),
					Scale(1.0, 0.1, 0.1),
				],
				
				Material(t, diffuse=colors[cindex%len(colors)]),
				parent=rotator)
			cindex += 1

	loc = Location(0, 4, 1)
	spheres = Object(t, [loc], parent=scene)
	for i in range(-4, 4):
		Sphere(t, [
			Location(0, i * 2, 0),
			Rotation(1, 1, 1),
			Scale(0.25, 0.25, 0.25)],
			Material(t, diffuse=(0, 0.25, 1), specular=(1, 1, 1), shininess=20,
				alpha=0.2,
				ior=1.5,
				reflectiveness=0.5), parent=spheres)

	Light(t, [Location(4, 4, 4)])

	animation = Animation(
		t,
		[
			LinearAnimator(rot, "x", [
				LinearKeyframe(0, 0),
				LinearKeyframe(240, 4 * math.pi)
			]),
			LinearAnimator(loc, "y", [
				LinearKeyframe(0, 4),
				LinearKeyframe(240, -4)
			])
		],
		out_file=out_file,
		start_frame=0,
		end_frame=239)
	if render_window:
		t.trace()
	else:
		animation.render()


