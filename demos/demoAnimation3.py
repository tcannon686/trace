
from raytrace import *

def run(title, out_file, render_basic=False, code_path=None, render_window=False):

	print(title)
	t = Trace(out_file=out_file, sky_color=(0, 0, 0), quality=32,
		basic=render_basic, code_path=code_path, render_window=render_window)

	rot = Rotation(math.pi, 0, 0)
	Sphere(t, [Location(0, 0, -4), rot], Material(t, tex_diffuse=
		Texture2d(t, ImageFile(t, 'demos/earth.png'))))
	Light(t, [Location(-5, 5, 0)])

	animation = Animation(
		t,
		[
			LinearAnimator(
				rot, 'x',
				[
					LinearKeyframe(0, 0),
					LinearKeyframe(120, math.pi * 2)
				])
		],
		out_file=out_file,
		start_frame=0, end_frame=119)

	
	if render_window:
		t.trace()
	else:
		animation.render()


