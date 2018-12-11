#!/usr/bin/python3

from trace import *

print("-Demo Animation 1-")

trace = Trace(sky_color=(0.25, 0, 0.25), quality=32)

rotation = Rotation(0, math.pi / 8, 0)
scene = Object(trace, [Location(0, 0, -8), rotation])
Sphere(trace, material=Material(trace, (1, 0, 0)), parent=scene)
Plane(trace, [Location(0, -1, 0), Scale(4, 4, 4)], material=Material(trace, (0, 0, 1), reflectiveness=0.5), parent=scene)

for i in range(0, 3):
    Light(trace, [Rotation(i * 2 * math.pi / 3, 0, 0), Location(0, 4, -5)])


animation = Animation(
	trace,
	[
		LinearAnimator(
			rotation, 'x',
			[
				LinearKeyframe(0, 1),
				LinearKeyframe(60, 1 + 2 * math.pi)
			])
	],
	start_frame=0,
	end_frame = 59,
	out_file="demoAnimation1.mp4")

animation.render()

