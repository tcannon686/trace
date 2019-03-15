#!/usr/bin/python3
"""
This file runs all of the demo files.
"""

import sys
sys.path.append('demos/')

def dsimple():
	import simple
	simple.run('Simple', 'simple.png', render_basic, 'simple.ray' if output_code else None, render_window)

def d1():
	import demo1
	demo1.run('Demo 1', 'demo1.png', render_basic, 'demo1.ray' if output_code else None, render_window)

def d2():
	import demo2
	demo2.run('Demo 2', 'demo2.png', render_basic, 'demo2.ray' if output_code else None, render_window)

def d3():
	import demo3
	demo3.run('Demo 3', 'demo3.png', render_basic, 'demo3.ray' if output_code else None, render_window)

def d4():
	import demo4
	demo4.run('Demo 4', 'demo4.png', render_basic, 'demo4.ray' if output_code else None, render_window)

def da1():
	import demoAnimation1
	demoAnimation1.run('Demo Animation 1', 'demoAnimation1.mp4', render_basic, 'demoAnimation1.ray' if output_code else None, render_window)

def da2():
	import demoAnimation2
	demoAnimation2.run('Demo Animation 2', 'demoAnimation2.mp4', render_basic, 'demoAnimation2.ray' if output_code else None, render_window)

def da3():
	import demoAnimation3
	demoAnimation3.run('Demo Animation 3', 'demoAnimation3.mp4', render_basic, 'demoAnimation3.ray' if output_code else None, render_window)

demos = {
	"demo1" : d1,
	"demo2" : d2,
	"demo3" : d3,
	"demo4" : d4,
	"demoAnimation1" : da1,
	"demoAnimation2" : da2,
	"demoAnimation3" : da3,
	"simple" : dsimple
}

demos_sorted = demos.keys()

i = 0
start_index = 1
render_basic = False
output_code = False
render_window = False

if start_index == len(sys.argv):
	print("Usage: " + sys.argv[0] + " " + "[basic] [code] [window] [all] [list] {[demo1] ... [demon]}")
else:
	if sys.argv[start_index] == 'basic':
		render_basic = True
		start_index += 1

	if sys.argv[start_index] == 'code':
		output_code = True
		start_index += 1

	if sys.argv[start_index] == 'window':
		render_window = True
		start_index += 1

	if sys.argv[start_index] == 'all':
		for name in demos_sorted:
			demos[name]()
		start_index += 1
	if sys.argv[start_index] == 'list':
		print('The following demos are available:')
		for name in demos_sorted:
			print('    ' + name)
		start_index += 1

	for key in sys.argv:
		if i >= start_index:
			if key in demos:
				demos[key]()
			else:
				print("No demo \"" + key + "\"!")
		i += 1




