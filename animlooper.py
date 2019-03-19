#!/usr/bin/python3

# This program will make an animation loop using ffmpeg.

import sys
import subprocess
import os

if len(sys.argv) != 4:
	print("Usage: " + sys.argv[0] + " <repeats> <in_file> <out_file>")
else:

	l = open('list.txt', 'w')
	for i in range(0, int(sys.argv[1])):
		l.write('file ')
		l.write(sys.argv[2])
		l.write('\n')
	l.close()

	subprocess.run(
		['ffmpeg', '-f', 'concat', '-i', 'list.txt', '-c', 'copy', sys.argv[3]])

	os.remove('list.txt')

