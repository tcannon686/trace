#ifdef INCLUDE_GUI

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <X11/Xlib.h>
#include <unistd.h>

#include "matrix.h"
#include "hashtable.h"
#include "trace.h"
#include "xwin.h"

#define NIL (0)

typedef struct
{
	Display *display;
	Colormap cmap;
	GC gc;
	Window window;
} callback_data_t;

void RenderWindowCallback(int x, int y, vector_t color, void *data)
{
	callback_data_t *d = (callback_data_t *)data;
	/*XColor xcolor;
	xcolor.red = 0xFFFF * color.x;
	xcolor.green = 0xFFFF * color.y;
	xcolor.blue = 0xFFFF * color.z;
	xcolor.flags = DoRed | DoGreen | DoBlue;
	XAllocColor(d->display, d->cmap, &xcolor);*/
	
	unsigned int c =
		((int) (color.x * 0xFF) << 16) |
		((int) (color.y * 0xFF) << 8) |
		(int) (color.z * 0xFF);
	
	XSetForeground(d->display, d->gc, c);
	XDrawPoint(d->display, d->window, d->gc, x, y);
	
	XFlush(d->display);
}

int WindowMain(render_settings_t *rs, int width, int height)
{
	XInitThreads();
	Display *display = XOpenDisplay(NULL);
	
	
	int blackColor = BlackPixel(display, DefaultScreen(display));
	int whiteColor = WhitePixel(display, DefaultScreen(display));
	
	Window window = XCreateSimpleWindow(
		display, DefaultRootWindow(display), 0, 0, width, height,
		0, blackColor, blackColor);
	
	XSelectInput(display, window, StructureNotifyMask);
	
	XMapWindow(display, window);
	
	Atom wmDeleteMessage = XInternAtom(display, "WM_DELETE_WINDOW", False);
	XSetWMProtocols(display, window, &wmDeleteMessage, 1);
	
	Colormap cmap = DefaultColormap(display, DefaultScreen(display));
	
	GC gc = XCreateGC(display, window, 0, NIL);
	XSetForeground(display, gc, whiteColor);
	
	vecc_t aspect = (vecc_t)height / (vecc_t)width;
	vecc_t scale_x = 2 / (vecc_t)width;
	vecc_t scale_y = -aspect * 2 / (vecc_t)height;
	
	printf("info: generating tree\n");
	CleanTriangles(&rs->triangles_ptr);
	kd_tree_t *tree_ptr = GenerateTree(rs->triangles_ptr, 0, NULL);
	
	callback_data_t cd;
	cd.display = display;
	cd.cmap = cmap;
	cd.gc = gc;
	cd.window = window;
	
	Render(
		RenderWindowCallback,
		(void *)&cd,
		width, height,
		tree_ptr,
		rs->lights_ptr,
		rs->sky,
		rs->section_size,
		aspect,
		scale_x, scale_y,
		rs->focal_length,
		rs->samples,
		rs->shadow_samples,
		rs->max_iterations,
		rs->threads);
	
	XEvent e;
	while(1)
	{
		XNextEvent(display, &e);
		if(e.type == MapNotify)
		{
			
		}
		else if(e.xclient.data.l[0] == wmDeleteMessage)
		{
			break;
		}
	}
	
	FreeTree(tree_ptr);
	XDestroyWindow(display, window);
	XCloseDisplay(display);
	return 1;
}


#endif

