#ifdef INCLUDE_GUI

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <unistd.h>
#include <time.h>

#include "matrix.h"
#include "hashtable.h"
#include "linkedlist.h"
#include "trace.h"
#include "xwin.h"

#define NIL (0)

void RenderWindowCallback(int x, int y, vector_t color, void *data);

typedef struct
{
	Display *display;
	Colormap cmap;
	GC gc;
	Window window;
	int level;
	int pix_width, pix_height;
	int width, height;
	int last_x, last_y;
	int last_mouse_x, last_mouse_y;
	render_settings_t *rs;
	kd_tree_t *tree_ptr;
	vecc_t rot_x, rot_y;
	vector_t cam_position;
	Atom wmDeleteMessage;
} callback_data_t;

static matrix_t GetTransform(callback_data_t *cd)
{
	return
		MatrixTimesMatrix(
			NewTranslateMatrix(cd->cam_position),
			MatrixTimesMatrix(
				NewRotateMatrix(cd->rot_x, vec3(0, -1, 0)),
				NewRotateMatrix(cd->rot_y, vec3(-1, 0, 0))));
}

static void RenderWindow(callback_data_t *cd)
{
	int width = cd->width;
	int height = cd->height;
	
	vecc_t aspect = (vecc_t)height / (vecc_t)width;
	vecc_t scale_x = 2 / (vecc_t)width;
	vecc_t scale_y = -aspect * 2 / (vecc_t)height;
	
	cd->level = 0;
	cd->last_x = 0;
	cd->last_y = 0;
	
	int size_pow_2;
	if(width > height)
		for(size_pow_2 = 1; size_pow_2 < width; size_pow_2 *= 2);
	else
		for(size_pow_2 = 1; size_pow_2 < height; size_pow_2 *= 2);
	
	cd->pix_width = size_pow_2 / 2;
	cd->pix_height = size_pow_2 / 2;
	
	Render(
		RenderWindowCallback,
		(void *)cd,
		width, height,
		GetTransform(cd),
		cd->tree_ptr,
		cd->rs->lights_ptr,
		cd->rs->scene_ptr,
		10000,	// Override section size.
		aspect,
		scale_x, scale_y,
		cd->rs->focal_length,
		cd->rs->samples,
		cd->rs->max_iterations,
		cd->rs->threads);
}

void RenderWindowCallback(int x, int y, vector_t color, void *data)
{
	
	callback_data_t *d = (callback_data_t *)data;
	
	// Cancel rendering if there are any events to handle.
	if(XPending(d->display))
	{
		CancelRender();
		return;
	}
	/*XColor xcolor;
	xcolor.red = 0xFFFF * color.x;
	xcolor.green = 0xFFFF * color.y;
	xcolor.blue = 0xFFFF * color.z;
	xcolor.flags = DoRed | DoGreen | DoBlue;
	XAllocColor(d->display, d->cmap, &xcolor);*/
	
	// When a level is finished.
	if(x < d->last_x)
	{
		d->level += 1;
		d->pix_width /= 2;
		d->pix_height /= 2;
		XFlush(d->display);
	}
	
	d->last_x = x;
	d->last_y = y;
	
	unsigned int c =
		((int) (color.x * 0xFF) << 16) |
		((int) (color.y * 0xFF) << 8) |
		(int) (color.z * 0xFF);
	
	
	XSetForeground(d->display, d->gc, c);
	if(d->pix_width == 1 && d->pix_height == 1)
	{
		XDrawPoint(d->display, d->window, d->gc, x, y);
	}
	else
		XFillRectangle(d->display, d->window, d->gc, x, y, d->pix_width, d->pix_height);
}

static void *EventLoop(void *data)
{
	callback_data_t *cd = (callback_data_t *)data;
	XEvent e;
	while(1)
	{
		XNextEvent(cd->display, &e);
		if(e.type == MapNotify)
		{
			
		}
		else if(e.type == MotionNotify)
		{
			XMotionEvent me = e.xmotion;
			vecc_t dx, dy;
			dx = me.x - cd->last_mouse_x;
			dy = me.y - cd->last_mouse_y;
			
			cd->rot_x += dx / 100;
			cd->rot_y += dy / 100;
			
			cd->last_mouse_x = me.x;
			cd->last_mouse_y = me.y;
			RenderWindow(cd);
		}
		else if(e.type == ButtonPress)
		{
			XButtonEvent be = e.xbutton;
			cd->last_mouse_x = be.x;
			cd->last_mouse_y = be.y;
		}
		else if(e.type == KeyPress)
		{
			XKeyEvent ke = e.xkey;
			KeySym ks = XLookupKeysym(&ke, 0);
			vector_t cam_velocity = vec4(0, 0, 0, 0);
			if(ks == XK_a)
				cam_velocity.x = -1;
			else if(ks == XK_d)
				cam_velocity.x = 1;
			else if(ks == XK_w)
				cam_velocity.z = -1;
			else if(ks == XK_s)
				cam_velocity.z = 1;
			cd->cam_position = 
				VectorPlusVector(
					cd->cam_position,
					MatrixTimesVector(
						GetTransform(cd),
						cam_velocity));
			RenderWindow(cd);
		}
		else if(e.type == Expose)
		{
			int width = e.xexpose.width;
			int height = e.xexpose.height;
			
			cd->width = width;
			cd->height = height;
			RenderWindow(cd);
		}
		else if(e.xclient.data.l[0] == cd->wmDeleteMessage)
		{
			break;
		}
	}
	
	FreeTree(cd->tree_ptr);
	XDestroyWindow(cd->display, cd->window);
	XCloseDisplay(cd->display);
	free(cd);
	return NULL;
}

int WindowMain(render_settings_t *rs, int width, int height)
{
	XInitThreads();
	Display *display = XOpenDisplay(NULL);
	
	int blackColor = BlackPixel(display, DefaultScreen(display));
	
	Window window = XCreateSimpleWindow(
		display, DefaultRootWindow(display), 0, 0, width, height,
		0, blackColor, blackColor);
	XStoreName(display, window, "Tom's Raytracer");
	XSelectInput(display, window,
		StructureNotifyMask | ExposureMask
		| ButtonMotionMask | ButtonPressMask
		| KeyPressMask);
	
	XMapWindow(display, window);
	
	Atom wmDeleteMessage = XInternAtom(display, "WM_DELETE_WINDOW", False);
	XSetWMProtocols(display, window, &wmDeleteMessage, 1);
	
	Colormap cmap = DefaultColormap(display, DefaultScreen(display));
	
	printf("info: generating tree\n");
	CleanPrimitives(rs, rs->primitives_ptr);
	kd_tree_t *tree_ptr = GenerateTree(rs->primitives_ptr, 0, NULL);
	
	callback_data_t *cd = malloc(sizeof(callback_data_t));
	cd->display = display;
	cd->cmap = cmap;
	
	cd->window = window;
	cd->cam_position = vec3(0, 0, 0);
	cd->rot_x = 0;
	cd->rot_y = 0;
	
	cd->tree_ptr = tree_ptr;
	cd->wmDeleteMessage = wmDeleteMessage;
	cd->rs = rs;
	
	cd->gc = XCreateGC(display, window, 0, NULL);
	EventLoop(cd);
	
	return 1;
}


#endif

