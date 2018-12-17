#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "hashtable.h"
#include "linkedlist.h"
#include "trace.h"
#include "lodepng.h"
#include "shade.h"
#include "cmd.h"
#include "texture.h"

void ImageFree(void *ptr)
{
	image_t *image_ptr = (image_t *)ptr;
	free(image_ptr->data);
	free(image_ptr);
}

void Texture2dFree(void *ptr)
{
	free(ptr);
}

vector_t Texture2dSampleNearestNeighbor(texture2d_t *tex_ptr, vector_t co)
{
	image_t *image_ptr = tex_ptr->image_ptr;
	int indexX = (int)(co.x * image_ptr->width) % image_ptr->width;
	int indexY = (int)(co.y * image_ptr->height) % image_ptr->height;
	
	if(indexX < 0)
		indexX += image_ptr->width;
	if(indexY < 0)
		indexY += image_ptr->height;
	
	unsigned int pixelIndex = (indexX + indexY * image_ptr->width) * 4;
	
	return vec4(
		((vecc_t)image_ptr->data[pixelIndex + 0]) / 255,
		((vecc_t)image_ptr->data[pixelIndex + 1]) / 255,
		((vecc_t)image_ptr->data[pixelIndex + 2]) / 255,
		((vecc_t)image_ptr->data[pixelIndex + 3]) / 255);
}

vector_t Texture2dSampleBilinear(texture2d_t *tex_ptr, vector_t co)
{
	image_t *image_ptr = tex_ptr->image_ptr;
	
	vecc_t offX = 1.0 / image_ptr->width;
	vecc_t offY = 1.0 / image_ptr->height;
	
	vecc_t x1 = (vecc_t)((int)(co.x * image_ptr->width)) / (vecc_t)image_ptr->width;
	vecc_t y1 = (vecc_t)((int)(co.y * image_ptr->height)) / (vecc_t)image_ptr->height;
	vecc_t x2 = x1 + offX;
	vecc_t y2 = y1 + offY;
	vecc_t x = co.x;
	vecc_t y = co.y;
	
	matrix_t mat =
		(MatrixInverse(
			NewMatrix(
				1, x1, y1, x1 * y1,
				1, x1, y2, x1 * y2,
				1, x2, y1, x2 * y1,
				1, x2, y2, x2 * y2)));
	vector_t Q11 = Texture2dSampleNearestNeighbor(tex_ptr, vec2(x1, y1));
	vector_t Q12 = Texture2dSampleNearestNeighbor(tex_ptr, vec2(x1, y2));
	vector_t Q21 = Texture2dSampleNearestNeighbor(tex_ptr, vec2(x2, y1));
	vector_t Q22 = Texture2dSampleNearestNeighbor(tex_ptr, vec2(x2, y2));
	
	vector_t r = MatrixTimesVector(mat, vec4(Q11.x, Q12.x, Q21.x, Q22.x));
	vector_t g = MatrixTimesVector(mat, vec4(Q11.y, Q12.y, Q21.y, Q22.y));
	vector_t b = MatrixTimesVector(mat, vec4(Q11.z, Q12.z, Q21.z, Q22.z));
	vector_t a = MatrixTimesVector(mat, vec4(Q11.w, Q12.w, Q21.w, Q22.w));
	return vec4(
		r.x + r.y * x + r.z * y + r.w * x * y,
		g.x + g.y * x + g.z * y + g.w * x * y,
		b.x + b.y * x + b.z * y + b.w * x * y,
		a.x + a.y * x + a.z * y + a.w * x * y);
}

vector_t Texture2dSample(texture2d_t *tex_ptr, vector_t co)
{
	if(tex_ptr->filter == BILINEAR)
		return Texture2dSampleBilinear(tex_ptr, co);
	else
		return Texture2dSampleNearestNeighbor(tex_ptr, co);
}

// CMDS

/*
 * Makes a texture from the given image index.
 */
int CmdMakeTexture2d(render_settings_t *rs)
{
	int image_index;
	if(fscanf(rs->input, "%i", &image_index) == 0)
	{
		fprintf(stderr, "error: 'make_texture_2d' not enough arguments.\n");
		return 1;
	}
	texture2d_t *tex_ptr = (texture2d_t *)malloc(sizeof(texture2d_t));
	memset(tex_ptr, 0, sizeof(texture2d_t));
	tex_ptr->image_ptr = (image_t *)ListGet(rs->images_ptr, image_index).data_ptr;
	tex_ptr->filter = BILINEAR;
	ListAppend(rs->texture2ds_ptr, (list_data_t)((void *)tex_ptr));
	rs->texture2d_ptr = tex_ptr;
	return 1;
}

/*
 * Makes an image with the specified width and height.
 */
int CmdMakeImage(render_settings_t *rs)
{
	int count, width, height;
	count = fscanf(rs->input, "%i %i", &width, &height);
	if(count == 0)
	{
		fprintf(stderr, "error: 'make_image' not enough arguments.\n");
		return 1;
	}
	image_t *image_ptr = (image_t *)malloc(sizeof(image_t));
	image_ptr->data = malloc(4 * image_ptr->width * image_ptr->height);
	image_ptr->width = width;
	image_ptr->height = height;
	memset(&image_ptr->data, 0, 4 * image_ptr->width * image_ptr->height);
	
	ListAppend(rs->images_ptr, (list_data_t)((void *)image_ptr));
	rs->image_ptr = image_ptr;
	
	return 1;
}

/*
 * Loads an image from the specified file.
 */
int CmdLoadImage(render_settings_t *rs)
{
	int count, error;
	char filename[64];
	count = fscanf(rs->input, "%s", filename);
	
	if(count == 0)
	{
		fprintf(stderr, "error: 'load_image' not enough arguments.\n");
		return 1;
	}
	
	unsigned int width, height;
	unsigned char *data;
	
	error = lodepng_decode32_file(&data, &width, &height, filename);
	if(error)
	{
		fprintf(stderr, "lodepng: error %u: %s\n", error, lodepng_error_text(error));
		return 1;
	}
	
	image_t *image_ptr = (image_t *)malloc(sizeof(image_t));
	image_ptr->width = width;
	image_ptr->height = height;
	image_ptr->data = data;
	
	ListAppend(rs->images_ptr, (list_data_t)((void *)image_ptr));
	rs->image_ptr = image_ptr;
	return 1;
}

void CommandsSetTexture(hashtable_t *table_cmds)
{
	SetCommand(table_cmds, "make_texture_2d", CmdMakeTexture2d);
	SetCommand(table_cmds, "make_image", CmdMakeImage);
	SetCommand(table_cmds, "load_image", CmdLoadImage);
}



