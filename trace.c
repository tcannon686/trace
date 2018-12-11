#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _WIN32
#	include <Windows.h>
#else
#	include <pthread.h>
#endif
#include <time.h>

#ifndef min
#define min(a, b)			a < b ? a : b
#endif

#ifndef max
#define max(a, b)			a > b ? a : b
#endif

#include "matrix.h"
#include "hashtable.h"
#include "trace.h"
#include "lodepng.h"
#include "shade.h"
#include "cmd.h"

#ifdef INCLUDE_GUI
#include "xwin.h"
#endif

#define DEFAULT_OUTPUT "output.png"
#define THRESHOLD 0.00001

vector_t RandomVector(int axes)
{
	vector_t ret = vec2(0, 0);
	for(int i = 0; i < axes; i ++)
	{
		ret.m[i] = 2 * ((vecc_t)rand() / RAND_MAX) - 1;
	}
	return ret;
}

vector_t Barycentric(vector_t v, triangle_t tri)
{
	vector_t ret;
	vecc_t area_inv = 1.0 / tri.area;
	
	vector_t cross;
	
	cross = VectorCrossVector(
			VectorMinusVector(v, tri.a),
			tri.ba);
	
	if(VectorDotVector(tri.normal, cross) < 0)
	{
		ret.x = -1;
		return ret;
	}
	else
		ret.z = VectorMagnitude(cross) * area_inv;
	
	cross = VectorCrossVector(
			VectorMinusVector(v, tri.b),
			tri.cb);
	
	if(VectorDotVector(tri.normal, cross) < 0)
	{
		ret.x = -1;
		return ret;
	}
	else
		ret.x = VectorMagnitude(cross) * area_inv;
	
	cross = VectorCrossVector(
			VectorMinusVector(v, tri.c),
			tri.ac);
	
	if(VectorDotVector(tri.normal, cross) < 0)
	{
		ret.x = -1;
		return ret;
	}
	else
		ret.y = VectorMagnitude(cross) * area_inv;
	
	ret.w = 0;
	return ret;
}

int RayTriangle(hit_t *hit_ptr, ray_t *ray, triangle_t tri)
{
	vecc_t ldotn = VectorDotVector(ray->d, tri.normal);
	if(ldotn == 0)
		return 0;
	vecc_t d = VectorDotVector(VectorMinusVector(tri.v0, ray->o), tri.normal);
	
	if(d * ldotn <= 0)
		return 0;
	d /= ldotn;
	if(d > hit_ptr->t)
		return 0;
	
	
	
	vector_t pos = VectorPlusVector(ray->o, VectorTimesScalar(ray->d, d));
	vector_t bary = Barycentric(pos, tri);
	
	if(bary.x < 0
		|| bary.y < 0
		|| bary.z < 0)
	{
		return 0;
	}
	
	
	hit_ptr->t = d;
	hit_ptr->normal = tri.normal;
	hit_ptr->material_ptr = tri.material_ptr;
	hit_ptr->position = pos;
	hit_ptr->ray = *ray;
	
	hit_ptr->normal =
		VectorNormalize(VectorPlusVector(
			VectorTimesScalar(tri.n0, bary.x),
			VectorPlusVector(
				VectorTimesScalar(tri.n1, bary.y),
				VectorTimesScalar(tri.n2, bary.z))));
	
	return 1;
}

int RayTriangles(hit_t *hit_ptr, ray_t *ray, tri_list_t *triangles_ptr)
{
	tri_list_t *current_ptr = triangles_ptr;
	int ret = 0;
	while(current_ptr != NULL)
	{
		ret = ret || RayTriangle(hit_ptr, ray, current_ptr->current);
		current_ptr = current_ptr->next_ptr;
	}
	return ret;
}

bounding_box_t BoundingBox(triangle_t tri)
{
	bounding_box_t ret;
	ret.a = tri.v[0];
	ret.b = tri.v[1];
	for(int i = 0; i < 3; i ++)
	{
		for(int j = 0; j < 3; j ++)
		{
			if(tri.v[i].m[j] < ret.a.m[j])
				ret.a.m[j] = tri.v[i].m[j];
			if(tri.v[i].m[j] > ret.b.m[j])
				ret.b.m[j] = tri.v[i].m[j];
		}
	}
	
	ret.origin = VectorMean(ret.a, ret.b);
	return ret;
}

bounding_box_t BoundingBoxTriList(tri_list_t *triangles_ptr)
{
	bounding_box_t ret;
	ret.a = triangles_ptr->current.v[0];
	ret.b = triangles_ptr->current.v[1];
	
	while(triangles_ptr != NULL)
	{
		triangle_t tri = triangles_ptr->current;
		for(int i = 0; i < 3; i ++)
		{
			for(int j = 0; j < 3; j ++)
			{
				if(tri.v[i].m[j] < ret.a.m[j])
					ret.a.m[j] = tri.v[i].m[j];
				if(tri.v[i].m[j] > ret.b.m[j])
					ret.b.m[j] = tri.v[i].m[j];
			}
		}
		triangles_ptr = triangles_ptr->next_ptr;
	}
	
	ret.origin = VectorMean(ret.a, ret.b);
	return ret;
}

bounding_box_t BoundingBoxTriListOrigins(tri_list_t *triangles_ptr)
{
	bounding_box_t ret;
	ret.a = triangles_ptr->current.origin;
	ret.b = triangles_ptr->current.origin;
	
	while(triangles_ptr != NULL)
	{
		triangle_t tri = triangles_ptr->current;
		for(int j = 0; j < 3; j ++)
		{
			if(tri.origin.m[j] < ret.a.m[j])
				ret.a.m[j] = tri.origin.m[j];
			if(tri.origin.m[j] > ret.b.m[j])
				ret.b.m[j] = tri.origin.m[j];
		}
		triangles_ptr = triangles_ptr->next_ptr;
	}
	
	ret.origin = VectorMean(ret.a, ret.b);
	return ret;
}


bounding_box_t Expand(bounding_box_t box)
{
	box.a.x -= THRESHOLD;
	box.a.y -= THRESHOLD;
	box.a.z -= THRESHOLD;
	
	box.b.x += THRESHOLD;
	box.b.y += THRESHOLD;
	box.b.z += THRESHOLD;
	return box;
}

void PrintVector(vector_t v)
{
	const char *string = VectorToString(v);
	printf("%s", string);
	free((void *)string);
}

vector_t VectorMix(vector_t left, vector_t right, vecc_t a)
{
	return VectorPlusVector(VectorTimesScalar(left, 1 - a), VectorTimesScalar(right, a));
}


// Remove triangles with area of 0.
void CleanTriangles(tri_list_t **triangles_ptr)
{
	tri_list_t *current = *triangles_ptr;
	while(current != NULL)
	{
		tri_list_t *next_ptr = current->next_ptr;
		if(current->current.area == 0)
		{
			if(current->last_ptr != NULL)
				current->last_ptr->next_ptr = next_ptr;
			else
				*triangles_ptr = next_ptr;
			if(next_ptr != NULL)
				next_ptr->last_ptr = current->last_ptr;
		}
		current = next_ptr;
	}
}

/*
 * Generates a tree and stores the maximum depth in depth_result.
 */
kd_tree_t *GenerateTree(tri_list_t *triangles_ptr, int depth, int *depth_result) {

	// Return nothing if no triangles are provided
	if(triangles_ptr == NULL)
		return NULL;
	
	if(depth_result != NULL)
	{
		if(depth >= *depth_result)
			*depth_result = depth;
	}
	
	kd_tree_t *tree_ptr = (kd_tree_t *)malloc(sizeof(kd_tree_t));
	memset(tree_ptr, 0, sizeof(kd_tree_t));
	tree_ptr->box = Expand(BoundingBoxTriList(triangles_ptr));
	tree_ptr->triangles_ptr = triangles_ptr;

	bounding_box_t box_origins = BoundingBoxTriListOrigins(triangles_ptr);

	vector_t origin = vec3(0, 0, 0);
	vector_t last = triangles_ptr->current.origin;
	tri_list_t *current_ptr;
	int count = 0, same = 1;

	for(current_ptr = triangles_ptr; current_ptr != NULL; current_ptr = current_ptr->next_ptr)
	{
		origin = VectorPlusVector(origin, current_ptr->current.origin);
		if(same && !VectorEqualsVector(last, current_ptr->current.origin))
			same = 0;
		count ++;
	}

	origin = VectorTimesScalar(origin, 1.0 / count);

	vector_t dimensions = VectorMinusVector(box_origins.b, box_origins.a);
	if(dimensions.x > dimensions.y && dimensions.x > dimensions.z)
		tree_ptr->axis_index = 0;
	else if(dimensions.y > dimensions.x && dimensions.y > dimensions.z)
		tree_ptr->axis_index = 1;
	else if(dimensions.z > dimensions.x && dimensions.z > dimensions.y)
		tree_ptr->axis_index = 2;
	else
		tree_ptr->axis_index = depth % 3;
	
	tree_ptr->axis = vec4(0, 0, 0, 0);
	tree_ptr->axis.m[tree_ptr->axis_index] = 1;
	tree_ptr->median = origin;

	// Break if this is a leaf node (i.e. all the triangles have the same origin)
	if(same == 1)
		return tree_ptr;
	
	// Otherwise split the triangles based on left or right

	// There should be no triangles at this node
	tree_ptr->triangles_ptr = NULL;

	tri_list_t
		*left_triangles_ptr = NULL,
		*right_triangles_ptr = NULL,
		*left_ptr = NULL,
		*right_ptr = NULL;
	tri_list_t *next_ptr;

	for(
		current_ptr = triangles_ptr; current_ptr != NULL; current_ptr = next_ptr)
	{
		// Make next_ptr be the next triangle, so that iteration continues
		// even if triangles are added or removed from the list
		next_ptr = current_ptr->next_ptr;

		vector_t direction = VectorMinusVector(current_ptr->current.origin, origin);
		vecc_t dot = VectorDotVector(direction, tree_ptr->axis);

		// Left side
		if(dot <= 0)
		{
			// If there are no triangles yet
			if(left_triangles_ptr == NULL)
			{
				// Remove the triangle from the list
				if(current_ptr->last_ptr != NULL)
					current_ptr->last_ptr->next_ptr = current_ptr->next_ptr;
				if(current_ptr->next_ptr != NULL)
					current_ptr->next_ptr->last_ptr = current_ptr->last_ptr;

				// Move the triangle to the left list
				left_triangles_ptr = current_ptr;
				left_ptr = current_ptr;
				current_ptr->last_ptr = NULL;
				current_ptr->next_ptr = NULL;
			}
			// If there are triangles
			else
			{
				// Remove the triangle from the list
				if(current_ptr->last_ptr != NULL)
					current_ptr->last_ptr->next_ptr = current_ptr->next_ptr;
				if(current_ptr->next_ptr != NULL)
					current_ptr->next_ptr->last_ptr = current_ptr->last_ptr;

				// Move the triangle to the left list
				left_ptr->next_ptr = current_ptr;
				current_ptr->last_ptr = left_ptr;
				current_ptr->next_ptr = NULL;
				left_ptr = current_ptr;
			}
		}
		// Right side
		else
		{
			// If there are no triangles yet
			if(right_triangles_ptr == NULL)
			{
				// Remove the triangle from the list
				if(current_ptr->last_ptr != NULL)
					current_ptr->last_ptr->next_ptr = current_ptr->next_ptr;
				if(current_ptr->next_ptr != NULL)
					current_ptr->next_ptr->last_ptr = current_ptr->last_ptr;

				// Move the triangle to the right list
				right_triangles_ptr = current_ptr;
				right_ptr = current_ptr;
				current_ptr->last_ptr = NULL;
				current_ptr->next_ptr = NULL;
			}
			// If there are triangles
			else
			{
				// Remove the triangle from the list
				if(current_ptr->last_ptr != NULL)
					current_ptr->last_ptr->next_ptr = current_ptr->next_ptr;
				if(current_ptr->next_ptr != NULL)
					current_ptr->next_ptr->last_ptr = current_ptr->last_ptr;

				// Move the triangle to the right list
				right_ptr->next_ptr = current_ptr;
				current_ptr->last_ptr = right_ptr;
				current_ptr->next_ptr = NULL;
				right_ptr = current_ptr;
			}
		}
	}

	tree_ptr->left_ptr = GenerateTree(left_triangles_ptr, depth + 1, depth_result);
	if(tree_ptr->left_ptr != NULL)
		tree_ptr->left_ptr->parent_ptr = tree_ptr;
	tree_ptr->right_ptr = GenerateTree(right_triangles_ptr, depth + 1, depth_result);
	if(tree_ptr->right_ptr != NULL)
		tree_ptr->right_ptr->parent_ptr = tree_ptr;
	
	return tree_ptr;
}

void FreeTriList(tri_list_t *list_ptr) {
	tri_list_t *current_ptr;
	tri_list_t *next_ptr;
	current_ptr = list_ptr;
	while(current_ptr != NULL) {
		next_ptr = current_ptr->next_ptr;
		free(current_ptr);
		current_ptr = next_ptr;
	}
}

void FreeTree(kd_tree_t *tree_ptr) {
	if(tree_ptr == NULL)
		return;
	if(tree_ptr->triangles_ptr != NULL)
		FreeTriList(tree_ptr->triangles_ptr);
	if(tree_ptr->left_ptr != NULL)
		FreeTree(tree_ptr->left_ptr);
	if(tree_ptr->right_ptr != NULL)
		FreeTree(tree_ptr->right_ptr);
}

void FreeMaterialList(mat_list_t *list_ptr) {
	mat_list_t *next_ptr;
	while(list_ptr != NULL) {
		next_ptr = list_ptr->next_ptr;
		HashTableFree(list_ptr->current.table);
		free(list_ptr);
		list_ptr = next_ptr;
	}
}

void FreeLightList(light_list_t *list_ptr) {
	light_list_t *next_ptr;
	while(list_ptr != NULL) {
		next_ptr = list_ptr->next_ptr;
		HashTableFree(list_ptr->current.table);
		free(list_ptr);
		list_ptr = next_ptr;
	}
}

vecc_t RayBox(ray_t *ray, bounding_box_t box)
{
	vecc_t tmin, tmax, tymin, tymax, tzmin, tzmax;
	if(ray->d_inverse.x >= 0)
	{
		tmin = (box.a.x - ray->o.x) * ray->d_inverse.x;
		tmax = (box.b.x - ray->o.x) * ray->d_inverse.x;
	}
	else
	{
		tmin = (box.b.x - ray->o.x) * ray->d_inverse.x;
		tmax = (box.a.x - ray->o.x) * ray->d_inverse.x;
	}
	
	if(ray->d_inverse.y >= 0)
	{
		tymin = (box.a.y - ray->o.y) * ray->d_inverse.y;
		tymax = (box.b.y - ray->o.y) * ray->d_inverse.y;
	}
	else
	{
		tymin = (box.b.y - ray->o.y) * ray->d_inverse.y;
		tymax = (box.a.y - ray->o.y) * ray->d_inverse.y;
	}
	
	if(tmin > tymax || tymin > tmax)
		return INFINITY;
	
	if(tymin > tmin)
		tmin = tymin;
	if(tymax < tmax)
		tmax = tymax;
	
	
	if(ray->d_inverse.z >= 0)
	{
		tzmin = (box.a.z - ray->o.z) * ray->d_inverse.z;
		tzmax = (box.b.z - ray->o.z) * ray->d_inverse.z;
	}
	else
	{
		tzmin = (box.b.z - ray->o.z) * ray->d_inverse.z;
		tzmax = (box.a.z - ray->o.z) * ray->d_inverse.z;
	}
	if(tmin > tzmax || tzmin > tmax)
		return INFINITY;
	
	
	if(tzmin > tmin)
		tmin = tzmin;
	if(tzmax < tmax)
		tmax = tzmax;

	if(tmax > 0)
		return tmin;
	else
		return INFINITY;
}

ray_t NewRay(vector_t o, vector_t d)
{
	ray_t ret;
	ret.o = o;
	ret.d = d;
	ret.d_inverse.x = 1.0 / ret.d.x;
	ret.d_inverse.y = 1.0 / ret.d.y;
	ret.d_inverse.z = 1.0 / ret.d.z;
	return ret;
}


int RayTreeRecursive(hit_t *hit_ptr, ray_t *ray, kd_tree_t *tree_ptr)
{
	if(tree_ptr == NULL)
		return 0;
	vecc_t box_hit = RayBox(ray, tree_ptr->box);
	if(box_hit < INFINITY && (box_hit < hit_ptr->t || hit_ptr->t >= INFINITY))
	{
		if(tree_ptr->triangles_ptr != NULL)
		{
			int hit = RayTriangles(hit_ptr, ray, tree_ptr->triangles_ptr);
			if(hit)
				hit_ptr->node_ptr = tree_ptr;
			return hit;
		}
		
		int left = 0, right = 0;
		if(ray->o.m[tree_ptr->axis_index] > tree_ptr->median.m[tree_ptr->axis_index])
		{
			right = RayTreeRecursive(hit_ptr, ray, tree_ptr->right_ptr);
			left = RayTreeRecursive(hit_ptr, ray, tree_ptr->left_ptr);
		}
		else
		{
			left = RayTreeRecursive(hit_ptr, ray, tree_ptr->left_ptr);
			right = RayTreeRecursive(hit_ptr, ray, tree_ptr->right_ptr);
		}
		return left || right;
		
	}
	return 0;
}

int RayTree(hit_t *hit_ptr, ray_t ray, kd_tree_t *tree_ptr)
{
	return RayTreeRecursive(hit_ptr, &ray, tree_ptr);
}

void PrintTree(kd_tree_t *tree_ptr, int indent)
{
	const char *boxA, *boxB, *median, *vertex;
	boxA = VectorToString(tree_ptr->box.a);
	boxB = VectorToString(tree_ptr->box.b);
	median = VectorToString(tree_ptr->median);
	
	for(int i = 0; i < indent; i++)
		printf(" ");
	if(tree_ptr->triangles_ptr != NULL)
		printf("leaf node ");
	printf("{\n");
	for(int i = 0; i < indent; i++)
		printf(" ");
	printf("    Box A: %s\n", boxA);
	for (int i = 0; i < indent; i++)
		printf(" ");
	printf("    Box B: %s\n", boxB);
	for (int i = 0; i < indent; i++)
		printf(" ");
	printf("    Median: %s\n", median);
	for (tri_list_t *current = tree_ptr->triangles_ptr; current != NULL; current = current->next_ptr) {
		for (int i = 0; i < indent; i++)
			printf(" ");
		printf("    Triangle %p:\n", current);
		for (int i = 0; i < 3; i++) {
			vertex = VectorToString(current->current.v[i]);
			for (int j = 0; j < indent; j++)
				printf(" ");
			printf("        Vertex: %s\n", median);
			free((void *)vertex);
		}
	}
	free((void *)boxA);
	free((void *)boxB);
	free((void *)median);
	
	if(tree_ptr->left_ptr != NULL)
		PrintTree(tree_ptr->left_ptr, indent + 4);
	if(tree_ptr->right_ptr != NULL)
		PrintTree(tree_ptr->right_ptr, indent + 4);
	
	for(int i = 0; i < indent; i++)
		printf(" ");
	printf("}\n");
}

void PrintTriangles(tri_list_t *triangles_ptr) {
	tri_list_t *current_ptr = triangles_ptr;
	while(current_ptr != NULL)
	{
		const char *v0, *v1, *v2, *n, *o;
		v0 = VectorToString(current_ptr->current.v0);
		v1 = VectorToString(current_ptr->current.v1);
		v2 = VectorToString(current_ptr->current.v2);
		n = VectorToString(current_ptr->current.normal);
		o = VectorToString(current_ptr->current.origin);
		printf("triangle: {\n    %s\n    %s\n    %s\n\n    normal: %s\n    origin: %s\n}\n",
			v0, v1, v2, n, o);
		free((void *)v0);
		free((void *)v1);
		free((void *)v2);
		free((void *)n);
		free((void *)o);
		current_ptr = current_ptr->next_ptr;
	}
}


void RenderSection(render_params_t *rp_ptr)
{
	int sampleWidth = (int)sqrt(rp_ptr->samples);
	vecc_t sampleSize = 1.0 / sampleWidth;
	kd_tree_t *last_hit = NULL;
	vector_t sample_offset;
	hit_t result;
	vector_t sample;
	for(int i = rp_ptr->x0; i < rp_ptr->x1; i ++)
	{
		for(int j = rp_ptr->y0; j < rp_ptr->y1; j ++)
		{
			vector_t pixel = vec3(0, 0, 0);
			for(int k = 0; k < rp_ptr->samples; k ++)
			{
				sample_offset = 
					vec2((k % sampleWidth) * sampleSize, (k / sampleWidth) * sampleSize);
				ray_t ray = NewRay(
					vec3(0, 0, 0),
					vec3(rp_ptr->scale_x * (((vecc_t)i - rp_ptr->width / 2.0) + sample_offset.x),
						rp_ptr->scale_y * (((vecc_t)j - rp_ptr->height / 2.0) + sample_offset.y),
						-rp_ptr->focal_length));
				result.t = INFINITY;
				result.tree_ptr = rp_ptr->tree_ptr;
				sample = rp_ptr->sky_color;
				if(last_hit != NULL)
					RayTree(&result, ray, last_hit);
				if(RayTree(&result, ray, rp_ptr->tree_ptr))
				{
					sample = ShadeHit(result, rp_ptr, 0);
					last_hit = result.node_ptr;
				}
				VectorPlusVectorP(&pixel, &pixel, &sample);
			}
			VectorTimesScalarP(&pixel, &pixel, 1.0 / rp_ptr->samples);
			rp_ptr->callback(i, j, pixel, rp_ptr->callback_data);
		}
	}
}

#ifdef _WIN32
DWORD WINAPI RenderThread(LPVOID lpParameter)
{
	render_params_list_t *rps = (render_params_list_t *)lpParameter;
	for(render_params_list_t *cur = rps; cur != NULL; cur = cur->next_ptr)
	{
		if(cur->is_rendering)
			continue;
		cur->is_rendering = 1;
		RenderSection(&cur->current);
		printf("info: rendered section (%i, %i), (%i, %i).\n",
			cur->current.x0, cur->current.y0, cur->current.x1, cur->current.y1);
	}
	return 0;
}
#else
void *RenderThread(void *lpParameter)
{
	render_params_list_t *rps = (render_params_list_t *)lpParameter;
	for(render_params_list_t *cur = rps; cur != NULL; cur = cur->next_ptr)
	{
		if(cur->is_rendering)
			continue;
		cur->is_rendering = 1;
		RenderSection(&cur->current);
		printf("info: rendered section (%i, %i), (%i, %i).\n",
			cur->current.x0, cur->current.y0, cur->current.x1, cur->current.y1);
	}
	return NULL;
}
#endif


void Render(
	pixel_traced_callback_t callback,
	void *callback_data,
	int width, int height,
	kd_tree_t *tree_ptr,
	light_list_t *lights_ptr,
	vector_t sky_color,
	int section_size,
	vecc_t aspect,
	vecc_t scale_x, vecc_t scale_y,
	vecc_t focal_length,
	int samples,
	int shadow_samples,
	int max_iterations,
	int num_threads)
{
#ifdef _WIN32
	HANDLE threads[8];
#else
	pthread_t threads[8];
#endif
	render_params_list_t *rps = NULL;
	render_params_list_t *last_ptr = NULL;
	for(int i = 0; i < width; i += section_size)
	{
		for(int j = 0; j < height; j += section_size)
		{
			render_params_list_t *rp = (render_params_list_t *)malloc(sizeof(render_params_list_t));
			rp->is_rendering = 0;
			rp->next_ptr = NULL;
			rp->current.width = width;
			rp->current.height = height;
			rp->current.tree_ptr = tree_ptr;
			rp->current.lights_ptr = lights_ptr;
			rp->current.sky_color = sky_color;
			rp->current.x0 = i;
			rp->current.y0 = j;
			rp->current.x1 = min(i + section_size, width);
			rp->current.y1 = min(j + section_size, height);
			rp->current.aspect = aspect;
			rp->current.scale_x = scale_x;
			rp->current.scale_y = scale_y;
			rp->current.focal_length = focal_length;
			rp->current.samples = samples;
			rp->current.shadow_samples = shadow_samples;
			rp->current.max_iterations = max_iterations;
			rp->current.callback = callback;
			rp->current.callback_data = callback_data;
			if(last_ptr != NULL)
				last_ptr->next_ptr = rp;
			else
				rps = rp;
			last_ptr = rp;
		}
	}
	
	if(num_threads > 1)
	{
	
		for(int i = 0; i < num_threads; i++)
		{
#ifdef _WIN32
			threads[i] = CreateThread(
				NULL, 0, RenderThread, (LPVOID)rps, 0, NULL);
			
			if(threads[i] == NULL)
			{
				fprintf(stderr, "error: failed to create thread.\n");
				return;
			}
#else
			if(pthread_create(&(threads[i]), NULL, &RenderThread, (void *)rps)
				!= 0)
			{
				fprintf(stderr, "error: failed to create thread.\n");
				return;
			}
#endif
		}

#ifdef _WIN32
		WaitForMultipleObjects(num_threads, threads, TRUE, INFINITE);
#else
		for(int i = 0; i < num_threads; i++)
		{
			pthread_join(threads[i], NULL);
		}
#endif
	
		last_ptr = rps;
		while(last_ptr != NULL)
		{
			render_params_list_t *next_ptr = last_ptr->next_ptr;
			free(last_ptr);
			last_ptr = next_ptr;
		}
	
	}
	else
	{
		for(int i = 0; i < num_threads; i++)
		{
#ifdef _WIN32
			RenderThread((LPVOID)rps);
#else
			RenderThread((void *)rps);
#endif
		}
		
		last_ptr = rps;
		while(last_ptr != NULL)
		{
			render_params_list_t *next_ptr = last_ptr->next_ptr;
			free(last_ptr);
			last_ptr = next_ptr;
		}
	}
	
}

void RenderImageCallback(int x, int y, vector_t color, void *data)
{
	image_t *image_ptr = (image_t *)data;
	int pixel_index = (x + y * image_ptr->width) * 4;
	image_ptr->data[pixel_index + 0] = (unsigned char)(color.m[0] * 255);
	image_ptr->data[pixel_index + 1] = (unsigned char)(color.m[1] * 255);
	image_ptr->data[pixel_index + 2] = (unsigned char)(color.m[2] * 255);
	image_ptr->data[pixel_index + 3] = (unsigned char)(255);
}

void SetCommand(hashtable_t *table, char *key, cmd_t cmd)
{
	hashtable_value_t value;
	value.pointer = (void *)cmd;
	HashTableSet(HashTableGetOrInsert(table, key), value);
}


int main(int argc, char **argv)
{
	hashtable_t *table_cmds = HashTableNewDefault();
	CommandsSetStandard(table_cmds);
	
	render_settings_t render_settings;
	render_settings_t *rs = &render_settings;
	rs->transform_ptr = rs->transform_stack;
	IdentityMatrixP(rs->transform_ptr);
	
	rs->triangles_ptr = NULL;
	rs->triangle_ptr = NULL;
	rs->materials_ptr = NULL;
	rs->material_ptr = NULL;
	rs->current_material_ptr = NULL;
	rs->lights_ptr = NULL;
	rs->light_ptr = NULL;
	rs->current_point = 0;
	rs->current_normal = 0;
	
	rs->focal_length = 1;
	rs->samples = 1;
	rs->shadow_samples = 1;
	rs->section_size = 256;
	rs->threads = 4;
	rs->max_iterations = 10;
	rs->input = stdin;
	memcpy(rs->out_file, DEFAULT_OUTPUT, sizeof(DEFAULT_OUTPUT));
	
	if(argc > 1)
		rs->input = fopen(argv[1], "r");

	while(1)
	{
	
		// INPUT
		char cmd_str[64];
		while(!fscanf(rs->input, "%s", cmd_str));
		hashtable_entry_t *entry = HashTableGet(table_cmds, cmd_str);
		if (entry != NULL)
		{
			cmd_t cmd = (cmd_t)entry->value.pointer;
			if (cmd(rs) == 0)
				break;
		}
		else
		{
			printf("error: unknown command \"%s\"\n", cmd_str);
		}
	}
	printf("info: freeing memory...\n");
	FreeMaterialList(rs->materials_ptr);
	FreeLightList(rs->lights_ptr);
	HashTableFree(table_cmds);
	printf("info: goodbye.\n");
	return 0;
}
