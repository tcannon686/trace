#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Windows.h>
#include <time.h>

#include "matrix.h"
#include "hashtable.h"
#include "trace.h"
#include "lodepng.h"
#include "shade.h"

#define DEFAULT_OUTPUT "output.png"
#define THRESHOLD 0.00001

int halton_index = 20;

vecc_t HaltonSequence(int i, int b)
{
	vecc_t f = 1.0;
	vecc_t r = 0.0;
	while(i > 0)
	{
		f = f / b;
		r = r + f * (i % b);
		i = floor(i / b);
	}
	
	return r;
}

vecc_t NextHalton(int b)
{
	halton_index ++;
	return HaltonSequence(halton_index, b);
}

void SeedHalton(int index)
{
	halton_index = index;
}

vector_t RandomVector(int axes)
{
	halton_index ++;
	vector_t ret = vec2(0, 0);
	for(int i = 0; i < axes; i ++)
	{
		ret.m[i] = 2 * HaltonSequence(halton_index, 2 + i) - 1;
	}
	return ret;
}

inline vector_t Barycentric(vector_t v, triangle_t tri)
{
	vector_t ret;
	vecc_t area_inv = 1.0 / tri.area;
	
	vector_t cross;
	
	cross = VectorCrossVector(
			VectorMinusVector(v, tri.a),
			tri.ba);
	
	if(VectorDotVector(tri.normal, cross) < 0)
	{
		ret.z = -1;
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
		ret.y = -1;
		return ret;
	}
	else
		ret.y = VectorMagnitude(cross) * area_inv;
	
	ret.w = 0;
	return ret;
}

inline int RayTriangle(hit_t *hit_ptr, ray_t *ray, triangle_t tri)
{
	vecc_t ldotn = VectorDotVector(ray->d, tri.normal);
	if(ldotn == 0)
		return 0;
	vecc_t d = VectorDotVector(VectorMinusVector(tri.v0, ray->o), tri.normal)
		/ ldotn;
	
	if(d > hit_ptr->t)
		return 0;
	
	if(d <= 0)
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

kd_tree_t *GenerateTree(tri_list_t *triangles_ptr, int depth) {

	// Return nothing if no triangles are provided
	if(triangles_ptr == NULL)
		return NULL;

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

	tree_ptr->left_ptr = GenerateTree(left_triangles_ptr, depth + 1);
	if(tree_ptr->left_ptr != NULL)
		tree_ptr->left_ptr->parent_ptr = tree_ptr;
	tree_ptr->right_ptr = GenerateTree(right_triangles_ptr, depth + 1);
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
		free(list_ptr);
		list_ptr = next_ptr;
	}
}

inline vecc_t RayBox(ray_t *ray, bounding_box_t box)
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
			return RayTriangles(hit_ptr, ray, tree_ptr->triangles_ptr);
		}
		
		int left, right;
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
	for(int i = rp_ptr->x0; i < rp_ptr->x1; i ++)
	{
		for(int j = rp_ptr->y0; j < rp_ptr->y1; j ++)
		{
			vector_t pixel = vec3(0, 0, 0);
			for(int k = 0; k < rp_ptr->samples; k ++)
			{
				vector_t sample_offset = 
					vec2((k % sampleWidth) * sampleSize, (k / sampleWidth) * sampleSize);
				ray_t ray = NewRay(
					vec3(0, 0, 0),
					vec3(rp_ptr->scale_x * (((vecc_t)i - rp_ptr->image.width / 2.0) + sample_offset.x),
						rp_ptr->scale_y * (((vecc_t)j - rp_ptr->image.height / 2.0) + sample_offset.y),
						-rp_ptr->focal_length));
				
				hit_t result;
				result.t = INFINITY;
				result.tree_ptr = rp_ptr->tree_ptr;
				vector_t sample = rp_ptr->sky_color;
				if(RayTree(&result, ray, rp_ptr->tree_ptr))
					sample = VectorClamp(result.material_ptr->shader(
						result,
						rp_ptr->lights_ptr, rp_ptr, 0, rp_ptr->max_iterations));
				VectorPlusVectorP(&pixel, &pixel, &sample);
			}
			VectorTimesScalarP(&pixel, &pixel, 1.0 / rp_ptr->samples);
			
			int pixel_index = (i + j * rp_ptr->image.width) * 4;
			rp_ptr->image.data[pixel_index + 0] = (unsigned char)(pixel.m[0] * 255);
			rp_ptr->image.data[pixel_index + 1] = (unsigned char)(pixel.m[1] * 255);
			rp_ptr->image.data[pixel_index + 2] = (unsigned char)(pixel.m[2] * 255);
			rp_ptr->image.data[pixel_index + 3] = (unsigned char)(255);
		}
	}
}

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


void Render(
	image_t image,
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
	HANDLE threads[64];
	render_params_list_t *rps = NULL;
	render_params_list_t *last_ptr = NULL;
	for(int i = 0; i < image.width; i += section_size)
	{
		for(int j = 0; j < image.height; j += section_size)
		{
			render_params_list_t *rp = (render_params_list_t *)malloc(sizeof(render_params_list_t));
			rp->is_rendering = 0;
			rp->next_ptr = NULL;
			rp->current.image = image;
			rp->current.tree_ptr = tree_ptr;
			rp->current.lights_ptr = lights_ptr;
			rp->current.sky_color = sky_color;
			rp->current.x0 = i;
			rp->current.y0 = j;
			rp->current.x1 = min(i + section_size, image.width);
			rp->current.y1 = min(j + section_size, image.height);
			rp->current.aspect = aspect;
			rp->current.scale_x = scale_x;
			rp->current.scale_y = scale_y;
			rp->current.focal_length = focal_length;
			rp->current.samples = samples;
			rp->current.shadow_samples = shadow_samples;
			rp->current.max_iterations = max_iterations;
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
			threads[i] = CreateThread(
				NULL, 0, RenderThread, (LPVOID)rps, 0, NULL);
			if(threads[i] == NULL)
			{
				fprintf(stderr, "error: failed to create thread.\n");
				return;
			}
		}
	
		WaitForMultipleObjects(num_threads, threads, TRUE, INFINITE);
	
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
			RenderThread((LPVOID)rps);
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

vecc_t MatGetNumber(material_t *material_ptr, char *key)
{
	return HashTableGet(material_ptr->table, key)->value.real;
}

int MatGetInteger(material_t *material_ptr, char *key)
{
	return HashTableGet(material_ptr->table, key)->value.integer;
}

vector_t MatGetVector(material_t *material_ptr, char *key)
{
	return HashTableGet(material_ptr->table, key)->value.vector;
}

void MatSetNumber(material_t *material_ptr, char *key, vecc_t value)
{
	HashTableGetOrInsert(material_ptr->table, key)->value.real = value;
}

void MatSetInteger(material_t *material_ptr, char *key, int value)
{
	HashTableGetOrInsert(material_ptr->table, key)->value.integer = value;
}

void MatSetVector(material_t *material_ptr, char *key, vector_t value)
{
	HashTableGetOrInsert(material_ptr->table, key)->value.vector = value;
}

int CmdQuit(render_settings_t *rs)
{
	return 0;
}

int CmdVertex(render_settings_t *rs)
{
	vecc_t x, y, z;
	vector_t vertex;

	int count = fscanf(rs->input, "%lf %lf %lf", &x, &y, &z);
	if (count == 0)
	{
		fprintf(stderr, "error: 'vertex' not enough arguments\n");
		return 1;
	}

	vertex = vec4(x, y, z, 1);

	MatrixTimesVectorP(
		&rs->triangle_ptr->current.v[rs->current_point],
		rs->transform_ptr,
		&vertex);

	rs->current_point++;
	if (rs->current_point > 2)
	{
		vector_t normal;
		normal = VectorNormalize(VectorCrossVector(
			VectorMinusVector(rs->triangle_ptr->current.v2, rs->triangle_ptr->current.v0),
			VectorMinusVector(rs->triangle_ptr->current.v1, rs->triangle_ptr->current.v0)));
		rs->triangle_ptr->current.normal = normal;
		rs->triangle_ptr->current.ba = VectorMinusVector(
			rs->triangle_ptr->current.b,
			rs->triangle_ptr->current.a);
		rs->triangle_ptr->current.cb = VectorMinusVector(
			rs->triangle_ptr->current.c,
			rs->triangle_ptr->current.b);
		rs->triangle_ptr->current.ac = VectorMinusVector(
			rs->triangle_ptr->current.a,
			rs->triangle_ptr->current.c);
		rs->triangle_ptr->current.area = VectorMagnitude(VectorCrossVector(
			rs->triangle_ptr->current.ba,
			rs->triangle_ptr->current.cb));
		/*if(current_normal == 0)
		{
		triangle_ptr->current.n[0] = normal;
		triangle_ptr->current.n[1] = normal;
		triangle_ptr->current.n[2] = normal;
		}*/

		rs->triangle_ptr->current.origin =
			VectorTimesScalar(
				VectorPlusVector(rs->triangle_ptr->current.v0,
					VectorPlusVector(rs->triangle_ptr->current.v1, rs->triangle_ptr->current.v2)),
				1.0 / 3.0);
		rs->triangle_ptr->current.material_ptr = rs->current_material_ptr;
		rs->current_point = 0;
	}
	return 1;
}

int CmdMakeFace(render_settings_t *rs)
{
	if (rs->triangle_ptr->current.area != 0)
	{
		if (rs->triangle_ptr->last_ptr != NULL)
			rs->triangle_ptr->last_ptr->next_ptr = rs->triangle_ptr;
		rs->triangle_ptr->next_ptr = malloc(sizeof(tri_list_t));
		memset(rs->triangle_ptr->next_ptr, 0, sizeof(tri_list_t));
		rs->triangle_ptr->next_ptr->last_ptr = rs->triangle_ptr;
		rs->triangle_ptr = rs->triangle_ptr->next_ptr;
		rs->triangle_ptr->last_ptr->next_ptr = NULL;
	}
	return 1;
}

int CmdNormal(render_settings_t *rs)
{
	vecc_t x, y, z;
	vector_t normal;

	int count = fscanf(rs->input, "%lf %lf %lf", &x, &y, &z);
	if (count == 0)
	{
		fprintf(stderr, "error: 'normal' not enough arguments.\n");
		return 1;
	}
	normal = vec4(x, y, z, 0);

	MatrixTimesVectorP(
		&rs->triangle_ptr->current.n[rs->current_normal],
		rs->transform_ptr,
		&normal);

	normal = VectorNormalize(normal);
	rs->current_normal++;

	if (rs->current_normal > 2)
		rs->current_normal = 0;
	return 1;
}

int CmdSky(render_settings_t *rs)
{
	vecc_t x, y, z;

	int count = fscanf(rs->input, "%lf %lf %lf", &x, &y, &z);
	if (count == 0)
	{
		fprintf(stderr, "error: 'sky' not enough arguments.\n");
		return 1;
	}

	rs->sky = vec4(x, y, z, 1);
	return 1;
}

int CmdCamFov(render_settings_t *rs)
{
	vecc_t fov = 3.14159265 / 2;

	int count = fscanf(rs->input, "%lf", &fov);
	if (count == 0)
	{
		fprintf(stderr, "error: 'cam_fov' not enough arguments.\n");
		return 1;
	}

	rs->focal_length = 1.0 / (2.0 * tan(fov / 4.0));
	return 1;
}

int CmdMatIndex(render_settings_t *rs)
{
	int index = 0;
	mat_list_t *current_ptr = rs->materials_ptr;
	int count = fscanf(rs->input, "%i", &index);
	if (count == 0)
	{
		fprintf(stderr, "error: 'mat_index' not enough arguments.\n");
		return 1;
	}
	for (int i = 0; i < index; i++)
		current_ptr = current_ptr->next_ptr;
	rs->current_material_ptr = &(current_ptr->current);
	return 1;
}

int CmdMatSetNumber(render_settings_t *rs)
{
	char key[32];
	vecc_t value;
	int count = fscanf(rs->input, "%s %lf", key, &value);
	if (count == 0)
	{
		fprintf(stderr, "error: 'mat_set_number' not enough arguments.\n");
		return 1;
	}
	MatSetNumber(&(rs->material_ptr->current), key, value);
	return 1;
}

int CmdMatSetInteger(render_settings_t *rs)
{
	char key[32];
	int value;
	int count = fscanf(rs->input, "%s %i", key, &value);
	if (count == 0)
	{
		fprintf(stderr, "error: 'mat_set_integer' not enough arguments.\n");
		return 1;
	}
	MatSetInteger(&(rs->material_ptr->current), key, value);
	return 1;
}

int CmdMatSetVector(render_settings_t *rs)
{
	char key[32];
	vector_t value;
	int count = fscanf(rs->input, "%s %lf %lf %lf", key, &value.x, &value.y, &value.z);
	if (count == 0)
	{
		fprintf(stderr, "error: 'mat_set_vector' not enough arguments.\n");
		return 1;
	}
	MatSetVector(&(rs->material_ptr->current), key, value);
	return 1;
}

int CmdMatShader(render_settings_t *rs)
{
	char key[32];
	int count = fscanf(rs->input, "%s", key);
	if (count == 0)
	{
		fprintf(stderr, "error: 'mat_set_vector' not enough arguments.\n");
		return 1;
	}
	if (strcmp(key, "phong"))
		rs->material_ptr->current.shader = PhongShader;
	return 1;
}

int CmdLightPosition(render_settings_t *rs)
{
	vecc_t x, y, z;
	vector_t position;
	int count = fscanf(rs->input, "%lf %lf %lf", &x, &y, &z);
	if (count == 0)
	{
		fprintf(stderr, "error: 'light_position' not enough arguments.\n");
		return 1;
	}
	position = vec4(x, y, z, 1);
	MatrixTimesVectorP(&rs->light_ptr->current.position, rs->transform_ptr, &position);
	return 1;
}

int CmdLightColor(render_settings_t *rs)
{
	vecc_t x, y, z;

	int count = fscanf(rs->input, "%lf %lf %lf", &x, &y, &z);
	if (count == 0)
	{
		fprintf(stderr, "error: 'light_color' not enough arguments.\n");
		return 1;
	}

	rs->light_ptr->current.color = vec3(x, y, z);
	return 1;
}

int CmdLightDistance(render_settings_t *rs)
{
	vecc_t x;

	int count = fscanf(rs->input, "%lf", &x);
	if (count == 0)
	{
		fprintf(stderr, "error: 'light_distance' not enough arguments.\n");
		return 1;
	}

	rs->light_ptr->current.distance = x;
	return 1;
}

int CmdLightEnergy(render_settings_t *rs)
{
	vecc_t x;

	int count = fscanf(rs->input, "%lf", &x);
	if (count == 0)
	{
		fprintf(stderr, "error: 'light_energy' not enough arguments.\n");
		return 1;
	}

	rs->light_ptr->current.energy = x;
	return 1;
}

int CmdMakeMaterial(render_settings_t *rs)
{
	rs->material_ptr->next_ptr = malloc(sizeof(mat_list_t));
	memset(rs->material_ptr->next_ptr, 0, sizeof(mat_list_t));
	rs->material_ptr->next_ptr->last_ptr = rs->material_ptr;
	rs->material_ptr = rs->material_ptr->next_ptr;
	rs->material_ptr->current.shader = PhongShader;
	rs->material_ptr->current.table = HashTableNewDefault();
	rs->current_material_ptr = &rs->material_ptr->current;
	return 1;
}

int CmdMakeLight(render_settings_t *rs)
{
	rs->light_ptr->next_ptr = malloc(sizeof(light_list_t));
	memset(rs->light_ptr->next_ptr, 0, sizeof(light_list_t));
	rs->light_ptr->next_ptr->last_ptr = rs->light_ptr;
	rs->light_ptr = rs->light_ptr->next_ptr;
	return 1;
}

int CmdInFile(render_settings_t *rs)
{
	char file[255];
	fscanf(rs->input, "%s", file);
	rs->input = fopen(file, "r");
	if (rs->input == NULL) {
		printf("error: file not found.\n");
		rs->input = stdin;
	}
	return 1;
}

int CmdOutFile(render_settings_t *rs)
{
	fscanf(rs->input, "%s", rs->out_file);
	return 1;
}

int CmdInStdIn(render_settings_t *rs)
{
	rs->input = stdin;
	return 1;
}

int CmdRenderSamples(render_settings_t *rs)
{
	int count = fscanf(rs->input, "%i", &rs->samples);
	if (count == 0)
	{
		fprintf(stderr, "error: 'render_samples' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdRenderShadowSamples(render_settings_t *rs)
{
	int count = fscanf(rs->input, "%i", &rs->shadow_samples);
	if (count == 0)
	{
		fprintf(stderr, "error: 'render_samples' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdRenderSectionSize(render_settings_t *rs)
{
	int count = fscanf(rs->input, "%i", &rs->section_size);
	if (count == 0)
	{
		fprintf(stderr, "error: 'render_section_size' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdRenderThreads(render_settings_t *rs)
{
	int count = fscanf(rs->input, "%i", &rs->threads);
	if (count == 0)
	{
		fprintf(stderr, "error: 'render_threads' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdRenderIterations(render_settings_t *rs)
{
	int count = fscanf(rs->input, "%i", &rs->max_iterations);
	if (count == 0)
	{
		fprintf(stderr, "error: 'render_iterations' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdRender(render_settings_t *rs)
{
	unsigned int width;
	unsigned int height;
	image_t image;
	int count = fscanf(rs->input, "%u %u", &width, &height);
	if (count < 2)
	{
		fprintf(stderr, "error: 'render' not enough arguments.\n");
		return 1;
	}

	vecc_t aspect = (vecc_t)height / (vecc_t)width;
	vecc_t scale_x = 2 / (vecc_t)width;
	vecc_t scale_y = -aspect * 2 / (vecc_t)height;
	image.width = width;
	image.height = height;
	image.data = (unsigned char *)malloc(width * height * 4);

	clock_t start_time, end_time;
	float elapsed_secs;

	start_time = clock();

	printf("info: generating tree\n");
	kd_tree_t *tree_ptr = GenerateTree(rs->triangles_ptr, 0);
	end_time = clock();
	elapsed_secs = (float)(end_time - start_time) / CLOCKS_PER_SEC;
	start_time = end_time;
	printf("info: elapsed time %f seconds.\n", elapsed_secs);
	printf("info: rendering\n");
	/*for(int i = 0; i < width; i ++)
	{
	for(int j = 0; j < height; j ++)
	{
	int index = (i + j * width) * 4;
	vector_t average = vec4(0, 0, 0, 0);
	for(int k = 0; k < samples; k ++)
	{
	for(int l = 0; l < samples; l ++)
	{
	vecc_t x = (vecc_t)(i - (int)width / 2) * scale_x
	+ k * scale_x / samples;
	vecc_t y = (vecc_t)((int)height / 2 - j) * scale_y
	+ l * scale_y / samples;
	vecc_t z = -focal_length;

	ray_t ray = NewRay(vec3(0, 0, 0), vec3(x, y, z));

	hit_t result = RayTree(ray, tree_ptr);
	if(result.t > 0)
	{
	vector_t color = result.material_ptr->shader(result, lights_ptr);
	color = VectorClamp(color);
	average = VectorPlusVector(average, color);
	}
	else
	{
	average = VectorPlusVector(average, sky);
	}
	}
	}
	average = VectorTimesScalar(average, 1.0f / (samples * samples));
	image[index + 0] = (unsigned char) (average.x * 255);
	image[index + 1] = (unsigned char) (average.y * 255);
	image[index + 2] = (unsigned char) (average.z * 255);
	image[index + 3] = (unsigned char) (255);
	}
	} */

	Render(
		image,
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

	end_time = clock();
	elapsed_secs = (float)(end_time - start_time) / CLOCKS_PER_SEC;
	printf("info: elapsed time %f seconds.\n", elapsed_secs);

	printf("info: saving file '%s'.\n", rs->out_file);
	unsigned error = lodepng_encode32_file(rs->out_file, image.data, width, height);

	/*if there's an error, display it*/
	if (error)
		printf("lodepng error %u: %s\n", error, lodepng_error_text(error));
	else
		printf("okay\n");
	FreeTree(tree_ptr);
	return 1;
}

int CmdPrintTriangles(render_settings_t *rs)
{
	printf("triangles: \n");
	PrintTriangles(rs->triangles_ptr);
	return 1;
}

int CmdTransformPop(render_settings_t *rs)
{
	rs->transform_ptr--;
	return 1;
}

int CmdTransformPush(render_settings_t *rs)
{
	rs->transform_ptr++;
	memcpy(rs->transform_ptr, rs->transform_ptr - 1, sizeof(matrix_t));
	return 1;
}

int CmdTransformTranslate(render_settings_t *rs)
{
	vector_t t;
	matrix_t translate;
	fscanf(rs->input, "%lf %lf %lf", &t.x, &t.y, &t.z);
	translate = NewTranslateMatrix(t);
	*rs->transform_ptr = MatrixTimesMatrix(*rs->transform_ptr, translate);
	return 1;
}

int CmdTransformScale(render_settings_t *rs)
{
	vector_t t;
	matrix_t scale;
	fscanf(rs->input, "%lf %lf %lf", &t.x, &t.y, &t.z);
	scale = NewScaleMatrix(t);
	*rs->transform_ptr = MatrixTimesMatrix(*rs->transform_ptr, scale);
	return 1;
}

int CmdTransformRotate(render_settings_t *rs)
{
	vecc_t angle;
	vector_t t;
	matrix_t rotate;
	fscanf(rs->input, "%lf, %lf %lf %lf", &angle, &t.x, &t.y, &t.z);
	rotate = NewRotateMatrix(angle, t);
	*rs->transform_ptr = MatrixTimesMatrix(*rs->transform_ptr, rotate);
	return 1;
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
	SetCommand(table_cmds, "quit", CmdQuit);
	SetCommand(table_cmds, "vertex", CmdVertex);
	SetCommand(table_cmds, "normal", CmdNormal);
	SetCommand(table_cmds, "make_face", CmdMakeFace);
	SetCommand(table_cmds, "sky", CmdSky);
	SetCommand(table_cmds, "cam_fov", CmdCamFov);
	SetCommand(table_cmds, "mat_index", CmdMatIndex);
	SetCommand(table_cmds, "mat_set_number", CmdMatSetNumber);
	SetCommand(table_cmds, "mat_set_integer", CmdMatSetInteger);
	SetCommand(table_cmds, "mat_set_vector", CmdMatSetVector);
	SetCommand(table_cmds, "mat_shader", CmdMatShader);
	SetCommand(table_cmds, "light_position", CmdLightPosition);
	SetCommand(table_cmds, "light_color", CmdLightColor);
	SetCommand(table_cmds, "light_energy", CmdLightEnergy);
	SetCommand(table_cmds, "light_distance", CmdLightDistance);
	SetCommand(table_cmds, "make_material", CmdMakeMaterial);
	SetCommand(table_cmds, "make_light", CmdMakeLight);
	SetCommand(table_cmds, "in_file", CmdInFile);
	SetCommand(table_cmds, "out_file", CmdOutFile);
	SetCommand(table_cmds, "in_stdin", CmdInStdIn);
	SetCommand(table_cmds, "render_samples", CmdRenderSamples);
	SetCommand(table_cmds, "render_shadow_samples", CmdRenderShadowSamples);
	SetCommand(table_cmds, "render_section_size", CmdRenderSectionSize);
	SetCommand(table_cmds, "render_threads", CmdRenderThreads);
	SetCommand(table_cmds, "render_iterations", CmdRenderIterations);
	SetCommand(table_cmds, "render", CmdRender);
	SetCommand(table_cmds, "print_triangles", CmdPrintTriangles);
	SetCommand(table_cmds, "transform_pop", CmdTransformPop);
	SetCommand(table_cmds, "transform_push", CmdTransformPush);
	SetCommand(table_cmds, "transform_translate", CmdTransformTranslate);
	SetCommand(table_cmds, "transform_scale", CmdTransformScale);
	SetCommand(table_cmds, "transform_rotate", CmdTransformRotate);
	
	render_settings_t render_settings;
	render_settings_t *rs = &render_settings;
	rs->transform_ptr = rs->transform_stack;
	IdentityMatrixP(rs->transform_ptr);
	
	rs->triangles_ptr = (tri_list_t *)malloc(sizeof(tri_list_t));
	rs->triangle_ptr = rs->triangles_ptr;
	rs->materials_ptr = (mat_list_t *)malloc(sizeof(mat_list_t));
	rs->material_ptr = rs->materials_ptr;
	rs->current_material_ptr = &rs->material_ptr->current;
	rs->lights_ptr = (light_list_t *)malloc(sizeof(light_list_t));
	rs->light_ptr = rs->lights_ptr;
	rs->current_point = 0;
	rs->current_normal = 0;

	memset(rs->triangles_ptr, 0, sizeof(tri_list_t));
	memset(rs->materials_ptr, 0, sizeof(mat_list_t));
	memset(rs->lights_ptr, 0, sizeof(light_list_t));

	rs->current_material_ptr->table = HashTableNewDefault();
	rs->focal_length = 1;
	rs->samples = 1;
	rs->shadow_samples = 1;
	rs->section_size = 256;
	rs->threads = 4;
	rs->max_iterations = 10;
	rs->input = stdin;
	memcpy(rs->out_file, DEFAULT_OUTPUT, sizeof(DEFAULT_OUTPUT));
	
	rs->material_ptr->current.shader = PhongShader;
	
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
