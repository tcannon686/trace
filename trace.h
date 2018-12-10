
#ifndef TRACE_H
#define TRACE_H

struct kd_tree;
typedef struct kd_tree kd_tree_t;

struct material;
typedef struct material material_t;

typedef struct
{
	vector_t o;
	vector_t d;
	vector_t d_inverse;
} ray_t;

typedef struct
{
	vecc_t t;
	vector_t position;
	vector_t normal;
	material_t *material_ptr;
	ray_t ray;
	kd_tree_t *tree_ptr;
	kd_tree_t *node_ptr;
} hit_t;


typedef struct
{
	vector_t position;
	hashtable_t *table;
} light_t;

struct light_list;
typedef struct light_list light_list_t;
typedef struct light_list
{
	light_t current;
	light_list_t *last_ptr;
	light_list_t *next_ptr;
} light_list_t;

struct render_params_list;
typedef struct render_params_list render_params_list_t;


typedef struct
{
	unsigned char *data;
	unsigned int width, height;
} image_t;

typedef void (*pixel_traced_callback_t)(int x, int y, vector_t color, void *data);

typedef struct render_params
{
	kd_tree_t *tree_ptr;
	light_list_t *lights_ptr;
	vector_t sky_color;
	int x0, y0;
	int x1, y1;
	int width, height;
	vecc_t aspect;
	vecc_t scale_x, scale_y;
	vecc_t focal_length;
	int samples;
	int shadow_samples;
	int max_iterations;
	pixel_traced_callback_t callback;
	void *callback_data;
} render_params_t;

typedef vector_t (*shader_t)(
	void *shader_data,
	hit_t *hit,
	render_params_t *rp_ptr,
	int iteration,
	int max_iterations);

#define ShadeHit(hit, rp_ptr, iteration)\
	VectorClamp(hit.material_ptr->shader(\
		hit.material_ptr->shader_data, &hit,\
		rp_ptr, iteration, rp_ptr->max_iterations));

typedef struct material
{
	shader_t shader;
	void *shader_data;
	hashtable_t *table;
} material_t;

struct mat_list;
typedef struct mat_list mat_list_t;
typedef struct mat_list
{
	material_t current;
	mat_list_t *last_ptr;
	mat_list_t *next_ptr;
} mat_list_t;

typedef struct
{
	union
	{
		vector_t v[3];
		struct
		{
			vector_t v0;
			vector_t v1;
			vector_t v2;
		};
	
		struct
		{
			vector_t a;
			vector_t b;
			vector_t c;
		};
	};
	
	union
	{
		vector_t n[3];
		struct
		{
			vector_t n0;
			vector_t n1;
			vector_t n2;
		};
	
		struct
		{
			vector_t na;
			vector_t nb;
			vector_t nc;
		};
	};
	
	union
	{
		struct
		{
			vector_t ba;
			vector_t cb;
			vector_t ac;
		};
		
		struct
		{
			vector_t v1v0;
			vector_t v2v1;
			vector_t v0v2;
		};
	};
	
	vector_t origin;
	vector_t normal;
	
	vecc_t area;
	
	material_t *material_ptr;
} triangle_t;

struct tri_list;
typedef struct tri_list tri_list_t;

typedef struct tri_list
{
	triangle_t current;
	tri_list_t *last_ptr;
	tri_list_t *next_ptr;
} tri_list_t;


typedef struct
{
	vector_t a;
	vector_t b;
	vector_t origin;
} bounding_box_t;

typedef struct kd_tree
{
	vector_t median;
	vector_t axis;
	int axis_index;
	kd_tree_t *left_ptr;
	kd_tree_t *right_ptr;
	kd_tree_t *parent_ptr;
	tri_list_t *triangles_ptr;
	bounding_box_t box;
} kd_tree_t;

typedef struct render_params_list
{
	int is_rendering;
	render_params_t current;
	render_params_list_t *next_ptr;
} render_params_list_t;

typedef struct render_settings
{
	matrix_t transform_stack[32];
	matrix_t *transform_ptr;
	vector_t sky;
	tri_list_t *triangles_ptr;
	tri_list_t *triangle_ptr;
	mat_list_t *materials_ptr;
	mat_list_t *material_ptr;
	material_t *current_material_ptr;
	light_list_t *lights_ptr;
	light_list_t *light_ptr;
	int current_point;
	int current_normal;
	vecc_t focal_length;
	int samples;
	int shadow_samples;
	int section_size;
	int threads;
	int max_iterations;
	FILE *input;
	char out_file[255];
} render_settings_t;

// Returns 1 if program execution should continue, otherwise some other number.
typedef int(*cmd_t)(render_settings_t *rs);

void SetCommand(hashtable_t *table, char *key, cmd_t cmd);

vector_t RandomVector(int axes);

ray_t NewRay(vector_t o, vector_t d);

kd_tree_t *GenerateTree(tri_list_t *triangles_ptr, int depth);
void CleanTriangles(tri_list_t **triangles_ptr);

int RayTriangle(hit_t *hit_ptr, ray_t *ray, triangle_t tri);
int RayTriangles(hit_t *hit_ptr, ray_t *ray, tri_list_t *triangles_ptr);
int RayTree(hit_t *hit_ptr, ray_t ray, kd_tree_t *tree_ptr);

vecc_t RayBox(ray_t *ray, bounding_box_t box);
vector_t Barycentric(vector_t v, triangle_t tri);

void SeedHalton(int index);

vector_t VectorMix(vector_t left, vector_t right, vecc_t a);

void PrintTriangles(tri_list_t *triangles_ptr);

void FreeTriList(tri_list_t *list_ptr);
void FreeTree(kd_tree_t *tree_ptr);
void FreeMaterialList(mat_list_t *list_ptr);
void FreeLightList(light_list_t *list_ptr);

void RenderImageCallback(int x, int y, vector_t color, void *data);

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
	int num_threads);

void SetCommand(hashtable_t *table, char *key, cmd_t cmd);

#define PropGet(material_ptr, key)	HashTableGet((material_ptr)->table, key)->value
#define PropGetOrInsert(material_ptr, key)	HashTableGetOrInsert((material_ptr)->table, key)->value

#endif
