
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
} hit_t;


typedef struct
{
	vector_t position;
	vector_t color;
	vecc_t distance;
	vecc_t energy;
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

typedef struct render_params
{
	image_t image;
	kd_tree_t *tree_ptr;
	light_list_t *lights_ptr;
	vector_t sky_color;
	int x0, y0;
	int x1, y1;
	vecc_t aspect;
	vecc_t scale_x, scale_y;
	vecc_t focal_length;
	int samples;
	int max_iterations;
} render_params_t;

typedef vector_t (*shader_t)(
	hit_t hit,
	light_list_t *lights_ptr,
	render_params_t *rp_ptr,
	int iteration,
	int max_iterations);

typedef struct material
{
	shader_t shader;
	/*vector_t specular;
	vector_t diffuse;
	vector_t ambient;
	vecc_t shininess;
	vecc_t reflectiveness;
	vecc_t alpha;
	vecc_t ior;
	int shadeless;*/
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

vector_t RandomVector(int axes);

ray_t NewRay(vector_t o, vector_t d);

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

vecc_t MatGetNumber(material_t *material_ptr, char *name);
void MatSetNumber(material_t *material_ptr, char *name, vecc_t value);
int MatGetInteger(material_t *material_ptr, char *name);
void MatSetInteger(material_t *material_ptr, char *name, int value);
vector_t MatGetVector(material_t *material_ptr, char *name);
void MatSetVector(material_t *material_ptr, char *name, vector_t value);

#endif