
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "matrix.h"
#include "hashtable.h"
#include "linkedlist.h"
#include "trace.h"
#include "lodepng.h"
#include "shade.h"
#include "xwin.h"
#include "cmd.h"
#include "halton.h"


#ifdef INCLUDE_GUI
int CmdRenderWindow(render_settings_t *rs)
{
	unsigned int width;
	unsigned int height;
	int count = fscanf(rs->input, "%u %u", &width, &height);
	if (count < 2)
	{
		fprintf(stderr, "error: 'render_window' not enough arguments.\n");
		return 1;
	}
	
	WindowMain(rs, width, height);
	return 1;
}
#endif

int CmdQuit(render_settings_t *rs)
{
	return 0;
}

int CmdVertex(render_settings_t *rs)
{
	vector_t vertex;

	if (!ReadVec3(rs, &vertex))
	{
		fprintf(stderr, "error: 'vertex' not enough arguments\n");
		return 1;
	}

	primitive_t *primitive;
	triangle_t *triangle_ptr;
	primitive = (primitive_t *)ListGetLast(rs->primitives_ptr).data_ptr;
	if(primitive->type != rs->triangle_class_ptr)
	{
		fprintf(stderr, "error: 'vertex' last primitive is not a triangle.");
		return 1;
	}

	triangle_ptr = (triangle_t *) primitive->data;

	MatrixTimesVectorP(
		&triangle_ptr->v[rs->current_point],
		rs->transform_ptr,
		&vertex);

	rs->current_point++;
	if (rs->current_point > 2)
	{
		vector_t normal;
		normal = VectorNormalize(VectorCrossVector(
			VectorMinusVector(triangle_ptr->v2, triangle_ptr->v0),
			VectorMinusVector(triangle_ptr->v1, triangle_ptr->v0)));
		triangle_ptr->normal = normal;
		triangle_ptr->ba = VectorMinusVector(
			triangle_ptr->b,
			triangle_ptr->a);
		triangle_ptr->cb = VectorMinusVector(
			triangle_ptr->c,
			triangle_ptr->b);
		triangle_ptr->ac = VectorMinusVector(
			triangle_ptr->a,
			triangle_ptr->c);
		triangle_ptr->area = VectorMagnitude(VectorCrossVector(
			triangle_ptr->ba,
			triangle_ptr->cb));
		/*if(current_normal == 0)
		{
		triangle_ptr->current.n[0] = normal;
		triangle_ptr->current.n[1] = normal;
		triangle_ptr->current.n[2] = normal;
		}*/

		triangle_ptr->origin =
			VectorTimesScalar(
				VectorPlusVector(triangle_ptr->v0,
					VectorPlusVector(triangle_ptr->v1, triangle_ptr->v2)),
				1.0 / 3.0);
		primitive->origin = triangle_ptr->origin;
		triangle_ptr->material_ptr = rs->current_material_ptr;
		rs->current_point = 0;
	}
	return 1;
}

int CmdMakeFace(render_settings_t *rs)
{
	primitive_t *primitive_ptr;
	primitive_ptr = (primitive_t *)malloc(sizeof(triangle_t));
	primitive_ptr->type = rs->triangle_class_ptr;
	primitive_ptr->data = malloc(sizeof(triangle_t));
	ListAppendPointer(rs->primitives_ptr, primitive_ptr);
	return 1;
}

int CmdNormal(render_settings_t *rs)
{
	vector_t normal;


	if (!ReadVec3(rs, &normal))
	{
		fprintf(stderr, "error: 'normal' not enough arguments.\n");
		return 1;
	}
	normal.w = 0;

	primitive_t *primitive;
	triangle_t *triangle_ptr;
	primitive = (primitive_t *)ListGetLast(rs->primitives_ptr).data_ptr;
	if(primitive->type != rs->triangle_class_ptr)
	{
		fprintf(stderr, "error: 'normal' last primitive is not a triangle.");
		return 1;
	}

	triangle_ptr = (triangle_t *) primitive->data;

	MatrixTimesVectorP(
		&triangle_ptr->n[rs->current_normal],
		rs->transform_ptr,
		&normal);

	normal = VectorNormalize(normal);
	rs->current_normal++;

	if (rs->current_normal > 2)
		rs->current_normal = 0;
	return 1;
}

int CmdTexCo2d(render_settings_t *rs)
{
	vector_t co;

	if (!ReadVec2(rs, &co))
	{
		fprintf(stderr, "error: 'texco2d' not enough arguments.\n");
		return 1;
	}
	
	primitive_t *primitive;
	triangle_t *triangle_ptr;
	primitive = (primitive_t *)ListGetLast(rs->primitives_ptr).data_ptr;
	if(primitive->type != rs->triangle_class_ptr)
	{
		fprintf(stderr, "error: 'vertex' last primitive is not a triangle.");
		return 1;
	}

	triangle_ptr = (triangle_t *) primitive->data;

	triangle_ptr->t[rs->current_texco] = co;
	
	rs->current_texco++;
	if (rs->current_texco > 2)
		rs->current_texco = 0;
	return 1;
}

int CmdSceneSkyColor(render_settings_t *rs)
{
	if (!ReadVec3(rs, &rs->scene_ptr->sky_color))
	{
		fprintf(stderr, "error: 'scene_sky_color' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdSceneSetVector(render_settings_t *rs)
{
	char key[32];
	if(!ReadName(rs, key, sizeof(key)))
	{
		fprintf(stderr, "error: 'scene_set_vector' expected name.\n");
		return 1;
	}

	if (!ReadVec3(rs, &PropGetOrInsert(rs->scene_ptr, key).vector))
	{
		fprintf(stderr, "error: 'scene_set_vector' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdSceneSetNumber(render_settings_t *rs)
{
	char key[32];

	if(!ReadName(rs, key, sizeof(key)))
	{
		fprintf(stderr, "error: 'scene_set_number' expected name.\n");
		return 1;
	}
	
	if (!ReadNumber(rs, &PropGetOrInsert(rs->scene_ptr, key).number))
	{
		fprintf(stderr, "error: 'scene_set_number' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdSceneSetInteger(render_settings_t *rs)
{
	char key[32];
	if(!ReadName(rs, key, sizeof(key)))
	{
		fprintf(stderr, "error: 'scene_set_integer' expected name.\n");
		return 1;
	}
	if (!ReadInteger(rs, &PropGetOrInsert(rs->scene_ptr, key).integer))
	{
		fprintf(stderr, "error: 'scene_set_integer' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdCamFov(render_settings_t *rs)
{
	vecc_t fov = 3.14159265 / 2;

	if(!ReadNumber(rs, &fov))
	{
		fprintf(stderr, "error: 'cam_fov' expected integer.\n");
		return 1;
	}

	rs->focal_length = 1.0 / (2.0 * tan(fov / 4.0));
	return 1;
}

int CmdMatIndex(render_settings_t *rs)
{
	int index = 0;
	if (!ReadInteger(rs, &index))
	{
		fprintf(stderr, "error: 'mat_index' not enough arguments.\n");
		return 1;
	}
	rs->current_material_ptr = (material_t *)ListGet(rs->materials_ptr, index).data_ptr;
	return 1;
}

int CmdMatSetNumber(render_settings_t *rs)
{
	char key[32];
	if(!ReadName(rs, key, sizeof(key)))
	{
		
		fprintf(stderr, "error: 'mat_set_number' expected name.\n");
		return 1;
	}
	
	if (!ReadNumber(rs, &PropGetOrInsert(rs->material_ptr, key).number))
	{
		fprintf(stderr, "error: 'mat_set_number' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdMatSetInteger(render_settings_t *rs)
{
	char key[32];
	if(!ReadName(rs, key, sizeof(key)))
	{
		fprintf(stderr, "error: 'mat_set_integer' expected name.\n");
		return 1;
	}
	if (!ReadInteger(rs, &PropGetOrInsert(rs->material_ptr, key).integer))
	{
		fprintf(stderr, "error: 'mat_set_integer' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdMatSetVector(render_settings_t *rs)
{
	char key[32];
	if(!ReadName(rs, key, sizeof(key)))
	{
		fprintf(stderr, "error: 'mat_set_vector' expected name.\n");
		return 1;
	}
	
	if (!ReadVec3(rs, &PropGetOrInsert(rs->material_ptr, key).vector))
	{
		fprintf(stderr, "error: 'mat_set_vector' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdMatSetTexture2d(render_settings_t *rs)
{
	char key[32];
	int value;
	
	if(!ReadName(rs, key, sizeof(key)))
	{
		fprintf(stderr, "error: 'mat_set_texture_2d' not enough arguments.\n");
		return 1;
	}

	if (!ReadInteger(rs, &value))
	{
		fprintf(stderr, "error: 'mat_set_texture_2d' not enough arguments.\n");
		return 1;
	}
	PropGetOrInsert(rs->material_ptr, key).pointer =
		ListGet(rs->texture2ds_ptr, value).data_ptr;
	return 1;
}

int CmdMatShader(render_settings_t *rs)
{
	char key[32];

	if (!ReadName(rs, key, sizeof(key)))
	{
		fprintf(stderr, "error: 'mat_shader' not enough arguments.\n");
		return 1;
	}
	if (strcmp(key, "phong"))
		rs->material_ptr->shader = PhongShader;
	else
	{
		fprintf(stderr, "error: '%s' shader not found.\n", key);
	}
	return 1;
}

int CmdLightPosition(render_settings_t *rs)
{
	vector_t position;

	if (!ReadVec3(rs, &position))
	{
		fprintf(stderr, "error: 'light_position' not enough arguments.\n");
		return 1;
	}
	MatrixTimesVectorP(&rs->light_ptr->position, rs->transform_ptr, &position);
	return 1;
}

int CmdLightSetVector(render_settings_t *rs)
{
	char key[32];

	if (!ReadName(rs, key, sizeof(key)))
	{
		fprintf(stderr, "error: 'light_set_vector' expected name.\n");
		return 1;
	}

	if (!ReadVec3(rs, &PropGetOrInsert(rs->light_ptr, key).vector))
	{
		fprintf(stderr, "error: 'light_set_vector' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdLightSetNumber(render_settings_t *rs)
{
	char key[32];
	
	if (!ReadName(rs, key, sizeof(key)))
	{
		fprintf(stderr, "error: 'light_set_number' expected name.\n");
		return 1;
	}

	if (!ReadNumber(rs, &PropGetOrInsert(rs->light_ptr, key).number))
	{
		fprintf(stderr, "error: 'light_set_number' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdLightSetInteger(render_settings_t *rs)
{
	char key[32];
	
	if (!ReadName(rs, key, sizeof(key)))
	{
		fprintf(stderr, "error: 'light_set_integers' expected name.\n");
		return 1;
	}

	if (!ReadInteger(rs, &PropGetOrInsert(rs->light_ptr, key).integer))
	{
		fprintf(stderr, "error: 'light_set_integers' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdMakeMaterial(render_settings_t *rs)
{
	material_t *material_ptr = malloc(sizeof(material_t));
	material_ptr->shader = PhongShader;
	material_ptr->shader_data = NULL;
	material_ptr->table = HashTableNewDefault();
	rs->material_ptr = material_ptr;
	rs->current_material_ptr = material_ptr;
	ListAppendPointer(rs->materials_ptr, (void *)material_ptr);
	/*if(rs->material_ptr == NULL)
	{
		material_t *material_ptr = malloc(sizeof(mat_list_t));
		memset(rs->material_ptr, 0, sizeof(mat_list_t));
		material_ptr->current.shader = PhongShader;
		material_ptr->current.shader_data = NULL;
		material_ptr->current.table = HashTableNewDefault();
		rs->current_material_ptr = &rs->material_ptr->current;
		rs->materials_ptr = rs->material_ptr;
	}
	else
	{
		rs->material_ptr->next_ptr = malloc(sizeof(mat_list_t));
		memset(rs->material_ptr->next_ptr, 0, sizeof(mat_list_t));
		rs->material_ptr->next_ptr->last_ptr = rs->material_ptr;
		rs->material_ptr = rs->material_ptr->next_ptr;
		rs->material_ptr->current.shader = PhongShader;
		rs->material_ptr->current.table = HashTableNewDefault();
		rs->current_material_ptr = &rs->material_ptr->current;
	}*/
	return 1;
}

int CmdMakeLight(render_settings_t *rs)
{
	light_t *light_ptr;
	light_ptr = (light_t *)malloc(sizeof(light_t));
	memset(light_ptr, 0, sizeof(light_t));
	light_ptr->table = HashTableNewDefault();
	rs->light_ptr = light_ptr;
	ListAppendPointer(rs->lights_ptr, (void *)light_ptr);
	return 1;
}

int CmdInFile(render_settings_t *rs)
{
	char file[255];
	if(!ReadString(rs, file, sizeof(file)))
	{
		fprintf(stderr, "error: 'in_file' expected string.\n");
		/* Exit the program. */
		return 0;
	}
	rs->input = fopen(file, "r");
	if (rs->input == NULL) {
		fprintf(stderr, "error: file not found.\n");
		rs->input = stdin;
	}
	return 1;
}

int CmdInStdIn(render_settings_t *rs)
{
	rs->input = stdin;
	return 1;
}

int CmdRenderSamples(render_settings_t *rs)
{
	if (!ReadInteger(rs, &rs->samples))
	{
		fprintf(stderr, "error: 'render_samples' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdRenderSectionSize(render_settings_t *rs)
{
	if (!ReadInteger(rs, &rs->section_size))
	{
		fprintf(stderr, "error: 'render_section_size' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdRenderThreads(render_settings_t *rs)
{
	if (!ReadInteger(rs, &rs->threads))
	{
		fprintf(stderr, "error: 'render_threads' not enough arguments.\n");
		return 1;
	}
	return 1;
}

int CmdRenderIterations(render_settings_t *rs)
{
	if (!ReadInteger(rs, &rs->max_iterations))
	{
		fprintf(stderr, "error: 'render_iterations' not enough arguments.\n");
		return 1;
	}
	return 1;
}


int CmdRender(render_settings_t *rs)
{
	int image_index;
	image_t *image_ptr;
	if (!ReadInteger(rs, &image_index))
	{
		fprintf(stderr, "error: 'render' not enough arguments.\n");
		return 1;
	}
	
	image_ptr = (image_t *)ListGet(rs->images_ptr, image_index).data_ptr;
	unsigned int width = image_ptr->width, height=image_ptr->height;

	vecc_t aspect = (vecc_t)height / (vecc_t)width;
	vecc_t scale_x = 2 / (vecc_t)width;
	vecc_t scale_y = -aspect * 2 / (vecc_t)height;
	
	

	clock_t start_time, end_time;
	float elapsed_secs;
	
	int max_depth = 0;

	start_time = clock();

	printf("info: generating tree.\n");
	CleanPrimitives(rs, rs->primitives_ptr);
	int tri_count = ListSize(rs->primitives_ptr);
	
	kd_tree_t *tree_ptr = GenerateTree(rs->primitives_ptr, 0, &max_depth);
	end_time = clock();
	elapsed_secs = (float)(end_time - start_time) / CLOCKS_PER_SEC;
	start_time = end_time;
	printf("info: elapsed time %f seconds.\n", elapsed_secs);
	printf("info: tree depth %i.\n", max_depth);
	printf("info: primitive count %i.\n", tri_count);
	printf("info: rendering\n");
	
	Render(
		RenderImageCallback,
		(void *)image_ptr,
		width, height,
		IdentityMatrix(),
		tree_ptr,
		rs->lights_ptr,
		rs->scene_ptr,
		rs->section_size,
		aspect,
		scale_x, scale_y,
		rs->focal_length,
		rs->samples,
		rs->max_iterations,
		rs->threads);

	end_time = clock();
	elapsed_secs = (float)(end_time - start_time) / CLOCKS_PER_SEC;
	printf("info: elapsed time %f seconds.\n", elapsed_secs);
	
	FreeTree(tree_ptr);
	ListFree(rs->materials_ptr);
	ListFree(rs->lights_ptr);
	ListFree(rs->texture2ds_ptr);
	rs->texture2ds_ptr = NULL;
	rs->materials_ptr = NULL;
	rs->lights_ptr = NULL;
	return 1;
}

static int _halton_i = 0;

int CmdHalton(render_settings_t *rs)
{
	int image_index, number;
	image_t *image_ptr;
	if (!ReadInteger(rs, &image_index) || !ReadInteger(rs, &number))
	{
		fprintf(stderr, "error: 'halton' not enough arguments.\n");
		return 1;
	}
	
	image_ptr = (image_t *)ListGet(rs->images_ptr, image_index).data_ptr;
	
	for(int i = 0; i < number; i ++)
	{
		double x = Halton(2, _halton_i);
		double y = Halton(3, _halton_i);
		_halton_i ++;
		
		int index = (
			(int)(x * image_ptr->width) +
			(int)(y * image_ptr->height) * image_ptr->width) * 4;
		
		image_ptr->data[index] = 0xFF;
		image_ptr->data[index + 1] = 0xFF;
		image_ptr->data[index + 2] = 0xFF;
		image_ptr->data[index + 3] = 0xFF;
	}
	return 1;
}

int CmdPrintPrimitives(render_settings_t *rs)
{
	printf("primitives: \n");
	PrintPrimitives(rs->primitives_ptr);
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
	if(!ReadVec3(rs, &t))
	{
		fprintf(stderr, "error: 'transform_translate' expected vec3.\n");
		return 1;
	}
	translate = NewTranslateMatrix(t);
	*rs->transform_ptr = MatrixTimesMatrix(*rs->transform_ptr, translate);
	return 1;
}

int CmdTransformScale(render_settings_t *rs)
{
	vector_t t;
	matrix_t scale;
	if(!ReadVec3(rs, &t))
	{
		fprintf(stderr, "error: 'transform_scale' expected vec3.\n");
		return 1;
	}
	scale = NewScaleMatrix(t);
	*rs->transform_ptr = MatrixTimesMatrix(*rs->transform_ptr, scale);
	return 1;
}

int CmdTransformRotate(render_settings_t *rs)
{
	vecc_t angle;
	vector_t t;
	matrix_t rotate;
	if(!ReadNumber(rs, &angle))
	{
		fprintf(stderr, "error: 'transform_rotate' expected number.\n");
	}
	if(!ReadVec3(rs, &t))
	{
		fprintf(stderr, "error: 'transform_rotate' expected vec3.\n");
		return 1;
	}
	rotate = NewRotateMatrix(angle, t);
	*rs->transform_ptr = MatrixTimesMatrix(*rs->transform_ptr, rotate);
	return 1;
}


void CommandsSetStandard(hashtable_t *table_cmds)
{
	SetCommand(table_cmds, "quit", CmdQuit);
	SetCommand(table_cmds, "vertex", CmdVertex);
	SetCommand(table_cmds, "normal", CmdNormal);
	SetCommand(table_cmds, "texco_2d", CmdTexCo2d);
	SetCommand(table_cmds, "make_face", CmdMakeFace);
	SetCommand(table_cmds, "scene_sky_color", CmdSceneSkyColor);
	SetCommand(table_cmds, "scene_set_number", CmdSceneSetNumber);
	SetCommand(table_cmds, "scene_set_integer", CmdSceneSetInteger);
	SetCommand(table_cmds, "scene_set_vector", CmdSceneSetVector);
	SetCommand(table_cmds, "cam_fov", CmdCamFov);
	SetCommand(table_cmds, "mat_index", CmdMatIndex);
	SetCommand(table_cmds, "mat_set_number", CmdMatSetNumber);
	SetCommand(table_cmds, "mat_set_integer", CmdMatSetInteger);
	SetCommand(table_cmds, "mat_set_vector", CmdMatSetVector);
	SetCommand(table_cmds, "mat_set_texture_2d", CmdMatSetTexture2d);
	SetCommand(table_cmds, "mat_shader", CmdMatShader);
	SetCommand(table_cmds, "light_position", CmdLightPosition);
	SetCommand(table_cmds, "light_set_number", CmdLightSetNumber);
	SetCommand(table_cmds, "light_set_integer", CmdLightSetInteger);
	SetCommand(table_cmds, "light_set_vector", CmdLightSetVector);
	SetCommand(table_cmds, "make_material", CmdMakeMaterial);
	SetCommand(table_cmds, "make_light", CmdMakeLight);
	SetCommand(table_cmds, "in_file", CmdInFile);
	SetCommand(table_cmds, "in_stdin", CmdInStdIn);
	SetCommand(table_cmds, "render_samples", CmdRenderSamples);
	SetCommand(table_cmds, "render_section_size", CmdRenderSectionSize);
	SetCommand(table_cmds, "render_threads", CmdRenderThreads);
	SetCommand(table_cmds, "render_iterations", CmdRenderIterations);
#ifdef INCLUDE_GUI
	SetCommand(table_cmds, "render_window", CmdRenderWindow);
#endif
	SetCommand(table_cmds, "render", CmdRender);
	SetCommand(table_cmds, "halton", CmdHalton);
	SetCommand(table_cmds, "print_primitives", CmdPrintPrimitives);
	SetCommand(table_cmds, "transform_pop", CmdTransformPop);
	SetCommand(table_cmds, "transform_push", CmdTransformPush);
	SetCommand(table_cmds, "transform_translate", CmdTransformTranslate);
	SetCommand(table_cmds, "transform_scale", CmdTransformScale);
	SetCommand(table_cmds, "transform_rotate", CmdTransformRotate);
}



