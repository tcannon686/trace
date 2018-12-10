
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "matrix.h"
#include "hashtable.h"
#include "trace.h"
#include "lodepng.h"
#include "shade.h"
#include "cmd.h"


#ifdef INCLUDE_GUI
int CmdRenderWindow(render_settings_t *rs)
{
	unsigned int width;
	unsigned int height;
	int count = fscanf(rs->input, "%u %u", &width, &height);
	if (count < 2)
	{
		fprintf(stderr, "error: 'render' not enough arguments.\n");
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
	if(rs->triangles_ptr == NULL)
	{
		rs->triangle_ptr = malloc(sizeof(tri_list_t));
		memset(rs->triangle_ptr, 0, sizeof(tri_list_t));
		rs->triangles_ptr = rs->triangle_ptr;
	}
	else
	{
		if (rs->triangle_ptr->last_ptr != NULL)
			rs->triangle_ptr->last_ptr->next_ptr = rs->triangle_ptr;
		rs->triangle_ptr->next_ptr = malloc(sizeof(tri_list_t));
		memset(rs->triangle_ptr->next_ptr, 0, sizeof(tri_list_t));
		rs->triangle_ptr->next_ptr->last_ptr = rs->triangle_ptr;
		rs->triangle_ptr = rs->triangle_ptr->next_ptr;
		rs->triangle_ptr->next_ptr = NULL;
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
	PropGetOrInsert(&(rs->material_ptr->current), key).number = value;
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
	PropGetOrInsert(&(rs->material_ptr->current), key).integer = value;
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
	PropGetOrInsert(&(rs->material_ptr->current), key).vector = value;
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

int CmdLightSetVector(render_settings_t *rs)
{
	char key[32];
	vecc_t x, y, z;

	int count = fscanf(rs->input, "%s %lf %lf %lf", key, &x, &y, &z);
	if (count == 0)
	{
		fprintf(stderr, "error: 'light_set_vector' not enough arguments.\n");
		return 1;
	}

	PropGetOrInsert(&(rs->light_ptr->current), key).vector = vec3(x, y, z);
	return 1;
}

int CmdLightSetNumber(render_settings_t *rs)
{
	char key[32];
	vecc_t x;

	int count = fscanf(rs->input, "%s %lf", key, &x);
	if (count == 0)
	{
		fprintf(stderr, "error: 'light_set_number' not enough arguments.\n");
		return 1;
	}

	PropGetOrInsert(&(rs->light_ptr->current), key).number = x;
	return 1;
}

int CmdLightSetInteger(render_settings_t *rs)
{
	char key[32];
	int x;

	int count = fscanf(rs->input, "%s %i", key, &x);
	if (count == 0)
	{
		fprintf(stderr, "error: 'light_set_integer' not enough arguments.\n");
		return 1;
	}

	PropGet(&(rs->light_ptr->current), key).integer = x;
	return 1;
}

int CmdMakeMaterial(render_settings_t *rs)
{
	if(rs->material_ptr == NULL)
	{
		rs->material_ptr = malloc(sizeof(mat_list_t));
		memset(rs->material_ptr, 0, sizeof(mat_list_t));
		rs->material_ptr->current.shader = PhongShader;
		rs->material_ptr->current.shader_data = NULL;
		rs->material_ptr->current.table = HashTableNewDefault();
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
	}
	return 1;
}

int CmdMakeLight(render_settings_t *rs)
{
	if(rs->light_ptr == NULL)
	{
		rs->light_ptr = malloc(sizeof(light_list_t));
		memset(rs->light_ptr, 0, sizeof(light_list_t));
		rs->light_ptr->current.table = HashTableNewDefault();
		rs->lights_ptr = rs->light_ptr;
	}
	else
	{
		rs->light_ptr->next_ptr = malloc(sizeof(light_list_t));
		memset(rs->light_ptr->next_ptr, 0, sizeof(light_list_t));
		rs->light_ptr->next_ptr->last_ptr = rs->light_ptr;
		rs->light_ptr = rs->light_ptr->next_ptr;
		rs->light_ptr->current.table = HashTableNewDefault();
	}
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
	CleanTriangles(&rs->triangles_ptr);
	kd_tree_t *tree_ptr = GenerateTree(rs->triangles_ptr, 0);
	end_time = clock();
	elapsed_secs = (float)(end_time - start_time) / CLOCKS_PER_SEC;
	start_time = end_time;
	printf("info: elapsed time %f seconds.\n", elapsed_secs);
	printf("info: rendering\n");
	
	Render(
		RenderImageCallback,
		(void *)&image,
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
	
	free((void *)image.data);
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


void CommandsSetStandard(hashtable_t *table_cmds)
{
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
	SetCommand(table_cmds, "light_set_number", CmdLightSetNumber);
	SetCommand(table_cmds, "light_set_integer", CmdLightSetInteger);
	SetCommand(table_cmds, "light_set_vector", CmdLightSetVector);
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
#ifdef INCLUDE_GUI
	SetCommand(table_cmds, "render_window", CmdRenderWindow);
#endif
	SetCommand(table_cmds, "render", CmdRender);
	SetCommand(table_cmds, "print_triangles", CmdPrintTriangles);
	SetCommand(table_cmds, "transform_pop", CmdTransformPop);
	SetCommand(table_cmds, "transform_push", CmdTransformPush);
	SetCommand(table_cmds, "transform_translate", CmdTransformTranslate);
	SetCommand(table_cmds, "transform_scale", CmdTransformScale);
	SetCommand(table_cmds, "transform_rotate", CmdTransformRotate);
}



