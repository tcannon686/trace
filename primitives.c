
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
#include "primitives.h"
#include "cmd.h"
#include "halton.h"

void BoundingBoxSphere(primitive_t *self)
{
    sphere_t *sphere = (sphere_t *)self->data;
    int i;

    for(i = 0; i < 3; i ++)
    {
        self->box.a.m[i] = sphere->origin.m[i] - sphere->radius * 2;
        self->box.b.m[i] = sphere->origin.m[i] + sphere->radius * 2;
    }

    self->origin = sphere->origin;
}

int RaySphere(primitive_t *self, hit_t *hit_ptr, ray_t *ray)
{
    sphere_t *sphere = (sphere_t *)self->data;

    vecc_t lDotOMinusC;
    vecc_t lenOMinusCSquared;
    vector_t oMinusC;
    vecc_t lenD = VectorMagnitude(ray->d);
    vector_t l = VectorTimesScalar(ray->d, 1.0 / lenD);
    vecc_t discriminant;
    vecc_t d, t;
	vecc_t sqrtDiscriminant;

    oMinusC = VectorMinusVector(ray->o, sphere->origin);

    lenOMinusCSquared = VectorDotVector(oMinusC, oMinusC);
    lDotOMinusC = VectorDotVector(l, oMinusC);

    discriminant = lDotOMinusC * lDotOMinusC - (lenOMinusCSquared - sphere->radius * sphere->radius);

    // Break if there was no intersection.
    if(discriminant < 0)
        return 0;

    // Otherwise return the near solution.

	sqrtDiscriminant = sqrt(discriminant);
	if(-lDotOMinusC - sqrtDiscriminant > 0)
		d = -lDotOMinusC - sqrt(discriminant);
	else
		d = -lDotOMinusC + sqrt(discriminant);
    t = d * lenD;
    
    // Don't return if the ray is behind something.
    if(t > hit_ptr->t || t <= 0)
        return 0;
    
    hit_ptr->t = t;
    hit_ptr->material_ptr = self->material_ptr;
    hit_ptr->position = 
        VectorPlusVector(
            VectorTimesScalar(l, d),
            ray->o);

    hit_ptr->normal =
        VectorTimesScalar(
            VectorMinusVector(hit_ptr->position, sphere->origin),
            1.0 / sphere->radius);
    hit_ptr->primitive_ptr = self;
    hit_ptr->texco = vec2(0, 0);

	hit_ptr->ray = *ray;

    return 1;
}

void PrintSphere(primitive_t *self)
{
    sphere_t *sphere = (sphere_t *)self->data;

    printf("sphere:\n");
    PrintVector(sphere->origin);
    printf("\n%lf\n", sphere->radius);
}

void InitPrimitives(render_settings_t *rs)
{
    rs->sphere_class.gen_box = BoundingBoxSphere;
    rs->sphere_class.hit_test = RaySphere;
    rs->sphere_class.name = "sphere";
    rs->sphere_class.print = PrintSphere;
}

int CmdMakeSphere(render_settings_t *rs)
{
    primitive_t *primitive_ptr;
	primitive_ptr = (primitive_t *)malloc(sizeof(primitive_t));
	primitive_ptr->type = &rs->sphere_class;
	primitive_ptr->data = malloc(sizeof(sphere_t));
    primitive_ptr->material_ptr = rs->current_material_ptr;
	ListAppendPointer(rs->primitives_ptr, primitive_ptr);
	return 1;
}

int CmdSphereRadius(render_settings_t *rs)
{
    vecc_t radius;

	if (!ReadNumber(rs, &radius))
	{
		fprintf(stderr, "error: 'sphere_radius' not enough arguments\n");
		return 1;
	}

	primitive_t *primitive;
	sphere_t *sphere_ptr;
	primitive = (primitive_t *)ListGetLast(rs->primitives_ptr).data_ptr;
	if(primitive->type != &rs->sphere_class)
	{
		fprintf(stderr, "error: 'sphere_radius' last primitive is not a sphere.");
		return 1;
	}

	sphere_ptr = (sphere_t *) primitive->data;
    sphere_ptr->radius = radius;
    return 1;
}

int CmdSphereOrigin(render_settings_t *rs)
{
    vector_t origin;

	if (!ReadVec3(rs, &origin))
	{
		fprintf(stderr, "error: 'sphere_origin' not enough arguments\n");
		return 1;
	}

	primitive_t *primitive;
	sphere_t *sphere_ptr;
	primitive = (primitive_t *)ListGetLast(rs->primitives_ptr).data_ptr;
	if(primitive->type != &rs->sphere_class)
	{
		fprintf(stderr, "error: 'sphere_origin' last primitive is not a sphere.");
		return 1;
	}


	sphere_ptr = (sphere_t *) primitive->data;
	MatrixTimesVectorP(&sphere_ptr->origin, rs->transform_ptr, &origin);
    return 1;
}

void CommandsSetPrimitives(hashtable_t *table)
{
    SetCommand(table, "make_sphere", CmdMakeSphere);
    SetCommand(table, "sphere_radius", CmdSphereRadius);
    SetCommand(table, "sphere_origin", CmdSphereOrigin);
}
