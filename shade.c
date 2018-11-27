#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"
#include "trace.h"

#define THRESHOLD			0.00001

vector_t PhongShader(hit_t hit, light_list_t *lights_ptr, render_params_t *rp_ptr, int iteration, int max_iterations)
{
	vector_t diffuse_color = vec3(0, 0, 0);
	vector_t specular_color = vec3(0, 0, 0);
	vector_t V = VectorNormalize(VectorNegate(hit.position));
	light_list_t *cur_light_ptr = lights_ptr;
	if(!hit.material_ptr->shadeless)
	{
		while(cur_light_ptr != NULL)
		{
			light_t l = cur_light_ptr->current;
			vector_t Id, Is;
			vector_t Lm = VectorMinusVector(l.position, hit.position);
			vecc_t distance = VectorMagnitude(Lm);
			vecc_t attenuation = l.energy * (l.distance / (l.distance + distance));
			Lm = VectorTimesScalar(Lm, 1.0 / distance);
			vector_t N = hit.normal;
			vector_t R = VectorReflect(Lm, N);
		
			float diffuse_intensity = fmin(
				fmax(VectorDotVector(Lm, N) * attenuation, 0),
				1);
			VectorTimesScalarP(
				&Id,
				&hit.material_ptr->diffuse,
				diffuse_intensity);
			
			float dot = VectorDotVector(R, V);
			float highlight = 0;

			if(dot > 0 && diffuse_intensity > 0)
				highlight = 
					fmin(fmax(pow(dot, hit.material_ptr->shininess), 0), 1) * attenuation;
			VectorTimesScalarP(
				&Is,
				&hit.material_ptr->specular,
				highlight);
		
		
			float shadow = 1.0f;
			vector_t ray_start = l.position;
			vector_t dir = VectorNormalize(
				VectorMinusVector(hit.position, ray_start));
		
			ray_t ray = NewRay(
				ray_start,
				dir);
			hit_t shadow_hit;
			shadow_hit.t = INFINITY;
			if(RayTree(&shadow_hit, ray, hit.tree_ptr) &&
				shadow_hit.t < distance - THRESHOLD)
				shadow -= 1;
		
			VectorTimesScalarP(&Id, &Id, shadow);
			VectorTimesScalarP(&Is, &Is, shadow);
			VectorPlusVectorP(&diffuse_color, &diffuse_color, &Id);
			VectorClampP(&diffuse_color, &diffuse_color);
			VectorPlusVectorP(&specular_color, &specular_color, &Is);
			VectorClampP(&specular_color, &specular_color);
			cur_light_ptr = cur_light_ptr->next_ptr;
		}
	}
	else
	{
		diffuse_color = hit.material_ptr->diffuse;
	}
	
	if(hit.material_ptr->reflectiveness > 0 && iteration < max_iterations)
	{
		hit_t reflect_hit;
		vector_t reflect_color = rp_ptr->sky_color;
		reflect_hit.t = INFINITY;
		reflect_hit.tree_ptr = hit.tree_ptr;
		vector_t R = VectorReflect(VectorNormalize(VectorNegate(hit.ray.d)), hit.normal);
		ray_t ray = NewRay(
			VectorPlusVector(hit.position, VectorTimesScalar(R, THRESHOLD)),
			R);
		if(RayTree(&reflect_hit, ray, reflect_hit.tree_ptr))
		{
			reflect_color = reflect_hit.material_ptr->shader(
				reflect_hit, lights_ptr, rp_ptr, iteration + 1, max_iterations);
		}
		VectorClampP(&reflect_color, &reflect_color);
		diffuse_color = VectorMix(diffuse_color, reflect_color, hit.material_ptr->reflectiveness);
	}
	
	if(hit.material_ptr->alpha < 1 && iteration < max_iterations)
	{
		hit_t refract_hit;
		vector_t refract_color = rp_ptr->sky_color;
		refract_hit.t = INFINITY;
		refract_hit.tree_ptr = hit.tree_ptr;
		vector_t d;
		float n1, n2;
		if(iteration % 2 == 0)
		{
			n1 = 1.0f;
			n2 = hit.material_ptr->ior;
		}
		else
		{
			n2 = 1.0f;
			n1 = hit.material_ptr->ior;
		}
		if(VectorDotVector(hit.normal, hit.ray.d) < 0)
		{
			
			d = VectorRefract(
				VectorNormalize(hit.ray.d),
				hit.normal,
				n1,
				n2);
		}
		else
		{
			d = VectorRefract(
				VectorNormalize(hit.ray.d),
				VectorNegate(hit.normal),
				n1,
				n2);
		}
		
		ray_t ray = NewRay(
			VectorPlusVector(hit.position, VectorTimesScalar(d, THRESHOLD)),
			d);
		
		if(RayTree(&refract_hit, ray, refract_hit.tree_ptr))
		{
			refract_color = refract_hit.material_ptr->shader(
				refract_hit, lights_ptr, rp_ptr, iteration + 1, max_iterations);
		}
		VectorClampP(&refract_color, &refract_color);
		diffuse_color = VectorMix(refract_color, diffuse_color, hit.material_ptr->alpha);
	}
	
	vector_t ret;
	VectorPlusVectorP(&ret, &diffuse_color, &specular_color);
	
	return ret;
}

