#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "matrix.h"
#include "hashtable.h"
#include "linkedlist.h"
#include "trace.h"
#include "texture.h"

#define THRESHOLD			0.0001

vector_t PhongShader(
	void *shader_data,
	hit_t *hit_ptr,
	render_params_t *rp_ptr,
	int iteration, int max_iterations)
{
	vector_t diffuse_color = vec3(0, 0, 0);
	vector_t specular_color = vec3(0, 0, 0);
	vector_t V = VectorNormalize(VectorNegate(hit_ptr->position));
	light_t *cur_light_ptr;

	vector_t diffuse = PropGet(hit_ptr->material_ptr, "diffuse").vector;
	vector_t specular = PropGet(hit_ptr->material_ptr, "specular").vector;
	vecc_t shininess = PropGet(hit_ptr->material_ptr, "shininess").number;
	vecc_t reflectiveness = PropGet(hit_ptr->material_ptr, "reflectiveness").number;
	vecc_t alpha = PropGet(hit_ptr->material_ptr, "alpha").number;
	vecc_t ior = PropGet(hit_ptr->material_ptr, "ior").number;
	int shadeless = PropGet(hit_ptr->material_ptr, "shadeless").integer;
	
	hashtable_entry_t *tex_diffuse =
		HashTableGet(hit_ptr->material_ptr->table, "tex_diffuse");
	if(tex_diffuse != NULL)
	{
		texture2d_t *tex_ptr = (texture2d_t *)tex_diffuse->value.pointer;
		diffuse = Texture2dSample(tex_ptr, hit_ptr->texco);
	}
	
	if(!shadeless)
	{
		ListIterate(rp_ptr->lights_ptr, &cur_light_ptr)
		{
			light_t l = *cur_light_ptr;
			vector_t Id, Is;
			vector_t Lm = VectorMinusVector(l.position, hit_ptr->position);
			vecc_t distance = VectorMagnitude(Lm);
			vecc_t l_distance = PropGet(&l, "distance").number;
			vecc_t l_energy = PropGet(&l, "energy").number;
			vecc_t attenuation = l_energy * (l_distance / (l_distance + distance));
			Lm = VectorTimesScalar(Lm, 1.0 / distance);
			vector_t N = hit_ptr->normal;
			vector_t R = VectorReflect(Lm, N);

			
		
			float diffuse_intensity = fmin(
				fmax(VectorDotVector(Lm, N) * attenuation, 0),
				1);
			VectorTimesScalarP(
				&Id,
				&diffuse,
				diffuse_intensity);
			
			float dot = VectorDotVector(R, V);
			float highlight = 0;

			if(dot > 0 && diffuse_intensity > 0)
				highlight = 
					fmin(fmax(pow(dot, shininess), 0), 1) * attenuation;
			VectorTimesScalarP(
				&Is,
				&specular,
				highlight);
		
		
			vecc_t shadow = 1.0;
			vecc_t shadow_sample_weight = 1.0 / rp_ptr->shadow_samples;
			vector_t ray_start = VectorPlusVector(hit_ptr->position, VectorTimesScalar(hit_ptr->normal, THRESHOLD));
			vector_t dir = VectorMinusVector(l.position, ray_start);
			vector_t d;
			
			if(rp_ptr->shadow_samples > 1)
			{
				for(int i = 0; i < rp_ptr->shadow_samples; i ++)
				{
					d = dir;
					ray_t ray = NewRay(
						ray_start,
						VectorPlusVector(d, VectorTimesScalar(RandomVector(3), 0.5f)));
					hit_t shadow_hit;
					shadow_hit.t = INFINITY;
			
					if(RayTree(&shadow_hit, ray, hit_ptr->tree_ptr) &&
						shadow_hit.t < 1.0 - THRESHOLD)
						shadow -= shadow_sample_weight;
				}
			}
			else if(rp_ptr->shadow_samples == 1)
			{
				d = dir;
				ray_t ray = NewRay(
					ray_start,
					d);
				hit_t shadow_hit;
				shadow_hit.t = INFINITY;
		
				if(RayTree(&shadow_hit, ray, hit_ptr->tree_ptr) &&
					shadow_hit.t < 1.0 - THRESHOLD)
					shadow -= shadow_sample_weight;
			}
		
			VectorTimesScalarP(&Id, &Id, shadow);
			VectorTimesScalarP(&Is, &Is, shadow);
			VectorPlusVectorP(&diffuse_color, &diffuse_color, &Id);
			VectorClampP(&diffuse_color, &diffuse_color);
			VectorPlusVectorP(&specular_color, &specular_color, &Is);
			VectorClampP(&specular_color, &specular_color);
		}
	}
	else
	{
		diffuse_color = diffuse;
	}
	
	if(reflectiveness > 0 && iteration < max_iterations)
	{
		hit_t reflect_hit;
		vector_t reflect_color = rp_ptr->sky_color;
		reflect_hit.t = INFINITY;
		reflect_hit.tree_ptr = hit_ptr->tree_ptr;
		vector_t R = VectorReflect(VectorNormalize(VectorNegate(hit_ptr->ray.d)), hit_ptr->normal);
		ray_t ray = NewRay(
			VectorPlusVector(hit_ptr->position, VectorTimesScalar(R, THRESHOLD)),
			R);
		if(RayTree(&reflect_hit, ray, reflect_hit.tree_ptr))
		{
			reflect_color = ShadeHit(reflect_hit, rp_ptr, iteration + 1);
		}
		VectorClampP(&reflect_color, &reflect_color);
		diffuse_color = VectorMix(diffuse_color, reflect_color, reflectiveness);
	}
	
	if(alpha < 1 && iteration < max_iterations)
	{
		hit_t refract_hit;
		vector_t refract_color = rp_ptr->sky_color;
		refract_hit.t = INFINITY;
		refract_hit.tree_ptr = hit_ptr->tree_ptr;
		vector_t d;
		float n1, n2;
		if(iteration % 2 == 0)
		{
			n1 = 1.0f;
			n2 = ior;
		}
		else
		{
			n2 = 1.0f;
			n1 = ior;
		}
		if(VectorDotVector(hit_ptr->normal, hit_ptr->ray.d) < 0)
		{
			
			d = VectorRefract(
				VectorNormalize(hit_ptr->ray.d),
				hit_ptr->normal,
				n1,
				n2);
		}
		else
		{
			d = VectorRefract(
				VectorNormalize(hit_ptr->ray.d),
				VectorNegate(hit_ptr->normal),
				n1,
				n2);
		}
		
		ray_t ray = NewRay(
			VectorPlusVector(hit_ptr->position, VectorTimesScalar(d, THRESHOLD)),
			d);
		
		if(RayTree(&refract_hit, ray, refract_hit.tree_ptr))
		{
			refract_color = ShadeHit(refract_hit, rp_ptr, iteration + 1);
		}
		VectorClampP(&refract_color, &refract_color);
		diffuse_color = VectorMix(refract_color, diffuse_color, alpha);
	}
	
	vector_t ret;
	VectorPlusVectorP(&ret, &diffuse_color, &specular_color);
	
	return ret;
}

