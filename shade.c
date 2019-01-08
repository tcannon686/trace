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
	int sample, int max_samples,
	int iteration, int max_iterations)
{
	vector_t diffuse_color = vec3(0, 0, 0);
	vector_t specular_color = vec3(0, 0, 0);
	vector_t ambient_color;
	vector_t V = VectorNormalize(VectorNegate(hit_ptr->position));
	light_t *cur_light_ptr;
	
	vector_t ambient = PropGetOrDefault(vector, rp_ptr->scene_ptr, "ambient", vec4(0, 0, 0, 0));
	int ao_samples = PropGetOrDefault(integer, rp_ptr->scene_ptr, "ao_samples", 0);
	vecc_t ao_distance = PropGetOrDefault(number, rp_ptr->scene_ptr, "ao_distance", 1.0);
	int shadow_samples = PropGetOrDefault(integer, rp_ptr->scene_ptr, "shadow_samples", 1);
	
	vector_t diffuse = PropGetOrDefault(vector, hit_ptr->material_ptr, "diffuse", vec3(1, 1, 1));
	vector_t m_ambient = PropGetOrDefault(vector, hit_ptr->material_ptr, "ambient", vec4(1, 1, 1, 1));
	vector_t specular = PropGetOrDefault(vector, hit_ptr->material_ptr, "specular", vec3(1, 1, 1));
	vecc_t shininess = PropGetOrDefault(number, hit_ptr->material_ptr, "shininess", 10.0);
	vecc_t reflectiveness = PropGetOrDefault(number, hit_ptr->material_ptr, "reflectiveness", 0.0);
	vecc_t alpha = PropGetOrDefault(number, hit_ptr->material_ptr, "alpha", 1.0);
	vecc_t ior = PropGetOrDefault(number, hit_ptr->material_ptr, "ior", 1.0);
	int shadeless = PropGetOrDefault(number, hit_ptr->material_ptr, "shadeless", 0);
	ambient_color = ambient;
	
	hashtable_entry_t *tex_diffuse =
		HashTableGet(hit_ptr->material_ptr->table, "tex_diffuse");
	if(tex_diffuse != NULL)
	{
		texture2d_t *tex_ptr = (texture2d_t *)tex_diffuse->value.pointer;
		diffuse = Texture2dSample(tex_ptr, hit_ptr->texco);
	}
	
	// A point slightly off the surface of the intersection.
	vector_t ray_start =
		VectorPlusVector(hit_ptr->position,
			VectorTimesScalar(hit_ptr->normal, THRESHOLD));
	
	if(!shadeless)
	{
		// Ambient occlusion
		if(ao_samples > 0 && (m_ambient.x > 0 || m_ambient.y > 0 || m_ambient.z > 0))
		{
			int occlude_count = 0;
			for(int i = 0; i < ao_samples; i ++)
			{
				vector_t d;
				// Find a random vector on the hemisphere.
				do d = RandomVector(3);
				while(VectorDotVector(d, hit_ptr->normal) < 0);
				
				ray_t ray = NewRay(
					ray_start,
					d);
				
				hit_t ambient_hit;
				ambient_hit.t = INFINITY;
				if(RayTree(&ambient_hit, ray, hit_ptr->tree_ptr)
					&& ambient_hit.t <= ao_distance)
					occlude_count ++;
			}
			ambient_color = VectorTimesScalar(ambient_color, 1.0 - (vecc_t) occlude_count / ao_samples);
		}
		ambient_color = VectorTimesVector(ambient_color, m_ambient);
		
		// For each light.
		ListIterate(rp_ptr->lights_ptr, &cur_light_ptr)
		{
			light_t l = *cur_light_ptr;
			vector_t Id, Is;
			vector_t Lm = VectorMinusVector(l.position, hit_ptr->position);
			vecc_t distance = VectorMagnitude(Lm);
			vecc_t l_distance = PropGetOrDefault(number, &l, "distance", 30.0);
			vecc_t l_energy = PropGetOrDefault(number, &l, "energy", 1.0);
			vecc_t l_soft_size = PropGetOrDefault(number, &l, "soft_size", 0.25);
			int l_enable_shadows = PropGetOrDefault(integer, &l, "enable_shadows", 1);
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
			vecc_t shadow_sample_weight = 1.0 / shadow_samples;
			vector_t dir = VectorMinusVector(l.position, ray_start);
			vector_t d;
			
			if(l_enable_shadows)
			{
				if(shadow_samples > 1 && l_soft_size > 0)
				{
					for(int i = 0; i < shadow_samples; i ++)
					{
						d = dir;
						ray_t ray = NewRay(
							ray_start,
							VectorPlusVector(
								d,
								VectorTimesScalar(RandomVector(3), l_soft_size)));
						hit_t shadow_hit;
						shadow_hit.t = INFINITY;
				
						if(RayTree(&shadow_hit, ray, hit_ptr->tree_ptr) &&
							shadow_hit.t < 1.0 - THRESHOLD)
							shadow -= shadow_sample_weight;
					}
				}
				else if(shadow_samples == 1)
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
			}
		
			VectorTimesScalarP(&Id, &Id, shadow);
			VectorTimesScalarP(&Is, &Is, shadow);
			VectorPlusVectorP(&diffuse_color, &diffuse_color, &Id);
			VectorClampP(&diffuse_color, &diffuse_color);
			VectorPlusVectorP(&specular_color, &specular_color, &Is);
			VectorClampP(&specular_color, &specular_color);
		}
		diffuse_color = VectorPlusVector(diffuse_color, ambient_color);
	}
	else
	{
		diffuse_color = diffuse;
	}
	
	if(reflectiveness > 0 && iteration < max_iterations)
	{
		hit_t reflect_hit;
		vector_t reflect_color = rp_ptr->scene_ptr->sky_color;
		reflect_hit.t = INFINITY;
		reflect_hit.tree_ptr = hit_ptr->tree_ptr;
		vector_t R = VectorReflect(
			VectorNormalize(VectorNegate(hit_ptr->ray.d)),
			hit_ptr->normal);
		ray_t ray = NewRay(
			VectorPlusVector(hit_ptr->position, VectorTimesScalar(R, THRESHOLD)),
			R);
		if(RayTree(&reflect_hit, ray, reflect_hit.tree_ptr))
		{
			reflect_color = ShadeHit(reflect_hit, rp_ptr, sample, iteration + 1);
		}
		VectorClampP(&reflect_color, &reflect_color);
		diffuse_color = VectorMix(diffuse_color, reflect_color, reflectiveness);
	}
	
	if(alpha < 1 && iteration < max_iterations)
	{
		hit_t refract_hit;
		vector_t refract_color = rp_ptr->scene_ptr->sky_color;
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
			refract_color = ShadeHit(refract_hit, rp_ptr, sample, iteration + 1);
		}
		VectorClampP(&refract_color, &refract_color);
		diffuse_color = VectorMix(refract_color, diffuse_color, alpha);
	}
	
	return VectorPlusVector(diffuse_color, specular_color);
}

