# Blender as a front end.

import bpy
from mathutils import Vector
import subprocess

renderer_path = '/home/tcannon/Projects/trace/trace.exe'

scene = bpy.context.scene

output = ''
output_path = None

camera_matrix = scene.camera.matrix_world.copy()
camera_matrix.invert_safe()

def vecstring(vector):
    return str(vector[0]) + ' ' + str(vector[1]) + ' ' + str(vector[2])

last_material = None

materials = []
for material in bpy.data.materials:
    output += 'make_material\n'
    output += 'mat_set_vector diffuse ' + vecstring(material.diffuse_color * material.diffuse_intensity) + '\n'
    output += 'mat_set_vector specular ' + vecstring(material.specular_color * material.specular_intensity) + '\n'
    output += 'mat_set_number shininess ' + str(material.specular_hardness) + '\n'
    output += 'mat_set_number reflectiveness ' + str(material.raytrace_mirror.reflect_factor) + '\n'
    output += 'mat_set_number alpha ' + str(material.alpha) + '\n'
    output += 'mat_set_number ior ' + str(material.raytrace_transparency.ior) + '\n'
    if material.use_shadeless:
        output += 'mat_set_integer shadeless 1\n'
    else:
        output += 'mat_set_integer shadeless 0\n'
    materials.append(material)

for object in scene.objects:
    if object.type == 'MESH':
        triangulate = object.modifiers.new('_triangulated', 'TRIANGULATE')
        mesh = object.to_mesh(scene, True, 'RENDER', True, False)
        object.modifiers.remove(triangulate)
        for polygon in mesh.polygons:
            fn = (camera_matrix * object.matrix_world).to_3x3() * polygon.normal
            mat_index = materials.index(mesh.materials[polygon.material_index])
            output += "make_face\n"
            if mat_index != last_material:
                output += 'mat_index ' + str(mat_index) + '\n';
                last_material = mat_index
            for vi in polygon.vertices:
                vertex = camera_matrix * (object.matrix_world * mesh.vertices[vi].co)
                if polygon.use_smooth:
                    normal = (camera_matrix * object.matrix_world).to_3x3() * mesh.vertices[vi].normal
                else:
                    normal = fn
                output += "vertex " + vecstring(vertex) + "\n"
                output += "normal " + vecstring(normal) + "\n"
        bpy.data.meshes.remove(mesh)
    elif object.type == 'LAMP':
        output += 'make_light\n'
        position = camera_matrix * (object.matrix_world * Vector((0, 0, 0)))
        output += 'light_position ' + vecstring(position) + '\n'
        output += 'light_set_vector color ' + vecstring(object.data.color) + '\n'
        output += 'light_set_number energy ' + str(object.data.energy) + '\n'
        output += 'light_set_number distance ' + str(object.data.distance) + '\n'

output += 'cam_fov ' + str(scene.camera.data.angle) + '\n'

output += 'out_file /home/tcannon/Projects/trace/output.png\n'
sky = scene.world.horizon_color
output += 'sky ' + vecstring(sky) + '\n'
output += 'render_iterations 10\n'
output += 'render_samples 4\n'
output += 'render_section_size 300\n'
output += 'render_threads 4\n'
output += 'render ' + str(scene.render.resolution_x) + ' ' + str(scene.render.resolution_y) + '\n'
output += 'quit\n'

if output_path != None:
    print('Saving to ' + output_path)
    f = open(output_path, 'w')
    f.write(output)
    f.close()

print(subprocess.run(
	renderer_path,
	input=bytes(output, 'utf-8'),
    timeout=60))

