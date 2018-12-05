
# Trace
Trace is a simple raytracer I created. Its current target is Windows, but it should be relatively easy to translate to Linux or Mac as the only platform specific code is for threading. It features a simple Python front end for creating images, but the C program can be used standalone. The program is not any more efficient than any other raytracer, nor any better, it was simply a learning experience to create something adequately speedy. It features some basic effects, including shadows, reflections, and refractions.

Trace uses [Lode PNG](https://lodev.org/lodepng/) for PNG loading and writing.

# Compiling
Simply open the solution in Visual Studio and compile, or call make from the root directory. I used MinGW and VC++ 2017 to compile it.

# Python Front End
The program includes a simple Python front end that starts the C program to do the rendering. You can import the trace.py file to try try it out. Demos are also included, simply run python on demos.py or demo1.py, demo2.py, etc. They will output demo1.png, demo2.png, etc. More documentation is provided in the trace.py file.

# C Console Application
The C console back end application is responsible for handling the actual rendering. The application reads from stdin and executes commands to take in geometry. It outputs "info: %s\n" for simple info about rendering and "error: %s\n" for errors. Below is the specification for the commands.

## General Commands
### quit
Exits the application.
### in_file
Load commands from an external file.
### out_file
The location of the PNG file to save to.
### render <width> <height>
Render an image of width and height and save it to out_file.

## Geometry
### vertex <x> <y> <z>
Creates a vertex.
### normal <x> <y> <z>
Sets the vertex normal at the current point in the triangle.
### make_face
Creates a new triangle.
### print_triangles
A debug tool to print all of the triangles currently created.

## Transform
### transform_push
Push the current transformation matrix onto the stack.
### transform_pop
Pop the current transformation matrix off the stack.
### transform_translate <x> <y> <z>
Multiply the current transformation matrix by a translation matrix.
### transform_scale <x> <y> <z>
Multiply the current transformation matrix by a scale matrix.
### transform_rotate <angle>, <x> <y> <z>
Multiply the current transformation matrix by a rotation matrix. Takes input in angle-axis form. The angle is in radians.


## Materials
### mat_index <index>
Sets the current material for following triangles.
### mat_diffuse <r> <g> <b>
Sets the diffuse color of the current material.
### mat_specular <r> <g> <b>
Sets the specular color of the current material.
### mat_ambient <r> <g> <b>
Sets the ambient color of the current material.
### mat_shininess <shininess>
Sets the shininess for the phong specular highlight of the current material.
### mat_reflectiveness <reflectiveness>
Sets the reflectiveness, from 0.0 to 1.0 of the current material.
### mat_alpha <alpha>
Sets the reflectiveness, from 0.0 to 1.0 of the current material.
### mat_ior <ior>
Sets the IOR of the current material.
### mat_shadeless {1|0}
Sets whether the current material will be shaded or not. Input of 1 will be shadeless, input of zero will be shaded.
### make_material
Create a new material and set it as the current material.

## Light
### light_position <x> <y> <z>
Sets the position of the current light.
### light_color <r> <g> <b>
Sets the color of the current light.
### light_distance <distance>
Sets the distance of the current light. At this distance, the intensity of the light will be at half its original value.
### light_energy <energy>
Sets the energy of the current light.
### make_light
Create a new light and set it as the current light.

## Camera
### cam_fov <angle>
Sets the vertical FOV of the camera. Takes radians.

## Render Options
### sky <r> <g> <b>
Sets the sky color.
### render_samples <samples>
Sets the number of MSAA samples per pixel.
### render_section_size <size>
Sets the size of each section for each thread to render.
### render_threads <count>
Sets the number of threads to start for the render.
### render_iteration <count>
Sets the maximum number of iterations for each reflection or refraction.
