
# Trace
Trace is a simple raytracer I created. Its current target is Linux and Windows, but it should be relatively easy to translate to Mac and other platforms as the only platform specific code is for threading and a very rudimentary GUI. It features a simple Python front end for creating images and videos, but the C program can be used standalone. The program does not aim to be any better than any other raytracer. It was simply a learning experience to create something moderately fast, although the whole process done on the CPU and it does not take advantage of the GPU whatsoever. It features some basic effects, including shadows, ambient occlusion, reflections, and refractions. My possible future plans for this project are to add more effects, create a script for Blender, switch to some kind of path tracing, and use OpenCL to speed up rendering.

Trace uses [LodePNG](https://lodev.org/lodepng/) for PNG loading and writing.

# Compiling
Simply open the solution in Visual Studio and compile, or run `make` from the root directory. I used MinGW, GCC, and VC++ 2017 to compile it. By default, the program does not compile with GUI support. Currently, only X is supported using the `render_window` command. To compile with the GUI components, run `make INCLUDE_GUI=true`. To build a release version, run `make RELEASE=true`.

# Python Front End
The program includes a simple Python front end that starts the C program to do the rendering. The python front end requires at least Python 3. You can import the raytrace.py file to try try it out. Demos are also included in the demo folder. You can manually run the demos by adding the trace folder to the path, or by running the `demos.py [demo name]` script. More documentation is provided in the raytrace.py and demos.py files.

# C Console Application
The C console back end application is responsible for handling the actual rendering. The application reads from stdin or a file and executes commands to take in geometry. It outputs `info: %s\n` for simple info about rendering and `error: %s\n` for errors. Below is the specification for the commands. The program can also be run with an optional argument to read from a file.

## General Commands
### quit
Exits the application.
### in_file
Load commands from an external file.
### out_file
The location of the PNG file to save to.
### render <image_index>
Render to the given image. This also effectively creates a new scene to be rendered, i.e. gets rid of all geometry, lights, textures, and materials.
### render_window <width> <height>
Currently only supported on X systems and disabled by default. To enable, run `make INCLUDE_GUI=true`. This program displays a render window. The window will display the rendered image with the specified dimensions. The user may click and drag on the window to rotate the view, and use WASD to move around. The program will continue to run, even after the window is closed.
### halton <image_index> <count>
Draws the halton sequence on the image at `image_index`. This is mostly for testing.

## Scene
### scene_sky_color <r> <g> <b>
Sets the color of the sky for the scene.
### scene_set_integer <key> <value>
Set an integer key in the current scene.
### scene_set_number <key> <value>
Set a number key in the current scene.
### scene_set_vector <key> <x> <y> <z>
Set a vector key in the current scene.

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
### mat_set_integer <key> <value>
Set an integer key in the current material.
### mat_set_number <key> <value>
Set a number key in the current material.
### mat_set_vector <key> <x> <y> <z>
Set a vector key in the current material.
### mat_set_texture <key> <texture_index>
See textures below. Sets a texture key in the current material.
### mat_shader phong
Set what shader to use for the current material. Currently, only a simple phong shader is available. Different shaders use different keys in the material. See shaders below.
### mat_index <index>
Set which material will be used for subsequent faces.
### make_material
Create a new material and set it as the current material.

## Images
Images can be written to, or read from and used as textures.
### image_load <filename>
Load an image.
### image_write <index> <filename>
Write the image at index to a PNG file.
### image_set_pixel <index> <x> <y> <r> <g> <b> <a>
Set a pixel on the given image, where `r`, `g`, `b`, and `a` are between 0.0 and 1.0 inclusive.
### make_image <width> <height>
Create an image with the specified width and height.

## Textures
Textures can be created and applied to objects.
### make_texture_2d <image_index>
Make a texture from the image at `image_index`.

## Light
### light_position <x> <y> <z>
Sets the position of the current light.
### light_set_integer <key> <value>
Set an integer key in the current light.
### light_set_number <key> <value>
Set a number key in the current light.
### light_set_vector <key> <x> <y> <z>
Set a vector key in the current light.
### make_light
Create a new light and set it as the current light.

## Camera
### cam_fov <angle>
Sets the vertical FOV of the camera in radians.

## Render Options
### sky <r> <g> <b>
Sets the sky color.
### render_samples <samples>
Sets the number of super samples per pixel.
### render_section_size <size>
Sets the size of each section for each thread to render. This should be a power of two. Non power of two sizes cause undefined results.
### render_threads <count>
Sets the number of threads to start for the render.
### render_iteration <count>
Sets the maximum number of iterations for each reflection or refraction. This number must be even for refractions to work properly.

## Shaders
A shader is a function that handles shading an object and emitting more rays. Shaders may request different material data and light data. Currently, only a simple phong shader is available.
### Phong Shader
This shader is a very simple, however it supports a few simple effects, such as shadows, reflections, and refractions. The shader uses the following material vector attributes: `diffuse`, `specular`, `ambient`. It uses these number attributes: `shininess`, `reflectiveness`, `alpha`, `ior`. It also uses the integer attribute `shadeless`, and optional texture attribute `tex_diffuse`. It uses the following light number attributes: `energy`, `distance`. It also uses the following scene vector attributes: `ambient`, `ao_samples`, `shadow_samples`.

