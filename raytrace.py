#!/usr/bin/python3

# Python binding
"""
A geometry generator for Tom's Raytracer, a simple raytracer.

Contains classes that generate geometry: Trace, Traceable, Rotation,
Location, Scale, Cube, Sphere, Plane, Material, Light. All generators
must be attached to a Trace object, which handles rendering. The file also 
includes a simple mechanism for rendering movies by running ffmpeg from the 
command line.
"""
import math
import subprocess
import sys
import os
import shutil

class Traceable:
    """A basic traceable object, from which all generators inherit."""
    def toCode(self, quality):
        """Converts self into code that can be passed to the raytracer executable."""
        return ''

class Rotation(Traceable):
    """A rotation transformation"""
    def __init__(self, x=0, y=0, z=0):
        """Constructor from euler angles."""
        self.x = x
        self.y = y
        self.z = z
    def toCode(self, quality):
        return 'transform_rotate ' + str(self.x) + ', 0 1 0\n'\
            'transform_rotate ' + str(self.y) + ', 1 0 0\n'\
            'transform_rotate ' + str(self.z) + ', 0 0 1\n'

class Scale(Traceable):
    """A scale transformation."""
    def __init__(self, x=1, y=1, z=1):
        self.x = x
        self.y = y
        self.z = z
    def toCode(self, quality):
        return 'transform_scale ' + str(self.x) + ' ' + str(self.y) + ' ' + str(self.z) + '\n'

class Location(Traceable):
    """A location transformation"""
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z
    def toCode(self, quality):
        return 'transform_translate ' + str(self.x) + ' ' + str(self.y) + ' ' + str(self.z) + '\n'

class Image(Traceable):
    def __init__(
        self,
        trace, width, height):
        self.trace = trace
        self.width = width
        self.height = height
        trace.images.append(self)
    
    def toCode(self, quality):
        return 'make_image ' + str(self.width) + ' ' + str(self.height) + '\n'

class ImageFile(Image):
    def __init__(
        self,
        trace,
        path):
        super().__init__(trace, 0, 0)
        self.path = path
    
    def toCode(self, quality):
        # no super()
        return 'image_load ' + self.path + '\n'

class Texture2d(Traceable):
    def __init__(
        self,
        trace,
        image):
        self.trace = trace
        self.image = image
        trace.texture2ds.append(self)
    
    def toCode(self, quality):
        ret = ''
        ret += 'make_texture_2d ' + str(self.trace.images.index(self.image)) + '\n'
        return ret

class Material(Traceable):
    """A material."""
    def __init__(
        self,
        trace,
        diffuse=(1, 1, 1),
        specular=(1, 1, 1),
        shininess=50,
        reflectiveness=0,
        alpha=1,
        ior=1,
        ambient=(1, 1, 1),
        shadeless=False,
        tex_diffuse=None,
        enable_shadows=True):
        """Create a material from the following values:
        diffuse, specular, shininess, reflectiveness, alpha,
        ior, shadeless. By default, the material is opaque white
        with shininess 50. """
        self.diffuse = diffuse
        self.specular = specular
        self.shininess = shininess
        self.reflectiveness = reflectiveness
        self.alpha = alpha
        self.ior = ior
        self.shadeless = shadeless
        self.enable_shadows = enable_shadows
        self.tex_diffuse = tex_diffuse
        self.ambient = ambient
        
        self.trace = trace
        trace.materials.append(self)
        
    def toCode(self, quality):
        output = ''
        output += 'make_material\n'
        output += 'mat_set_vector diffuse ' + str(self.diffuse[0]) + ' ' + str(self.diffuse[1]) + ' ' + str(self.diffuse[2]) + '\n'
        output += 'mat_set_vector specular ' + str(self.specular[0]) + ' ' + str(self.specular[1]) + ' ' + str(self.specular[2]) + '\n'
        output += 'mat_set_vector ambient ' + str(self.ambient[0]) + ' ' + str(self.ambient[1]) + ' ' + str(self.ambient[2]) + '\n'
        output += 'mat_set_number shininess ' + str(self.shininess) + '\n'
        output += 'mat_set_number reflectiveness ' + str(self.reflectiveness) + '\n'
        output += 'mat_set_number alpha ' + str(self.alpha) + '\n'
        output += 'mat_set_number ior ' + str(self.ior) + '\n'
        if self.shadeless:
            output += 'mat_set_integer shadeless 1\n'
        else:
            output += 'mat_set_integer shadeless 0\n'
        if self.enable_shadows:
            output += 'mat_set_integer enable_shadows 1\n'
        else:
            output += 'mat_set_integer enable_shadows 0\n'
        if self.tex_diffuse != None:
            output += 'mat_set_texture_2d tex_diffuse ' + str(self.trace.texture2ds.index(self.tex_diffuse)) + '\n'
        return output

class Object(Traceable):
    """ An object generator, from which all objects, e.g., meshes and
    lights, inherit."""
    def __init__(self, trace, transforms = [], parent = None):
        """Creates an object.
        
        All object's have the same constructor
        pattern, (trace, transforms, ..., parent). Objects can optionally
        be attatched to another object by setting the parent argument.
        By default, parent is none. Objects with a parent inherit their transform.
        transforms is a list of transforms that will be applied to the object."""
        self.trace = trace
        self.transforms = transforms
        if parent == None:
            trace.objects.append(self)
        else:
            parent.__children.append(self)
        self.__children = []
    
    def primitiveToCode(self, quality):
        """ Converts the object's primitive data to code.

        Child classes should override this to provide functionality.
        For example, the cube class overrides it to create a unit cube.
        """
        return ''

    def toCode(self, quality):
        ret = ''
        ret += 'transform_push\n'
        for transform in self.transforms:
            ret += transform.toCode(quality)
        ret += self.primitiveToCode(quality)
        for child in self.__children:
            ret += child.toCode(quality)
        ret += 'transform_pop\n'

        return ret

class Triangle(Traceable):
    def __init__(
        self,
        trace,
        vertices=None,
        normals=None,
        texco_2ds=None):
        self.trace = trace
        self.vertices = vertices
        self.normals = normals
        self.texco_2ds = texco_2ds
    
    def toCode(self, quality):
        lines = ['make_face']
        for i in range(0, 3):
            if self.vertices:
                lines += [' '.join(['vertex'] + [str(v) for v in self.vertices[i]])]
            if self.normals:
                lines += [' '.join(['normal'] + [str(v) for v in self.normals[i]])]
            if self.texco_2ds:
                lines += [' '.join(['texco_2d'] + [str(v) for v in self.texco_2ds[i]])]
        return '\n'.join(lines)

class Quad(Traceable):
    def __init__(
        self,
        trace,
        vertices=None,
        normals=None,
        texco_2ds = None):
        self.trace = trace
        self.vertices = vertices
        self.normals = normals
        self.texco_2ds = texco_2ds
    
    def toCode(self, quality):
        lines = ['make_face']
        for i in range(0, 3):
            if self.vertices:
                lines += [' '.join(['vertex'] + [str(v) for v in self.vertices[i]])]
            if self.normals:
                lines += [' '.join(['normal'] + [str(v) for v in self.normals[i]])]
            if self.texco_2ds:
                lines += [' '.join(['texco_2d'] + [str(v) for v in self.texco_2ds[i]])]
        lines += ['make_face']
        for i in range(2, 5):
            if self.vertices:
                lines += [' '.join(['vertex'] + [str(v) for v in self.vertices[i % 4]])]
            if self.normals:
                lines += [' '.join(['normal'] + [str(v) for v in self.normals[i % 4]])]
            if self.texco_2ds:
                lines += [' '.join(['texco_2d'] + [str(v) for v in self.texco_2ds[i % 4]])]
        
        return '\n'.join(lines)

class Mesh(Object):
    """A mesh object.

    Meshes are the same as regular objects except that they take an additional
    material parameter."""
    def __init__(
        self,
        trace,
        transforms = [],
        primitives = [],
        material=None, parent=None):
        super().__init__(trace, transforms, parent)
        self.material = material if material != None else trace.default_material
        self.primitives = primitives
    def primitiveToCode(self, quality):
        ret = 'mat_index ' + str(self.trace.materials.index(self.material)) + '\n'
        ret += '\n'.join([primitive.toCode(quality) for primitive in self.primitives])
        return ret

class MeshGenerator(Mesh):
    def __init__(
        self,
        trace,
        transforms = [],
        primitives = [],
        material=None, parent=None):
        super().__init__(trace, transforms, primitives, material, parent)
    
    def genPrimitives(self, quality):
        """Set the mesh's primitives to those generated with the specified
        quality."""
        pass
    
    def primitiveToCode(self, quality):
        self.genPrimitives(quality)
        return super().primitiveToCode(quality)

class Plane(MeshGenerator):
    """A plane mesh, perpendicular to the Y-axis."""
    def __init__(
        self,
        trace,
        transforms=[],
        material=None, parent=None):
        super().__init__(trace, transforms, None, material, parent)
    
    def genPrimitives(self, quality):
        self.primitives = [
            Quad(
                self.trace,
                vertices=[
                    (-1, 0, -1),
                    (1, 0, -1),
                    (1, 0, 1),
                    (-1, 0, 1)
                ],
                normals=[
                    (0, 1, 0),
                    (0, 1, 0),
                    (0, 1, 0),
                    (0, 1, 0)
                ])]

class Cube(MeshGenerator):
    """A cube mesh."""
    def __init__(self, trace, transforms = [], material = None, parent=None):
        super().__init__(trace, transforms,
            None,
            material, parent)
    
    def genPrimitives(self, quality):
        self.primitives = [
                Quad(
                    self.trace,
                    vertices=[
                    (-1, -1, 1),
                    (1, -1, 1),
                    (1, 1, 1),
                    (-1, 1, 1)],
                    normals=[
                    (0, 0, 1),
                    (0, 0, 1),
                    (0, 0, 1),
                    (0, 0, 1)]),
                Quad(
                    self.trace,
                    vertices=[
                    (-1, -1, -1),
                    (1, -1, -1),
                    (1, 1, -1),
                    (-1, 1, -1)],
                    normals=[
                    (0, 0, -1),
                    (0, 0, -1),
                    (0, 0, -1),
                    (0, 0, -1)]),
                Quad(
                    self.trace,
                    vertices=[
                    (1, -1, -1),
                    (1, 1, -1),
                    (1, 1, 1),
                    (1, -1, 1)],
                    normals=[
                    (1, 0, 0),
                    (1, 0, 0),
                    (1, 0, 0),
                    (1, 0, 0)]),
                Quad(
                    self.trace,
                    vertices=[
                    (-1, -1, -1),
                    (-1, 1, -1),
                    (-1, 1, 1),
                    (-1, -1, 1)],
                    normals=[
                    (-1, 0, 0),
                    (-1, 0, 0),
                    (-1, 0, 0),
                    (-1, 0, 0)]),
                Quad(
                    self.trace,
                    vertices=[
                    (-1, -1, -1),
                    (1, -1, -1),
                    (1, -1, 1),
                    (-1, -1, 1)],
                    normals=[
                    (0, -1, 0),
                    (0, -1, 0),
                    (0, -1, 0),
                    (0, -1, 0)]),
                Quad(
                    self.trace,
                    vertices=[
                    (-1, 1, -1),
                    (1, 1, -1),
                    (1, 1, 1),
                    (-1, 1, 1)],
                    normals=[
                    (0, 1, 0),
                    (0, 1, 0),
                    (0, 1, 0),
                    (0, 1, 0)])]

class Lathe(MeshGenerator):
    """A spun object."""
    def __init__(self, trace, transforms = [],
        material = None,
        x = lambda t : math.sin(t), y = lambda t : -math.cos(t),
        tmin = 0, tmax = math.pi,
        subdivideX=2, subdivideY=1,
        parent=None):
        super().__init__(trace, transforms, None, material, parent)
        self.x = x
        self.y = y
        self.tmin = tmin
        self.tmax = tmax
        self.subdivideX = subdivideX
        self.subdivideY = subdivideY
    
    def genPrimitives(self, quality):
        subdivisionsX = quality * self.subdivideX
        subdivisionsY = quality * self.subdivideY
        stepLon = 2 * math.pi / subdivisionsX
        stepT = self.tmax / subdivisionsY
        stepY = 2.0 / (subdivisionsY)
        self.primitives = []
        for i in range(0, subdivisionsX):
            for j in range(0, subdivisionsY):
                thetaX = i * stepLon
                t = self.tmin + j * stepT
                y0 = self.y(t)
                y1 = self.y(t + stepT)

                r0 = self.x(t)
                r1 = self.x(t + stepT)

                x00 = r0 * math.sin(thetaX)
                z00 = r0 * math.cos(thetaX)

                x01 = r1 * math.sin(thetaX)
                z01 = r1 * math.cos(thetaX)

                x10 = r0 * math.sin(thetaX + stepLon)
                z10 = r0 * math.cos(thetaX + stepLon)

                x11 = r1 * math.sin(thetaX + stepLon)
                z11 = r1 * math.cos(thetaX + stepLon)

                d1x = (y1 - y0)
                d1y = -(r1 - r0)
                d1 = math.sqrt(d1x * d1x + d1y * d1y)
                d1x /= d1
                d1y /= d1

                d0x = (self.y(t) - self.y(t - stepT))
                d0y = -(self.x(t) - self.x(t - stepT))
                d0 = math.sqrt(d0x * d0x + d0y * d0y)
                d0x /= d0
                d0y /= d0

                d2x = (self.y(t + stepT * 2) - self.y(t + stepT))
                d2y = -(self.x(t + stepT * 2) - self.x(t + stepT))
                d2 = math.sqrt(d2x * d2x + d2y * d2y)
                d2x /= d2
                d2y /= d2
                
                n00x = (d1x + d0x) * math.sin(thetaX) / 2
                n00z = (d1x + d0x) * math.cos(thetaX) / 2
                n10x = (d1x + d0x) * math.sin(thetaX + stepLon) / 2
                n10z = (d1x + d0x) * math.cos(thetaX + stepLon) / 2

                n01x = (d1x + d2x) * math.sin(thetaX) / 2
                n01z = (d1x + d2x) * math.cos(thetaX) / 2
                n11x = (d1x + d2x) * math.sin(thetaX + stepLon) / 2
                n11z = (d1x + d2x) * math.cos(thetaX + stepLon) / 2

                n0y = (d0y + d1y) / 2
                n1y = (d1y + d2y) / 2
                self.primitives.append(
                    Quad(
                        self.trace,
                        vertices=[
                            (x00, y0, z00),
                            (x01, y1, z01),
                            (x11, y1, z11),
                            (x10, y0, z10)],
                        normals=[
                            (n00x, n0y, n00z),
                            (n01x, n1y, n01z),
                            (n11x, n1y, n11z),
                            (n10x, n0y, n10z)],
                        texco_2ds=[
                            ((thetaX / 2) / math.pi, 1 - t / math.pi),
                            ((thetaX / 2) / math.pi, 1 - (t + stepT) / math.pi),
                            (((thetaX + stepLon) / 2) / math.pi, 1 - (t + stepT) / math.pi),
                            (((thetaX + stepLon) / 2) / math.pi, 1 - (t) / math.pi)]))

                """self.primitives.append(
                    Triangle(
                        self.trace,
                        vertices=[
                            (x01, y1, z01),
                            (x10, y0, z10),
                            (x11, y1, z11)],
                        normals=[
                            (n01x, n1y, n01z),
                            (n10x, n0y, n10z),
                            (n11x, n1y, n11z)],
                        texco_2ds=[
                            ((thetaX / 2) / math.pi, 1 - (t + stepT) / math.pi),
                            (((thetaX + stepLon) / 2) / math.pi, 1 - t / math.pi),
                            (((thetaX + stepLon) / 2) / math.pi, 1 - (t + stepT) / math.pi)]))"""

class Sphere(Lathe):
    """A unit UV sphere."""
    def __init__(
        self,
        trace,
        transforms = [],
        material = None, parent=None):
        super().__init__(trace, transforms, material,
            lambda t: math.sin(t), lambda t: -math.cos(t),
            0, math.pi, 2, 1, parent)

class Light(Object):
    """A light."""
    def __init__(self, trace, transforms = [],
        color=(1, 1, 1),
        energy=1,
        distance=30, enable_shadows=True, parent=None):
        """Creates a light from the following properties:
        color, energy, and distance. By default, the light
        will attenuate to half its value at distance=30."""
        super().__init__(trace, transforms, parent)
        self.color = color
        self.energy = energy
        self.distance = distance
        self.enable_shadows = enable_shadows
    
    def primitiveToCode(self, quality):
        ret = super().primitiveToCode(quality)
        ret += 'make_light\n'
        ret += 'light_set_vector color ' + str(self.color[0]) + ' ' + str(self.color[1]) + ' ' + str(self.color[2]) + '\n'
        ret += 'light_set_number energy ' + str(self.energy) + '\n'
        ret += 'light_set_number distance ' + str(self.distance) + '\n'
        if self.enable_shadows:
            ret += 'light_set_integer enable_shadows 1\n'
        else:
            ret += 'light_set_integer enable_shadows 0\n'
        ret += 'light_position 0 0 0\n'
        return ret


class Trace:
    """Renderer.

    Calls the actual rendering program, written in C.
    """
    def __init__(self,
        out_file = 'output.png', renderer_path = None,
        width=640, height=480, fov=70 * math.pi / 180, render_samples = 4,
        shadow_samples=1, render_threads=4, render_section_size=128,
        render_iterations=4, quality=16, sky_color=(0.25, 0.25, 0.25),
        ambient_color=(0, 0, 0), ao_samples=0, ao_distance=1.0,
        render_window=False, render_image=True, code_path=None, basic=False):
        """Creates a renderer from the following properties:
        renderer_path, out_file, width, height, fov, render_samples, render_threads,
        render_section_size, render_iterations, sky_color, code_path.
        
        renderer_path is the path to the rendering program, by default, trace.exe in the 
        same directory. out_file is the output png to be written to. fov is the fov on
        the y axis, in radians. render_samples is the number of MSAA samples. render_section_size
        is the size of each section for each thread to render. render_iterations is the maximum
        depth of recursion for each reflection/refraction. sky_color is the background color
        of the image. code_path is the optional location of where to write the code generated by
        the generators."""
        if renderer_path == None:
            if sys.platform == 'win32':
                self.renderer_path = 'trace.exe'
            else:
                self.renderer_path = './trace.exe'
        else:
            self.renderer_path = renderer_path
        self.out_file = out_file
        self.render_window = render_window
        self.fov = fov
        self.render_samples = render_samples
        self.shadow_samples = shadow_samples
        self.render_threads = render_threads
        
         # Check if section size is a power of 2.
        pow2 = 2
        while pow2 < render_section_size:
        	pow2 *= 2
        if pow2 != render_section_size:
        	print("error: render_section_size must be a power of two!")
        	render_section_size = 128
        self.render_section_size = render_section_size
        
        self.render_iterations = render_iterations
        self.sky_color = sky_color
        self.ambient_color = ambient_color
        self.ao_samples = ao_samples
        self.ao_distance = ao_distance
        self.materials = []
        self.objects = []
        self.texture2ds = []
        self.images = []
        self.code_path = code_path
        self.quality = quality
        self.default_material = Material(self)
        self.render_image = render_image
        self.render_target = Image(self, width, height)
        self.basic = basic
    
    def trace(self):
        """Renders the scene. This will save the rendered image to out_path."""
        output = ''.join([image.toCode(self.quality) for image in self.images])
        output += ''.join([texture.toCode(self.quality) for texture in self.texture2ds])
        output += ''.join([material.toCode(self.quality) for material in self.materials])
        output += ''.join([obj.toCode(self.quality) for obj in self.objects])
        output += 'cam_fov ' + str(self.fov) + '\n'
        output += 'render_samples ' + ('1' if self.basic else str(self.render_samples)) + '\n'
        output += 'render_threads ' + str(self.render_threads) + '\n'
        output += 'render_section_size ' + str(self.render_section_size) + '\n'
        output += 'render_iterations ' + ('2' if self.basic else str(self.render_iterations)) + '\n'
        
        output += 'scene_sky_color '\
        	+ str(self.sky_color[0])\
        	+ ' ' + str(self.sky_color[1])\
        	+ ' ' + str(self.sky_color[2]) + '\n'
        output += 'scene_set_vector ambient '\
        	+ str(self.ambient_color[0])\
        	+ ' ' + str(self.ambient_color[1])\
        	+ ' ' + str(self.ambient_color[2]) + '\n'
        output += 'scene_set_integer shadow_samples ' + ('1' if self.basic else str(self.shadow_samples)) + '\n'
        if not self.basic:
        	output += 'scene_set_integer ao_samples ' + str(self.ao_samples) + '\n'
        elif self.ao_samples >= 1:
        	output += 'scene_set_integer ao_samples 1\n'
        else:
        	output += 'scene_set_integer ao_samples 0\n'
        output += 'scene_set_number ao_distance ' + str(self.ao_distance) + '\n'
        if self.render_window:
            output += 'render_window ' + str(self.render_target.width) + ' ' + str(self.render_target.height) + '\n'
        else:
            output += 'render ' + str(self.images.index(self.render_target)) + '\n'
            output += 'image_write ' + str(self.images.index(self.render_target)) + ' ' + str(self.out_file) + '\n'
        output += 'quit\n'
        if self.code_path != None:
            print('Saving to ' + self.code_path)
            f = open(self.code_path, 'w')
            f.write(output)
            f.close()
        return subprocess.run(
            self.renderer_path,
            input=bytes(output, 'utf-8'))

class Animator:
    """ Animators can be used with the Animation class to render animations. 
This requires ffmpeg be installed. This class is an abstract class, use 
LinearAnimator to perform animations. target is the object to change the value
of var in."""
    def __init__(self, target, var):
        self.target = target
        self.var = var
    
    def getValue(self, frame):
        return 0
    
    def move(self, frame):
        self.target.__dict__[self.var] = self.getValue(frame)

class LinearKeyframe():
    def __init__(self, time, value):
        self.time = time
        self.value = value

class LinearAnimator(Animator):
    def __init__(self, target, var, frames):
        super().__init__(target, var)
        self.frames = frames
    
    def getValue(self, frame):
        nextKeyframe = None
        lastKeyframe = None
        # Find which keyframes this is between
        for keyframe in self.frames:
            if keyframe.time <= frame:
                if lastKeyframe != None:
                    if keyframe.time > lastKeyframe.time:
                        lastKeyframe = keyframe
                else:
                    lastKeyframe = keyframe
            elif keyframe.time >= frame:
                if nextKeyframe != None:
                    if keyframe.time < nextKeyframe.time:
                        nextKeyframe = keyframe
                else:
                    nextKeyframe = keyframe
        if nextKeyframe == None or nextKeyframe.time == lastKeyframe.time:
            return lastKeyframe.value
        else:
            return lastKeyframe.value\
                + (frame - lastKeyframe.time)\
                * (nextKeyframe.value - lastKeyframe.value)\
                / (nextKeyframe.time - lastKeyframe.time)
            

class Animation:
    """A simple system for rendering animations. It works by rendering a series 
of PNG files and then running ffmpeg to combine them into a move. This does 
require ffmpeg be installed, and in the path."""
    def __init__(self, trace, animators=[], start_frame=0, end_frame=0,
        framerate=30, out_file="animation.avi"):
        self.out_file = out_file
        self.trace = trace
        self.animators = animators
        self.start_frame = start_frame
        self.end_frame = end_frame
        self.framerate = framerate
        pass
    
    def render(self):
        width = 0
        height = 0
        dirname = self.out_file.replace('.', '_')
        os.mkdir(dirname)
        for frame in range(self.start_frame, self.end_frame):
            for animator in self.animators:
                animator.move(frame)
            self.trace.out_file = dirname + '/trace' + str(frame) + '.png'
            if self.trace.render_target.width > width:
                width = self.trace.render_target.width
            if self.trace.render_target.height > height:
                height = self.trace.render_target.height
            self.trace.trace()
        subprocess.run(['ffmpeg', '-f', 'image2',
            '-framerate', str(self.framerate),
            '-i', dirname + '/trace%d.png',
            '-s', str(width) + 'x' + str(height),
            self.out_file])
        shutil.rmtree(dirname)



