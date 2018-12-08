# Python binding
"""
A geometry generator for Tom's Raytracer, a simple raytracer.

Contains classes that generate geometry: Trace, Traceable, Rotation,
Location, Scale, Cube, Sphere, Plane, Material, Light. All generators
must be attached to a Trace object, which handles rendering.
"""
import math
import subprocess
import sys

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

class Material(Traceable):
    """A material"""
    def __init__(
        self,
        trace,
        diffuse=(1, 1, 1),
        specular=(1, 1, 1),
        shininess=50,
        reflectiveness=0,
        alpha=1,
        ior=1,
        shadeless=False):
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
        self.trace = trace
        trace.materials.append(self)
        
    def toCode(self, quality):
        output = ''
        output += 'mat_set_vector diffuse ' + str(self.diffuse[0]) + ' ' + str(self.diffuse[1]) + ' ' + str(self.diffuse[2]) + '\n'
        output += 'mat_set_vector specular ' + str(self.specular[0]) + ' ' + str(self.specular[1]) + ' ' + str(self.specular[2]) + '\n'
        output += 'mat_set_number shininess ' + str(self.shininess) + '\n'
        output += 'mat_set_number reflectiveness ' + str(self.reflectiveness) + '\n'
        output += 'mat_set_number alpha ' + str(self.alpha) + '\n'
        output += 'mat_set_number ior ' + str(self.ior) + '\n'
        if self.shadeless:
            output += 'mat_set_integer shadeless 1\n'
        else:
            output += 'mat_set_integer shadeless 0\n'
        output += 'make_material\n'
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

class Mesh(Object):
    """A mesh object.

    Meshes are the same as regular objects except that they take an additional
    material parameter."""
    def __init__(self, trace, transforms = [], material=None, parent=None):
        super().__init__(trace, transforms, parent)
        self.material = material if material != None else trace.defaultMaterial
    def primitiveToCode(self, quality):
        return 'mat_index ' + str(self.trace.materials.index(self.material)) + '\n'

class Plane(Mesh):
    """A plane mesh, perpendicular to the Y-axis."""
    def __init__(self, trace, transforms = [], material = None, parent=None):
        super().__init__(trace, transforms, material, parent)
    
    def primitiveToCode(self, quality):
        ret = super().primitiveToCode(quality)
        ret +=\
"""vertex -0.5 0 -0.5
normal 0 1 0
vertex 0.5 0 -0.5
normal 0 1 0
vertex 0.5 0 0.5
normal 0 1 0
make_face
vertex 0.5 0 0.5
normal 0 1 0
vertex -0.5 0 0.5
normal 0 1 0
vertex -0.5 0 -0.5
normal 0 1 0
make_face
"""
        return ret

class Cube(Mesh):
    """A unit cube mesh."""
    primitiveCode = """vertex -0.5 0.5 -0.5
normal 0 1 0
vertex 0.5 0.5 -0.5
normal 0 1 0
vertex 0.5 0.5 0.5
normal 0 1 0
make_face
vertex 0.5 0.5 0.5
normal 0 1 0
vertex -0.5 0.5 0.5
normal 0 1 0
vertex -0.5 0.5 -0.5
normal 0 1 0
make_face
vertex -0.5 -0.5 -0.5
normal 0 -1 0
vertex 0.5 -0.5 -0.5
normal 0 -1 0
vertex 0.5 -0.5 0.5
normal 0 -1 0
make_face
vertex 0.5 -0.5 0.5
normal 0 -1 0
vertex -0.5 -0.5 0.5
normal 0 -1 0
vertex -0.5 -0.5 -0.5
normal 0 -1 0
make_face
vertex -0.5 -0.5 -0.5
normal -1 0 0
vertex -0.5 -0.5 0.5
normal -1 0 0
vertex -0.5 0.5 0.5
normal -1 0 0
make_face
vertex -0.5 -0.5 -0.5
normal -1 0 0
vertex -0.5 0.5 -0.5
normal -1 0 0
vertex -0.5 0.5 0.5
normal -1 0 0
make_face
vertex 0.5 -0.5 -0.5
normal 1 0 0
vertex 0.5 -0.5 0.5
normal 1 0 0
vertex 0.5 0.5 0.5
normal 1 0 0
make_face
vertex 0.5 -0.5 -0.5
normal 1 0 0
vertex 0.5 0.5 -0.5
normal 1 0 0
vertex 0.5 0.5 0.5
normal 1 0 0
make_face
vertex -0.5 -0.5 -0.5
normal 0 0 -1
vertex 0.5 -0.5 -0.5
normal 0 0 -1
vertex 0.5 0.5 -0.5
normal 0 0 -1
make_face
vertex -0.5 -0.5 -0.5
normal 0 0 -1
vertex 0.5 0.5 -0.5
normal 0 0 -1
vertex -0.5 0.5 -0.5
normal 0 0 -1
make_face
vertex -0.5 -0.5 0.5
normal 0 0 1
vertex 0.5 -0.5 0.5
normal 0 0 1
vertex 0.5 0.5 0.5
normal 0 0 1
make_face
vertex -0.5 -0.5 0.5
normal 0 0 1
vertex 0.5 0.5 0.5
normal 0 0 1
vertex -0.5 0.5 0.5
normal 0 0 1
make_face
"""
    def __init__(self, trace, transforms = [], material = None, parent=None):
        super().__init__(trace, transforms, material, parent)
    def primitiveToCode(self, quality):
        ret = super().primitiveToCode(quality)
        ret += Cube.primitiveCode
        return ret

class Lathe(Mesh):
    """A spun object."""
    def __init__(self, trace, transforms = [], material = None,
        x = lambda t : math.sin(t), y = lambda t : -math.cos(t),
        tmin = 0, tmax = math.pi,
        subdivideX=2, subdivideY=1,
        parent=None):
        super().__init__(trace, transforms, material, parent)
        self.x = x
        self.y = y
        self.tmin = tmin
        self.tmax = tmax
        self.subdivideX = subdivideX
        self.subdivideY = subdivideY
    
    def primitiveToCode(self, quality):
        ret = super().primitiveToCode(quality)
        subdivisionsX = quality * self.subdivideX
        subdivisionsY = quality * self.subdivideY
        stepLon = 2 * math.pi / subdivisionsX
        stepT = self.tmax / subdivisionsY
        stepY = 2.0 / (subdivisionsY)
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

                #if j != 0:
                ret += 'vertex ' + str(x00) + ' ' + str(y0) + ' ' + str(z00) + '\n'
                ret += 'normal ' + str(n00x) + ' ' + str(n0y) + ' ' + str(n00z) + '\n'
                ret += 'vertex ' + str(x01) + ' ' + str(y1) + ' ' + str(z01) + '\n'
                ret += 'normal ' + str(n01x) + ' ' + str(n1y) + ' ' + str(n01z) + '\n'
                ret += 'vertex ' + str(x10) + ' ' + str(y0) + ' ' + str(z10) + '\n'
                ret += 'normal ' + str(n10x) + ' ' + str(n0y) + ' ' + str(n10z) + '\n'
                ret += 'make_face\n'

                #if j != subdivisionsY - 1:
                ret += 'vertex ' + str(x01) + ' ' + str(y1) + ' ' + str(z01) + '\n'
                ret += 'normal ' + str(n01x) + ' ' + str(n1y) + ' ' + str(n01z) + '\n'
                ret += 'vertex ' + str(x10) + ' ' + str(y0) + ' ' + str(z10) + '\n'
                ret += 'normal ' + str(n10x) + ' ' + str(n0y) + ' ' + str(n10z) + '\n'
                ret += 'vertex ' + str(x11) + ' ' + str(y1) + ' ' + str(z11) + '\n'
                ret += 'normal ' + str(n11x) + ' ' + str(n1y) + ' ' + str(n11z) + '\n'
                ret += 'make_face\n'
        return ret

class Sphere(Lathe):
    """A unit UV sphere."""
    def __init__(self, trace, transforms = [], material = None, parent=None):
        super().__init__(trace, transforms, material,
            lambda t: math.sin(t), lambda t: -math.cos(t),
            0, math.pi, 2, 1, parent)

class Light(Object):
    """A light."""
    def __init__(self, trace, transforms = [],
        color=(1, 1, 1),
        energy=1,
        distance=30, parent=None):
        """Creates a light from the following properties:
        color, energy, and distance. By default, the light
        will attenuate to half its value at distance=30."""
        super().__init__(trace, transforms, parent)
        self.color = color
        self.energy = energy
        self.distance = distance
    
    def primitiveToCode(self, quality):
        ret = super().primitiveToCode(quality)
        ret += 'light_color ' + str(self.color[0]) + ' ' + str(self.color[1]) + ' ' + str(self.color[2]) + '\n'
        ret += 'light_energy ' + str(self.energy) + '\n'
        ret += 'light_distance ' + str(self.distance) + '\n'
        ret += 'light_position 0 0 0\n'
        ret += 'make_light\n'
        return ret


class Trace:
    """Renderer.

    Calls the actual rendering program, written in C.
    """
    def __init__(self,
        out_file = 'output.png', renderer_path = None,
        width=640, height=480, fov=70 * math.pi / 180, render_samples = 4,
        render_shadow_samples=1, render_threads = 4, render_section_size=300,
        render_iterations=10, quality=16, sky_color = (0.25, 0.25, 0.25),
        code_path=None):
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
        self.width = width
        self.height = height
        self.fov = fov
        self.render_samples = render_samples
        self.render_shadow_samples = render_shadow_samples
        self.render_threads = render_threads
        self.render_section_size = render_section_size
        self.render_iterations = render_iterations
        self.sky_color = sky_color
        self.materials = []
        self.objects = []
        self.code_path = code_path
        self.quality = quality
        self.defaultMaterial = Material(self)
    
    def trace(self):
        """Renders the scene. This will save the rendered image to out_path."""
        output = ''
        output += ''.join([material.toCode(self.quality) for material in self.materials])
        output += ''.join([obj.toCode(self.quality) for obj in self.objects])
        output += 'cam_fov ' + str(self.fov) + '\n'
        output += 'render_samples ' + str(self.render_samples) + '\n'
        output += 'render_shadow_samples ' + str(self.render_shadow_samples) + '\n'
        output += 'render_threads ' + str(self.render_threads) + '\n'
        output += 'render_section_size ' + str(self.render_section_size) + '\n'
        output += 'render_iterations ' + str(self.render_iterations) + '\n'
        output += 'sky ' + str(self.sky_color[0]) + ' ' + str(self.sky_color[1]) + ' ' + str(self.sky_color[2]) + '\n'
        output += 'out_file ' + self.out_file + '\n'
        output += 'render ' + str(self.width) + ' ' + str(self.height) + '\n'
        output += 'quit\n'
        if self.code_path != None:
            print('Saving to ' + self.code_path)
            f = open(self.code_path, 'w')
            f.write(output)
            f.close()
        return subprocess.run(
            self.renderer_path,
            input=bytes(output, 'utf-8'))



