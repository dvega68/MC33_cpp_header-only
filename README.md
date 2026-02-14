### Feel free to use the MC33 c++ header-only library.

---

#### INFO:

The header-only mc33cpp.h file contains all the code for the [MC33++ library version 5.3](https://facyt-quimicomp.neocities.org/MC33_libraries.html), which is a C++ version of the code in the paper: 
Vega, D., Abache, J., Coll, D., [A Fast and Memory-Saving Marching Cubes 33 implementation with the correct interior test](http://jcgt.org/published/0008/03/01), *Journal of Computer Graphics Techniques (JCGT)*, vol. 8, no. 3, 1-18, 2019.

The MC33 library is an open source software. The distribution and use rights are under the terms of the [MIT license](https://opensource.org/licenses/MIT), described in the file "LICENSE.txt".

![FLTK example](https://repository-images.githubusercontent.com/469939412/decb05bb-c9dc-4019-96bc-11f1f6dee8c8 "Screenshot of the FLTK example")

---

#### USAGE:

1. Copy the header-only library (mc33cpp.h file) into your code folder. You need to define mc33cpp_implementation in only one of your code files before including mc33cpp.h (you can include mc33cpp.h in multiple files without issue):
	```c
	#define mc33cpp_implementation // only once in your project
	#include "mc33cpp.h"
	```
2. create a grid3D object, and read a data file (use the member functions read_grd, read_grd_binary, read_scanfiles, read_raw_file or read_dat_file):
	```c
	  grid3D G;
	  G.read_dat_file("filename.dat");

	  // or as a pointer

	  grid3D *Z = new grid3D;
	  Z->read_dat_file("filename.dat");
	```

3. create a MC33 object and assign it the grid3D
	```c
	  MC33 MC;
	  MC.set_grid3d(G); // or MC.set_grid3d(Z);
	```

4. Calculate an isosurface with the MC33 algorithm:
	```c
	  surface S;
	  MC.calculate_isosurface(S, isovalue);
	```

See the header-only library for the use of other functions.

---

#### CUSTOMIZING:

There are some options that can be modified. You have do it by editing the CUSTOMIZING section in the mc33cpp.h file:

1. To change the data type of the grid (the default value is float) define GRD_TYPE_SIZE and/or GRD_INTEGER. For example:
	```c
	//#define GRD_INTEGER
	#define GRD_TYPE_SIZE 8 // the data type is double

	#define GRD_INTEGER
	#define GRD_TYPE_SIZE 4 // the data type is unsigned int

	#define GRD_INTEGER
	#define GRD_TYPE_SIZE 2 // the data type is unsigned short int

	#define GRD_INTEGER
	#define GRD_TYPE_SIZE 1 // the data type is unsigned char
	```

2. By default the members of MC33 class are float. The member can be changed to double by defining MC33_DOUBLE_PRECISION to 1. This option is also enabled by defining GRD_TYPE_SIZE to 8. And it is disabled if GRD_TYPE_SIZE is defined to 1 or 2 when GRD_INTEGER is defined. The vertex array type of surface class is also modified when using this option.
	```c
	#define MC33_DOUBLE_PRECISION 1 // double type for MC33 class members
	```

3. If you do not use inclined grids, you can define GRD_ORTHOGONAL:
	```c
	#define GRD_ORTHOGONAL
	```

4. If you need to exchange the front and back surfaces, define MC33_NORMAL_NEG:
	```c
	#define MC33_NORMAL_NEG
	```

5. The default color of isosurfaces, can be changed:
	```c
	#define DEFAULT_SURFACE_COLOR 0xFF18A0C8// RGBA 0xAABBGGRR: red 200, green 160, blue 24
	```

---

#### COMPILING THE EXAMPLES:

In Debian terminal window, go to the FLTK_example or GLUT_example folder and write:
```sh
make -f makefiledebian.mak
```

Or in a msys2 MinGW64 Shell (Windows), write:
```sh
make -f makefileMinGW-w64.mak
```

For the FLTK example in any operating system you also can use the fltk-config script:
```sh
path/fltk-1.X.Y/fltk-config --use-gl --compile TestMC33.cpp
```

The makefiles use the -Ofast optimization option and the fltk-config script uses a lower optimization level.

In the GLUT example, the file containing the grid must be passed to the program on the command line, and no other grid files can be read from the running program. The grid file can be dragged and dropped into the executable in the Windows File Explorer. Examples of usage of the `generate_grid_from_fn` function of the grid3d class were included in this code, and are available if the grid file is not specified.

In the FLTK example, a new grid can be read from the running program, but all previous surfaces will be removed from memory. This code has keyboard shortcuts similar to the GLUT example, but using ctrl instead of alt. No usage examples for the `generate_grid_from_fn` function were included, but subgrid management was included.

---

#### GRIDS

The grid3D class has functions for building and managing subgrids that use the same data as the main grid. For example, the resolution can be reduced, and the generated isosurfaces will use less memory:
```c
  grid3D G;
  G.read_dat_file("filename.dat");
  const unsigned int *N = G.get_N();
  // build a subgrid with the half of resolution in each dimension:
  G.add_subgrid(0, 0, 0, N[0]/2, N[1]/2, N[2]/2, 2, 2, 2);

  MC33 MC;
  // The first subgrid uses one eighth of the data from the main grid
  MC.set_grid3d(G.get_subgrid(0));
```

By modifying the parameters of `add_subgrid` the grid can also be split.  
The subgrids can be deleted by using `del_subgrid(i)`, where `i` is the subgrid index.

A new function was added to the grid3D class. The `generate_grid_from_fn` function permits build a grid by using a scalar function `double fn(double x, double y, double z)`.

for example:
```c
// sphere function
double fs(double x, double y, double z) {
  const double radius = 1.0;
  const double cx = 2.0, cy = 2.0, cz = 2.0;
  x -= cx; y -= cy; z -= cz;
  return radius*radius - x*x - y*y - z*z;
}

  .
  .
  grid3D G;
  G.generate_grid_from_fn(0.5, 0.5, 0.5, // coordinates of the grid origin
                          3.5, 3.5, 3.5, // coordinates of the opposite corner
                          0.03, 0.03, 0.03, // steps
                          fs);
  MC33 MC;
  MC.set_grid3d(G);
  .
  .
```

If fn (the last argument of `generate_grid_from_fn`) is NULL, an empty grid will be created but with memory reserved for the data. The data can be filled using the `set_grid_value` function.

If you already have a data array of the same type as the data in the `grid3d` class, you can use the `set_data_pointer` function to set the internal pointers to the grid data. This avoids duplicating the data. When the `grid3d` object is destroyed, the external data will not be modified.

In this version of the library interpolation functions (trilinear and tricubic type so far) have been added. The default interpolation function is the trilinear type. It can be changed to tricubic using the `set_interpolation function`. The `interpolated_value(x, y, z)` function is used to get the interpolated value at the x, y, z position.

```c
  grid3d G;
  G.generate_grid_from_fn(0.5, 0.5, 0.5, 3.5, 3.5, 3.5, 0.03, 0.03, 0.03, fs);
  std::cout.precision(9);
  std::cout << "\nSphere: "  << fs(2.3, 2.3, 2.3); // value at a grid point
  std::cout << " " << fs(2.301, 2.302, 2.301); // point that is not on the grid
  //G.set_interpolation(1);
  std::cout << "\nLinear: " << G.interpolated_value(2.3, 2.3, 2.3);
  std::cout << " " << G.interpolated_value(2.301, 2.302, 2.301);
  G.set_interpolation(3);
  std::cout << "\n Cubic: "  << G.interpolated_value(2.3, 2.3, 2.3);
  std::cout << " " << G.interpolated_value(2.301, 2.302, 2.301);
  G.set_interpolation(0); // no valid value
  std::cout << "\n   Nan: "  << G.interpolated_value(2.3, 2.3, 2.3);
  std::cout << " " << G.interpolated_value(2.301, 2.302, 2.301);
```

For more information, see the `grid3d` class in the mc33cpp.h file.

---

#### OTHERS:

To display the surface, you can use the `draw()` and `drawdfaft()` functions of the `surface` class by defining `MC33_USE_DRAW_OPEN_GL` (in only one of the project files) before including the mc33cpp.h file in your code. These functions use the OpenGL library. Alternatively, you can implement your own drawing functions in your code after including mc33cpp.h.

To calculate the size (in bytes) of an isosurface, without calculating the isosurface, use:
```c
  size_t size = MC.size_of_isosurface(iso, nV, nT);
```
where iso is the isovalue (a `float` or `double`), nV and nT are unsigned integers that will contain the number of vertices and triangles, respectively. If you do not want to calculate nV and nT, use:
```c
  size_t size = MC.size_of_isosurface(iso);
```
See [this link](https://stackoverflow.com/questions/65066235/estimating-size-of-marching-cubes-output-geometry)

Two new funtions where added to save the surface: `surface::save_obj` and `surface::save_ply`, the first saves the surface data in a Wavefront .obj file, the other saves the data in a "Polygon File Format" (.ply) file.

---

#### ACKNOWLEDGEMENT:

Thanks to Miguel D&iacute;az (\*), for testing the library.  
And thanks again to Julien Hess (\*\*) for his helpful suggestions.

(\*) Miguel D&iacute;az e-mail: <mdiaz92@outlook.com>, IVIC doctoral student.  
(\*\*) Julien Hess e-mail: <julien.hess.ch@gmail.com>, consultant of [MATHICSE-Group](https://www.epfl.ch/labs/mathicse/)

---

See [MC33_libraries](https://facyt-quimicomp.neocities.org/MC33_libraries.html) web page.  
Mail to: <dvega@uc.edu.ve>
