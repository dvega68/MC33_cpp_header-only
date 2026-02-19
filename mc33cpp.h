/*
	File: MC33cpp.h
	Programmed by: David Vega - dvega@uc.edu.ve
	version: 5.3
	February 2026
	This library is the C++ version of the library described in the paper:
	Vega, D., Abache, J., Coll, D., A Fast and Memory Saving Marching Cubes 33
	implementation with the correct interior test, Journal of Computer Graphics
	Techniques (JCGT), vol. 8, no. 3, 1–17, 2019.
*/

#ifndef MC33cpp_h_
#define MC33cpp_h_

/********************************** USAGE ************************************/
/*
//1. Header
#define mc33cpp_implementation // only once in the project
#include "mc33cpp.h"

//2. Read a grid file.
	grid3d G;
	G.read_dat_file("filename.dat");

//3. create a MC33 object and assign it the grid3d.
	MC33 MC;
	MC.set_grid3d(G);

//4. calculate an isosurface.
	surface S;
	MC.calculate_isosurface(S, isovalue);

*/
/********************************CUSTOMIZING**********************************/
//The following defines can be changed:
//#define GRD_INTEGER // for dataset with integer type
#define GRD_TYPE_SIZE 4 // 1, 2, 4 or 8 (8 for double, if not defined GRD_INTEGER)
#define MC33_DOUBLE_PRECISION 0 // 1 means double type for MC33 class members, used only with double or size 8 integer grid data
//#define GRD_ORTHOGONAL // If defined, the library only works with orthogonal grids.
//#define MC33_NORMAL_NEG // the front and back surfaces are exchanged.
//#define DEFAULT_SURFACE_COLOR 0xFF80FF40// RGBA 0xAABBGGRR: red 64, green 255, blue 128
/*****************************************************************************/

#define MC33_VERSION_MAJOR 5
#define MC33_VERSION_MINOR 3

#ifdef GRD_INTEGER
#if GRD_TYPE_SIZE == 4
typedef unsigned int GRD_data_type; // variable type of the grid data, by default it is float.
#elif GRD_TYPE_SIZE == 2
#undef MC33_DOUBLE_PRECISION
#define MC33_DOUBLE_PRECISION 0
typedef unsigned short int GRD_data_type;
#elif GRD_TYPE_SIZE == 1
#undef MC33_DOUBLE_PRECISION
#define MC33_DOUBLE_PRECISION 0
typedef unsigned char GRD_data_type;
#else
#error "Incorrect size of the integer data type. GRD_TYPE_SIZE permitted values: 1, 2 or 4."
#endif
#elif GRD_TYPE_SIZE == 8
#undef MC33_DOUBLE_PRECISION
#define MC33_DOUBLE_PRECISION 1
typedef double GRD_data_type;
#elif MC33_DOUBLE_PRECISION
#undef GRD_TYPE_SIZE
#define GRD_TYPE_SIZE 8
typedef double GRD_data_type;
#else
typedef float GRD_data_type;
#undef GRD_TYPE_SIZE
#define GRD_TYPE_SIZE 4
#endif

#if MC33_DOUBLE_PRECISION
typedef double MC33_real;
#else
typedef float MC33_real;
#endif

#include <vector>
#include <functional>
class MC33;

/*
The class grid3d contains a 3D matrix (F[][][]) that stores the values of a function
evaluated at points of a 3D regularly spaced grid. N[] is the number of intervals in
each dimension. The grid contains (N[2] + 1)*(N[1] + 1)*(N[0] + 1) points. L[] is the
grid size in each dimension. r0[] are the coordinates of the first grid point. d[]
is the distance between adjacent points in each dimension (can be different for each
dimension), nonortho has a value of 1 when the grid is inclined else 0. _A and A_
are the matrices that transform from inclined to orthogonal coordinates and vice
versa, respectively. If the grid is periodic (is infinitely repeated along each
dimension) the flag periodic must be different from 0.

In this library, if GRD_ORTHOGONAL is defined, then nonortho, _A and A_ can be
removed from this structure, and it only works with orthogonal grids.
*/
class grid3d {
private:
	GRD_data_type ***F;
	int x_data;
	// if x_data is negative, the data of grid are external an they cannot be modified using
	// the member functions of this class. If it is positive, the grid data is managed by the
	// functions of this class. In the subgrids x_data is equal to the inner index step.
	int periodic; // to work with periodical grids.
	double r0[3], d[3]; // grid origin coordinates, distances between grid points
	unsigned int N[3];
	float L[3]; // size of the grid axes
	char title[160];
	grid3d **subgrid;
	unsigned int nsg, maxnsg; // for subgrids
	GRD_data_type (grid3d::*interpolation)(double*) const; // pointer to interpolation function
#ifndef GRD_ORTHOGONAL
	int nonortho;
	float Ang[3]; // angles between grid axes.
	double _A[3][3], A_[3][3];
	void update_matrices(); // set _A and A_ by using Ang
#endif
	int alloc_F(); // allocates memory for the grid data
	void free_F(); // release the allocated memory
	GRD_data_type bad_value(double*) const; // always returns nan
//Interpolation functions:
	GRD_data_type trilinear(double *r) const;
	GRD_data_type tricubic(double *r) const;
public:
//Get pointers to some class members:
	const unsigned int* get_N();
	const float* get_L();
	const double* get_r0();
	const double* get_d();
	const char* get_title();
#ifndef GRD_ORTHOGONAL
	const float* get_Ang();
	const double (*get__A())[3];
	const double (*get_A_())[3];
	int isnotorthogonal(); // returns nonortho value
#endif

//Generates an orthogonal grid from a function fn(x,y,z). xi and xf are the
//limits of the interval along the x axis, yi and yf along the y axis and zi
//and zf along the z axis. dx, dy and dz are the respective step sizes.
	int generate_grid_from_fn(double xi, double yi, double zi, double xf, double yf, double zf,
		double dx, double dy, double dz, double (*fn)(double x, double y, double z));
//Get a grid point value
	GRD_data_type get_grid_value(unsigned int i, unsigned int j, unsigned int k);
//Calculate a value at position (x, y, z) by interpolating the grid values.
//This function don't work with subgrids.
	GRD_data_type interpolated_value(double x, double y, double z);
//Select the interpolation type, 1 for trilinear and 3 for tricubic.
	int set_interpolation(int i);

//Modifying the grid parameters:
	/******************************************************************
	set_grid_dimensions sets the new dimensions of a grid data. It overrides the effect of
	the other functions that modify the grid parameters. Nx, Ny and Nz are the number of
	grid points in each dimension. It returns 0 if memory allocation was successful.*/
	int set_grid_dimensions(unsigned int Nx, unsigned int Ny, unsigned int Nz);
	void set_grid_value(unsigned int i, unsigned int j, unsigned int k, GRD_data_type value); // set a grid point value
	void set_ratio_aspect(double rx, double ry, double rz); // modifies d and L
	void set_r0(double x, double y, double z); // modifies r0
#ifndef GRD_ORTHOGONAL
	void set_Ang(float angle_bc, float angle_ca, float angle_ab); // modifies Ang
#endif
	void set_title(const char *s); // copy the c style string s to title
	void delete_grid_data(); // Delete the internal grid data

	/******************************************************************
	set_data_pointer creates internal pointers that point to the external data array. data
	must be stored with the nested inner loop running from i = 0 to Nx - 1 and the outer
	loop from k = 0 to Nz - 1. The data content cannot be modified using class functions.
	It returns 0 if memory allocation of internal pointers was successful.*/
	int set_data_pointer(unsigned int Nx, unsigned int Ny, unsigned int Nz, GRD_data_type* data);

// Managing subgrids:
	/******************************************************************
	The grid can contain several subgrids that use the same data. The add_subgrid function
	make a subgrid and stores it in an internal array (the subgrids do not contain data,
	only pointers to the main grid data). The input parameters of the function are the three
	indices of the origin grid point, the number of points in each dimension and the index
	step in each dimension. The function returns 0 if the subgrid was successfully added,
	-1 for wrong input parameters and -2 for error in memory allocation. */
	int add_subgrid(unsigned int Oi, unsigned int Oj, unsigned int Ok,
					unsigned int Ni, unsigned int Nj, unsigned int Nk,
					unsigned int Si, unsigned int Sj, unsigned int Sk);
	void clear_subgrid(); // erase all subgrids
	grid3d *get_subgrid(unsigned int i); // returns a pointer to the subgrid i
	void del_subgrid(unsigned int i); // delete the subgrid i
	unsigned int subgrid_size(); // returns the number of subgrids

//Reading grid data from files:
	/******************************************************************
	The functions returns zero when succeeds. If the file could not be opened or the data
	does not match the format, the return value is -1. If a memory error occurred, the
	return value is -2. A return value of -4 means that the data read may be incomplete. */
	/******************************************************************
	read_grd reads a *.grd file from the DMol3 program.*/
	int read_grd(const char *filename);

	/** read_grd_binary reads a file with an internal binary format */
	int read_grd_binary(const char *filename);
	/******************************************************************
	read_scanfiles reads a set of files that contain a slab of res*res scan data points,
	the data points are read as unsigned short int (if order is different from 0, the
	bytes of the unsigned short are exchanged). The filename must end with a number, and
	the function reads all files with end number greater or equal to filename.
	(Some datasets: http://www.graphics.stanford.edu/data/voldata/voldata.html)*/
	int read_scanfiles(const char *filename, unsigned int res, int order);

	/******************************************************************
	read_raw_file reads a file that contains integer (8, 16 or 32 bits) or float (or double)
	data points. byte is the number of bytes of the integer (1, 2 or 4), or of the float (4
	or 8). If the data is big endian, byte must be negative (only for integers). The vector
	n[3] contains the number of points in each dimension. The size of file must be
	abs(byte)*n[0]*n[1]*n[2].*/
	int read_raw_file(const char *filename, unsigned int *n, int byte, int isfloat = 0);

	/******************************************************************
	read_dat_file reads a dat file. The function returns 0 when succeeds.
	http://www.cg.tuwien.ac.at/research/vis/datasets/*/
	int read_dat_file(const char *filename);

	/******************************************************************
	save the grid points values in a raw data file of MC33_real (float or double). The
	function returns 0 when succeeds.*/
	int save_raw_file(const char *filename);

	grid3d();
	// Copy constructor
	grid3d(const grid3d &);
	~grid3d();
friend MC33;
};

template <class T>
struct MC33_v3 {
	T v[3];
};

/* The class surface contains the data of a isosurface. The isovalue is iso, the
number of surface points is nV, and the number of triangles is nT. The vector V
contains the vertex coordinates, the vector T contains the triangle indices, The
vector N contains the normal coordinates (one normal for each vertex), and the
color vector contains the color index of each point.*/
class surface {
private:
	unsigned int nV, nT;
	std::vector<MC33_v3<unsigned int>> T;
	std::vector<MC33_v3<MC33_real>> V;
	std::vector<MC33_v3<float>> N;
	std::vector<int> color;
	MC33_real iso;
	unsigned int sflag;
public:
	union {
		void *p;
		long long ul;
		int i[2];
		short si[4];
		char c[8];
		float f[2];
		double df;
	} user; // user data
	MC33_real get_isovalue(); // returns the isovalue
	unsigned int get_num_vertices(); // gets the number of vertices
	unsigned int get_num_triangles(); // gets the number of triangles
	const unsigned int *getTriangle(unsigned int n); // gets a pointer to indices of triangle n
	const MC33_real *getVertex(unsigned int n); // gets a pointer to coordinates of vertex n
	const float *getNormal(unsigned int n); // gets a pointer to the normal vector n
	void flipNormals(); // reverses the direction of all normals
	void flipTriangles(); // toggle CW / CCW vertex order in triangles
	const unsigned char *getColor(unsigned int n); // gets a pointer to the color of vertex n
	void setColor(unsigned int n, unsigned char *pcolor);

	/******************************************************************
	Saves all the surface *S data (in binary format) to a "filename" file. The
	return value is 0 if the call succeeds, else -1.*/
	int save_bin(const char *filename);

	/******************************************************************
	Saves all the surface *S data (in plain text format) to a "filename" file.
	The return value is 0 if the call succeeds, else -1.*/
	int save_txt(const char *filename);

	/******************************************************************
	Saves the surface *S data (without the color) to Wavefront .obj file.
	The return value is 0 if the call succeeds, else -1.*/
	int save_obj(const char *filename);

	/******************************************************************
	Saves the surface *S data to Polygon File Format .ply file.
	https://paulbourke.net/dataformats/ply/
	The return value is 0 if the call succeeds, else -1.*/
	int save_ply(const char *filename, const char* author = 0, const char* object = 0);

	/******************************************************************
	Reads (from a "filename" file) the surface data stored in binary format.
	The return value is 0 if the call succeeds, else -1.*/
	int read_bin(const char *filename);

	/* Draw the surface, this function can be implemented by the user.*/
	void draw();

	/* Draw some points of the surface, this function can be implemented by the user.*/
	void drawdraft();

	/* Clear all vector data */
	void clear();

	/* Correct the data vector sizes */
	void adjustvectorlenght();

	surface();
friend MC33;
};

/* Marching cubes 33 class.
The function member set_grid3d must be called once before calculate an isosurface.
If the grid3d object is modified or you want to use another grid3d object with the
same MC33 object, this function must be called again.
The function calculate_isosurface fill the surface object with the isosurface data.
*/
class MC33 {
private:
	static int DefaultColor;
	surface *S;
	//Auxiliary grid variables
	unsigned int nx, ny, nz;
	const GRD_data_type ***F;
	MC33_real MC_O[3], MC_D[3], ca, cb;
#ifndef GRD_ORTHOGONAL
	double _A[3][3], A_[3][3];
#endif
	//Assign memory for the vertex r[3], normal (r + 3)[3]. The return value is
	//the new vertex label.
	std::function<unsigned int(MC33_real*)> store_point;

	//Other auxiliary variables
	int memoryfault;
	unsigned int di; // for subgrids, index step for inner loop
	// temporary structures that store the indexes of triangle vertices:
	unsigned int **Dx, **Dy, **Ux, **Uy, **Lz;
	MC33_real *v;
	const unsigned short int table[2310]; // Triangle pattern look up table
	//Procedures
	int face_tests(int *, int) const;
	int face_test1(int) const;
	int interior_test(int, int) const;
	unsigned int surfint(unsigned int, unsigned int, unsigned int, MC33_real *);
	void find_case(unsigned int, unsigned int, unsigned int, unsigned int);
	void case_count(unsigned int, unsigned int, unsigned int, unsigned int);
	int init_temp_isosurface();
	void free_temp_D_U();
	void clear_temp_isosurface();
public:
	// Set the color of the next isosurface
	void set_default_surface_color(unsigned char *color);
	// set the grid parameters:
	int set_grid3d(grid3d *G);
	int set_grid3d(grid3d &G);
	// Calculate the isosurface with isovalue iso and store the data in the surface Sf:
	int calculate_isosurface(surface &Sf, MC33_real iso);
	/* Return the size in bytes of an isosurface with out calculate it (nV and nT are
	the number of vertices and triangles):*/
	std::size_t size_of_isosurface(MC33_real iso, unsigned int &nV, unsigned int &nT);
	// Return the size in bytes of an isosurface with out calculate it:
	std::size_t size_of_isosurface(MC33_real iso);
	MC33();
	~MC33();
};

#ifndef DEFAULT_SURFACE_COLOR
#define DEFAULT_SURFACE_COLOR 0xff5c5c5c; // grey red 92 green 92 blue 92
#endif

#ifndef GRD_ORTHOGONAL
//c = Ab, A is a 3x3 upper triangular matrix. If t != 0, A is transposed.
template<typename T> void T_multTSA_b(const double (*A)[3], T *b, T *c, int t) {
	if (t) {
		c[2] = A[0][2]*b[0] + A[1][2]*b[1] + A[2][2]*b[2];
		c[1] = A[0][1]*b[0] + A[1][1]*b[1];
		c[0] = A[0][0]*b[0];
	} else {
		c[0] = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
		c[1] = A[1][1]*b[1] + A[1][2]*b[2];
		c[2] = A[2][2]*b[2];
	}
}
//Performs the multiplication of the matrix A and the vector b: c = Ab. If t != 0, A is transposed.
template<typename T> void T_multA_b(const double (*A)[3], T *b, T *c, int t) {
	double u,v;
	if (t) {
		u = A[0][0]*b[0] + A[1][0]*b[1] + A[2][0]*b[2];
		v = A[0][1]*b[0] + A[1][1]*b[1] + A[2][1]*b[2];
		c[2] = A[0][2]*b[0] + A[1][2]*b[1] + A[2][2]*b[2];
	} else {
		u = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
		v = A[1][0]*b[0] + A[1][1]*b[1] + A[1][2]*b[2];
		c[2] = A[2][0]*b[0] + A[2][1]*b[1] + A[2][2]*b[2];
	}
	c[0] = u;
	c[1] = v;
}

extern void (*multAbf)(const double (*)[3], MC33_real *, MC33_real *, int);
extern void (*mult_TSAbf)(const double (*)[3], MC33_real *, MC33_real *, int);
extern void (*mult_Abf)(const double (*)[3], MC33_real *, MC33_real *, int);

#endif // GRD_ORTHOGONAL

/* Define MC33_USE_DRAW_OPEN_GL before include MC33.h (only once in your project), to
compile the following surface::draw and surface::drawdraft functions. */
#ifdef MC33_USE_DRAW_OPEN_GL
#if MC33_DOUBLE_PRECISION
#define GL_MC33_real GL_DOUBLE
#else
#define GL_MC33_real GL_FLOAT
#endif
void surface::draw() {
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glVertexPointer(3, GL_MC33_real, 0, &V[0]);
	glNormalPointer(GL_FLOAT, 0, &N[0]);
	glColorPointer(4, GL_UNSIGNED_BYTE, 0, &color[0]);//rgb alpha

	glDrawElements(GL_TRIANGLES, 3*nT, GL_UNSIGNED_INT, &T[0]);

	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
}

void surface::drawdraft() {
	glDisable(GL_LIGHTING);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glVertexPointer(3, GL_MC33_real, 12*sizeof(MC33_real), &V[0]);
	glColorPointer(3, GL_UNSIGNED_BYTE, 16, &color[0]);//rgb

	glDrawArrays(GL_POINTS, 0, nV>>2);

	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
	glEnable(GL_LIGHTING);
}
#endif // MC33_USE_DRAW_OPEN_GL

#ifdef mc33cpp_implementation

#ifdef _MSC_VER
#pragma warning( push )
// disable warning when a double value is assigned to float variable
#pragma warning( disable : 4244 )
#endif

#include <cmath>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef GRD_ORTHOGONAL
void (*multAbf)(const double (*)[3], MC33_real *, MC33_real *, int);
void (*mult_TSAbf)(const double (*)[3], MC33_real *, MC33_real *, int) = T_multTSA_b;
void (*mult_Abf)(const double (*)[3], MC33_real *, MC33_real *, int) = T_multA_b;

#if GRD_TYPE_SIZE == 8
void (*multAb)(const double (*)[3], double *, double *, int) = mult_Abf;
#else
void (*multAb)(const double (*)[3], double *, double *, int) = T_multA_b;
#endif

void setIdentMat3x3d(double (*A)[3]) {
	for (double *d = A[0] + 8; --d != A[0];)
		d[0] = 0.0;
	for (int i = 0; i != 3; i++)
		A[i][i] = 1.0;
}
#endif // GRD_ORTHOGONAL

#include <algorithm>

#if defined (__SSE__) || ((defined (_M_IX86) || defined (_M_X64)) && !defined(_CHPE_ONLY_))
// https://stackoverflow.com/questions/59644197/inverse-square-root-intrinsics
// faster than 1.0f/std::sqrt, but with little accuracy.
#include <immintrin.h>
inline float invSqrt(float f) {
	__m128 temp = _mm_set_ss(f);
	temp = _mm_rsqrt_ss(temp);
	return _mm_cvtss_f32(temp);
}
#else
inline float invSqrt(float f) {
	return 1.0/sqrt(f);
}
#endif

MC33::MC33() : F(0), table{
/* INDEX
The 11 least significant bits either contain the position of the triangle pattern
in this array (cases 1, 2, 5, 8, 9, 11, and 14) or are used to calculate it. The
12th bit is used to determine whether the order of the triangle vertices (and the
normal) should be reversed. The 4 most significant bits are related to the MC33 case.

     Case              Patern position
  __________     _______________________________
  F  E  D  C  B  A  9  8  7  6  5  4  3  2  1  0
              ^
           reverse
*/
 0x0000, 0x0885, 0x0886, 0x0895, 0x0883, 0x1816, 0x089D, 0x0943,//007
 0x0884, 0x0897, 0x1814, 0x0916, 0x0891, 0x094C, 0x091F, 0x048F,//015
 0x0882, 0x089B, 0x1808, 0x0934, 0x2803, 0x3817, 0x3814, 0x0525,//023
 0x180E, 0x0928, 0x4802, 0x049D, 0x3815, 0x0541, 0x6004, 0x0110,//031
 0x0881, 0x180A, 0x0899, 0x0922, 0x1806, 0x4800, 0x0913, 0x0499,//039
 0x2802, 0x380E, 0x3811, 0x053D, 0x380D, 0x6002, 0x0521, 0x010D,//047
 0x088B, 0x0946, 0x0937, 0x0493, 0x3813, 0x6003, 0x0545, 0x012E,//055
 0x380F, 0x0531, 0x600A, 0x011C, 0x5001, 0x3001, 0x3009, 0x0087,//063
 0x0880, 0x2801, 0x1804, 0x3807, 0x088F, 0x380B, 0x092B, 0x0539,//071
 0x1802, 0x3806, 0x4806, 0x6001, 0x0925, 0x051D, 0x0495, 0x010A,//079
 0x1812, 0x380A, 0x4805, 0x6009, 0x3816, 0x5002, 0x6008, 0x3010,//087
 0x4801, 0x6000, 0x7000, 0x4003, 0x6006, 0x3004, 0x4007, 0x1010,//095
 0x0889, 0x3808, 0x093A, 0x052D, 0x0949, 0x6005, 0x0491, 0x0131,//103
 0x380C, 0x5000, 0x600B, 0x3012, 0x0549, 0x3002, 0x0119, 0x008D,//111
 0x0907, 0x0535, 0x04A1, 0x013D, 0x0529, 0x3005, 0x0140, 0x0093,//119
 0x6007, 0x3000, 0x4004, 0x1000, 0x3003, 0x2000, 0x100C, 0x007f,//127

/*
  Vertices:            Edges:               Faces:
    3 ___________2        _____2______         ____________
   /|           /|      /|           /|      /|           /|
  / |          / |     B |          A |     / |    2     / |
7/___________6/  |    /_____6_____ /  |    /___________ /  |
|   |        |   |   |   3        |   1   |   |     4  |   |
|   |        |   |   |   |        |   |   | 3 |        | 1 |     z
|   0________|___1   |   |_____0__|___|   |   |_5______|___|     |
|  /         |  /    7  /         5  /    |  /         |  /      |____y
| /          | /     | 8          | 9     | /      0   | /      /
4/___________5/      |/_____4_____|/      |/___________|/      x

*/
/*
Vertices order in triangles:
       0     1-----2
      / \     \ x /   o: Front surface
     / o \     \ /    x: Back surface
    1-----2     0
*/
/* TRIANGLE PATERNS
Each short integer corresponds to one triangle. The most significant hexadecimal
digit is set to 0 for the last triangle in the pattern. The other 3 digits are
the cube edges where the triangle vertices are located.
*/
/*position*index#vertices*/
// Case 1 (128)
/* 0*127#0*/0x0380,
/* 1*064#1*/0x0109,
/* 2*032#2*/0x021A,
/* 3*016#3*/0x0B32,
/* 4*004#5*/0x0945,
/* 5*008#4*/0x0748,
/* 6*001#7*/0x07B6,
/* 7*002#6*/0x06A5,

// Case 2 (136)
/* 0*063#01*/0x1189,0x0138,
/* 1*096#12*/0x129A,0x0092,
/* 2*048#23*/0x1B3A,0x031A,
/* 3*111#03*/0x12B0,0x0B80,
/* 4*068#15*/0x1045,0x0105,
/* 5*012#45*/0x1975,0x0987,
/* 6*119#04*/0x1340,0x0374,
/* 7*003#67*/0x1BA5,0x07B5,
/* 8*009#47*/0x1486,0x0B68,
/* 9*034#26*/0x1615,0x0216,
/*10*017#37*/0x1726,0x0732,
/*11*006#56*/0x146A,0x094A,

// Case 3.1 (160)
/* 0*123#05*/0x1945,0x0038,
/* 1*072#14*/0x1109,0x0748,
/* 2*066#16*/0x16A5,0x0109,
/* 3*036#25*/0x1945,0x021A,
/* 4*018#36*/0x16A5,0x032B,
/* 5*033#27*/0x17B6,0x021A,
/* 6*126#07*/0x17B6,0x0038,
/* 7*024#34*/0x1B32,0x0748,
/* 8*095#02*/0x121A,0x0038,
/* 9*080#13*/0x1B32,0x0109,
/*10*010#46*/0x16A5,0x0874,
/*11*005#57*/0x1945,0x07B6,

// Case 3.2 (184)
/* 0*123#05*/0x1905,0x1035,0x1453,0x0843,
/* 1*072#14*/0x1974,0x1917,0x1871,0x0081,
/* 2*066#16*/0x1605,0x1950,0x116A,0x0106,
/* 3*036#25*/0x1A45,0x1942,0x124A,0x0192,
/* 4*018#36*/0x13A5,0x1B56,0x132A,0x0B35,
/* 5*033#27*/0x11A6,0x17B2,0x1217,0x0671,
/* 6*126#07*/0x1786,0x1806,0x1B60,0x03B0,
/* 7*024#34*/0x1834,0x1324,0x1742,0x0B72,
/* 8*095#02*/0x123A,0x138A,0x11A8,0x0018,
/* 9*080#13*/0x1129,0x12B9,0x109B,0x030B,
/*10*010#46*/0x14A5,0x1A86,0x148A,0x0768,
/*11*005#57*/0x1B65,0x1794,0x17B9,0x059B,

// Case 4.1.1 (232)
/* 0*125#06*/0x16A5,0x0380,
/* 1*065#17*/0x17B6,0x0109,
/* 2*040#24*/0x121A,0x0748,
/* 3*020#35*/0x1945,0x0B32,

//The numbers in parentheses are the diagonal (interior test)
// Case 4.1.2 (240)
/* 0*(06)*125#06*/0x10A5,0x1805,0x1863,0x1685,0x16A3,0x03A0,
/* 1*(17)*065#17*/0x1796,0x11B6,0x1169,0x1970,0x17B0,0x00B1,
/* 2*(24)*040#24*/0x174A,0x17A2,0x1872,0x1A41,0x1481,0x0182,
/* 3*(35)*020#35*/0x12B5,0x1B34,0x1943,0x15B4,0x1592,0x0293,

// Case 5 (264)
/* 0*112#123*/0x1B9A,0x1930,0x09B3,
/* 1*079#023*/0x180A,0x101A,0x0B8A,
/* 2*047#013*/0x189B,0x1B12,0x0B91,
/* 3*031#012*/0x128A,0x1382,0x09A8,
/* 4*038#256*/0x1246,0x1419,0x0421,
/* 5*011#467*/0x18A5,0x1485,0x0BA8,
/* 6*110#037*/0x1026,0x1067,0x0078,
/* 7*059#015*/0x1845,0x1538,0x0135,
/* 8*014#456*/0x176A,0x1A87,0x098A,
/* 9*035#267*/0x1715,0x11B2,0x017B,
/*10*076#145*/0x1175,0x1708,0x0710,
/*11*025#347*/0x1426,0x1832,0x0248,
/*12*070#156*/0x106A,0x1460,0x010A,
/*13*055#014*/0x1149,0x1174,0x0371,
/*14*103#034*/0x142B,0x174B,0x0024,
/*15*019#367*/0x12A5,0x1532,0x0735,
/*16*050#326*/0x1635,0x136B,0x0153,
/*17*098#126*/0x1695,0x1609,0x0206,
/*18*115#045*/0x1935,0x1390,0x0753,
/*19*118#047*/0x13B6,0x1360,0x0406,
/*20*007#576*/0x19BA,0x17B4,0x0B94,
/*21*049#237*/0x17A6,0x171A,0x0317,
/*22*100#125*/0x1A25,0x1245,0x0042,
/*23*013#457*/0x1965,0x19B6,0x08B9,

// Case 6.1.1 (336)
/* 0*121#056*/0x146A,0x1A94,0x0380,
/* 1*061#016*/0x16A5,0x1189,0x0138,
/* 2*109#036*/0x16A5,0x180B,0x002B,
/* 3*124#067*/0x1BA5,0x157B,0x0038,
/* 4*093#026*/0x1615,0x1621,0x0803,
/* 5*117#046*/0x16A5,0x1034,0x0374,
/* 6*073#147*/0x18B6,0x1648,0x0109,
/* 7*067#167*/0x1BA5,0x17B5,0x0109,
/* 8*097#127*/0x129A,0x17B6,0x0092,
/* 9*062#017*/0x17B6,0x1189,0x0138,
/*10*081#137*/0x1726,0x1732,0x0109,
/*11*069#157*/0x1045,0x17B6,0x0105,
/*12*104#124*/0x129A,0x1209,0x0748,
/*13*044#245*/0x1975,0x121A,0x0879,
/*14*041#247*/0x1486,0x121A,0x08B6,
/*15*056#234*/0x131A,0x1B3A,0x0874,
/*16*087#024*/0x121A,0x1374,0x0403,
/*17*042#246*/0x1615,0x1216,0x0874,
/*18*107#035*/0x1945,0x12B0,0x0B80,
/*19*052#235*/0x1945,0x1B3A,0x031A,
/*20*022#356*/0x146A,0x194A,0x032B,
/*21*028#345*/0x1975,0x1987,0x02B3,
/*22*084#135*/0x1045,0x1510,0x0B32,
/*23*021#357*/0x1945,0x1726,0x0732,

// Case 6.1.2 (408)
/* 0*(06)*121#056*/0x136A,0x190A,0x1094,0x1804,0x1684,0x1386,0x03A0,
/* 1*(06)*061#016*/0x1895,0x11A5,0x136A,0x1591,0x13A1,0x1863,0x0856,
/* 2*(06)*109#036*/0x10A5,0x1AB6,0x102A,0x1A2B,0x186B,0x1568,0x0058,
/* 3*(06)*124#067*/0x10A5,0x1785,0x187B,0x138B,0x1A3B,0x103A,0x0058,
/* 4*(06)*093#026*/0x1685,0x1236,0x1321,0x1031,0x1501,0x1805,0x0863,
/* 5*(06)*117#046*/0x10A5,0x1456,0x1376,0x1674,0x1054,0x13A0,0x0A36,
/* 6*(17)*073#147*/0x1496,0x1948,0x1098,0x1B08,0x110B,0x161B,0x0169,
/* 7*(17)*067#167*/0x1795,0x11A5,0x11BA,0x1915,0x1097,0x1B07,0x00B1,
/* 8*(17)*097#127*/0x19A6,0x16A2,0x1B62,0x10B2,0x17B0,0x1970,0x0796,
/* 9*(17)*062#017*/0x1796,0x113B,0x1B38,0x17B8,0x1978,0x1169,0x01B6,
/*10*(17)*081#137*/0x1796,0x1730,0x1032,0x1102,0x1612,0x1916,0x0970,
/*11*(17)*069#157*/0x1165,0x1745,0x1756,0x1704,0x1B61,0x10B1,0x0B07,
/*12*(24)*104#124*/0x149A,0x1208,0x1809,0x1489,0x174A,0x127A,0x0728,
/*13*(24)*044#245*/0x19A5,0x175A,0x11A9,0x1819,0x1218,0x1728,0x027A,
/*14*(24)*041#247*/0x14A6,0x18B2,0x12B6,0x1A26,0x11A4,0x1814,0x0182,
/*15*(24)*056#234*/0x141A,0x1B7A,0x17B3,0x1873,0x1183,0x1481,0x04A7,
/*16*(24)*087#024*/0x141A,0x1014,0x1103,0x1213,0x1723,0x1A27,0x04A7,
/*17*(24)*042#246*/0x1645,0x1415,0x1746,0x1276,0x1872,0x1182,0x0814,
/*18*(35)*107#035*/0x1925,0x1B84,0x1480,0x1940,0x1290,0x1B52,0x05B4,
/*19*(35)*052#235*/0x19A5,0x1319,0x191A,0x1B5A,0x145B,0x134B,0x0439,
/*20*(35)*022#356*/0x1B6A,0x1B46,0x12BA,0x192A,0x1329,0x1439,0x034B,
/*21*(35)*028#345*/0x12B5,0x1398,0x1387,0x1B37,0x15B7,0x1925,0x0293,
/*22*(35)*084#135*/0x1B45,0x1125,0x1210,0x1320,0x1430,0x1B34,0x0B52,
/*23*(35)*021#357*/0x1265,0x1574,0x1567,0x1347,0x1943,0x1293,0x0925,

// Case 6.2 (576)
/* 0*121#056*/0x136A,0x190A,0x13A0,0x1684,0x0386,
/* 1*061#016*/0x1685,0x136A,0x1895,0x113A,0x0863,
/* 2*109#036*/0x10A5,0x1856,0x102A,0x186B,0x0058,
/* 3*124#067*/0x10A5,0x1785,0x1058,0x1A3B,0x003A,
/* 4*093#026*/0x1685,0x1236,0x1863,0x1501,0x0805,
/* 5*117#046*/0x10A5,0x1A36,0x1054,0x1376,0x03A0,
/* 6*073#147*/0x1496,0x1169,0x1B08,0x110B,0x061B,
/* 7*067#167*/0x1795,0x11BA,0x10B1,0x1097,0x0B07,
/* 8*097#127*/0x19A6,0x1796,0x10B2,0x17B0,0x0970,
/* 9*062#017*/0x1796,0x113B,0x161B,0x1978,0x0169,
/*10*081#137*/0x1796,0x1126,0x1730,0x1970,0x0916,
/*11*069#157*/0x1165,0x1704,0x1B07,0x1B61,0x00B1,
/*12*104#124*/0x149A,0x1208,0x1728,0x174A,0x027A,
/*13*044#245*/0x1A75,0x127A,0x1819,0x1218,0x0728,
/*14*041#247*/0x14A6,0x18B2,0x1182,0x11A4,0x0814,
/*15*056#234*/0x174A,0x1B7A,0x1183,0x1481,0x0A41,
/*16*087#024*/0x141A,0x1014,0x1723,0x1A27,0x04A7,
/*17*042#246*/0x1415,0x1276,0x1814,0x1872,0x0182,
/*18*107#035*/0x1925,0x1B84,0x15B4,0x1290,0x0B52,
/*19*052#235*/0x1AB5,0x1319,0x1439,0x145B,0x034B,
/*20*022#356*/0x1B46,0x192A,0x134B,0x1329,0x0439,
/*21*028#345*/0x12B5,0x1398,0x1293,0x15B7,0x0925,
/*22*084#135*/0x12B5,0x1125,0x1430,0x1B34,0x05B4,
/*23*021#357*/0x1265,0x1925,0x1734,0x1943,0x0293,

// Case 7.1 (696)
/* 0*037#257*/0x1945,0x121A,0x07B6,
/* 1*088#134*/0x12B3,0x1874,0x0109,
/* 2*026#346*/0x16A5,0x1B32,0x0874,
/* 3*091#025*/0x1945,0x121A,0x0038,
/* 4*122#057*/0x1945,0x17B6,0x0038,
/* 5*082#136*/0x16A5,0x1B32,0x0109,
/* 6*074#146*/0x16A5,0x1109,0x0748,
/* 7*094#027*/0x121A,0x17B6,0x0380,

//The characters inside of the square bracket are face test results
// Case 7.2 (720)
/* 0*037#257*[.--..+]*/0x1B65,0x19B5,0x121A,0x1794,0x0B97,
/* 1*088#134*[-..-+.]*/0x1B92,0x1874,0x1B30,0x19B0,0x0912,
/* 2*026#346*[..--.+]*/0x14A5,0x1A86,0x1B32,0x1876,0x08A4,
/* 3*091#025*[--..+.]*/0x1945,0x138A,0x1801,0x1A81,0x0A23,
/* 4*122#057*[-..-.+]*/0x1B65,0x1380,0x1794,0x1B97,0x09B5,
/* 5*082#136*[.--.+.]*/0x16A5,0x1129,0x1B92,0x1B30,0x09B0,
/* 6*074#146*[--...+]*/0x14A5,0x1876,0x1109,0x18A4,0x0A86,
/* 7*094#027*[..--+.]*/0x17B6,0x123A,0x18A3,0x1801,0x0A81,
//(760)
/* 0*037#257*[.+-..-]*/0x1A45,0x17B6,0x124A,0x1219,0x0429,
/* 1*088#134*[-..+-.]*/0x1109,0x1483,0x1243,0x12B7,0x0427,
/* 2*026#346*[..-+.-]*/0x16A5,0x1274,0x12B7,0x1483,0x0243,
/* 3*091#025*[-+..-.]*/0x1A45,0x1038,0x1219,0x1429,0x024A,
/* 4*122#057*[-..+.-]*/0x1945,0x1806,0x1678,0x103B,0x060B,
/* 5*082#136*[.+-.-.]*/0x1605,0x1095,0x116A,0x132B,0x0061,
/* 6*074#146*[-+...-]*/0x1605,0x1095,0x116A,0x1748,0x0061,
/* 7*094#027*[..-+-.]*/0x1786,0x121A,0x103B,0x160B,0x0068,
//(800)
/* 0*037#257*[.-+..-]*/0x1945,0x11A6,0x1716,0x17B2,0x0172,
/* 1*088#134*[+..--.]*/0x132B,0x1749,0x1108,0x1718,0x0179,
/* 2*026#346*[..+-.-]*/0x13A5,0x1B56,0x1748,0x135B,0x032A,
/* 3*091#025*[+-..-.]*/0x1905,0x121A,0x1350,0x1384,0x0534,
/* 4*122#057*[+..-.-]*/0x1905,0x1534,0x17B6,0x1384,0x0350,
/* 5*082#136*[.-+.-.]*/0x13A5,0x1B56,0x1109,0x132A,0x035B,
/* 6*074#146*[+-...-]*/0x16A5,0x1749,0x1179,0x1108,0x0718,
/* 7*094#027*[..+--.]*/0x11A6,0x1380,0x17B2,0x1172,0x0716,

// Case 7.3 (840)
/* 0*037#257*[.+-..+]*/0x1C65,0x1C5A,0x17C4,0x1C7B,0x1C94,0x1C19,0x1C21,0x1CA2,0x0CB6,
/* 1*088#134*[-..++.]*/0x1C74,0x1CB7,0x1C2B,0x11C9,0x1C12,0x1C09,0x1C30,0x1C83,0x0C48,
/* 2*026#346*[..-+.+]*/0x1CA5,0x1C6A,0x1C32,0x1C83,0x1C48,0x1C54,0x1C76,0x1CB7,0x0C2B,
/* 3*091#025*[-+..+.]*/0x1AC5,0x1C38,0x1C23,0x1CA2,0x1C45,0x1C94,0x1C19,0x1C01,0x0C80,
/* 4*122#057*[-..+.+]*/0x1C65,0x19C5,0x1CB6,0x1C3B,0x1C03,0x1C80,0x17C4,0x1C78,0x0C94,
/* 5*082#136*[.+-.+.]*/0x16C5,0x1C6A,0x1C95,0x1C09,0x1C30,0x1CB3,0x1C2B,0x1C12,0x0CA1,
/* 6*074#146*[-+...+]*/0x1C95,0x1C6A,0x1C10,0x1CA1,0x1C48,0x1C76,0x1C87,0x1C54,0x0C09,
/* 7*094#027*[..-++.]*/0x1C1A,0x17C6,0x1C01,0x1C80,0x1C78,0x1CB6,0x1C3B,0x1C23,0x0CA2,
//(912)
/* 0*037#257*[.-+..+]*/0x1C65,0x1CA6,0x1C59,0x11C2,0x1CB2,0x17C4,0x1C7B,0x1C94,0x0C1A,
/* 1*088#134*[+..-+.]*/0x1C2B,0x11C9,0x1C12,0x1C49,0x1C74,0x1C87,0x1C08,0x1C30,0x0CB3,
/* 2*026#346*[..+-.+]*/0x1CA5,0x1C54,0x1C76,0x1C48,0x1C2A,0x1C32,0x1CB3,0x1C6B,0x0C87,
/* 3*091#025*[+-..+.]*/0x1C45,0x12CA,0x1C84,0x1C38,0x1C23,0x1C59,0x1C1A,0x1C01,0x0C90,
/* 4*122#057*[+..-.+]*/0x1C65,0x1C03,0x1C90,0x1C59,0x1CB6,0x17C4,0x1C7B,0x1C84,0x0C38,
/* 5*082#136*[.-+.+.]*/0x1CA5,0x1C56,0x1C09,0x1C30,0x1CB3,0x1C6B,0x1C2A,0x1C12,0x0C91,
/* 6*074#146*[+-...+]*/0x1CA5,0x1C6A,0x1C54,0x1C76,0x1C87,0x1C08,0x11C9,0x1C10,0x0C49,
/* 7*094#027*[..+-+.]*/0x1CA6,0x17C6,0x1C1A,0x1C01,0x1C80,0x1C38,0x1C23,0x1CB2,0x0C7B,
//(984)
/* 0*037#257*[.++..-]*/0x1AC5,0x1CA6,0x1C94,0x1C19,0x1C21,0x1CB2,0x1C7B,0x1C67,0x0C45,
/* 1*088#134*[+..+-.]*/0x1C2B,0x11C9,0x1C49,0x1C74,0x1CB7,0x1C32,0x1C83,0x1C08,0x0C10,
/* 2*026#346*[..++.-]*/0x1CA5,0x1C56,0x1C2A,0x1C32,0x1C83,0x1C48,0x1C74,0x1CB7,0x0C6B,
/* 3*091#025*[++..-.]*/0x1AC5,0x12CA,0x1C45,0x1C84,0x1C38,0x1C03,0x1C90,0x1C19,0x0C21,
/* 4*122#057*[+..+.-]*/0x1C45,0x1CB6,0x1C3B,0x1C03,0x1C90,0x1C59,0x1C84,0x1C78,0x0C67,
/* 5*082#136*[.++.-.]*/0x16C5,0x1C2A,0x13CB,0x1C6B,0x1C95,0x1C09,0x1C10,0x1CA1,0x0C32,
/* 6*074#146*[++...-]*/0x16C5,0x1C6A,0x1C95,0x1C74,0x17C8,0x1C08,0x1C10,0x1CA1,0x0C49,
/* 7*094#027*[..++-.]*/0x1CA6,0x1C80,0x1C78,0x1C67,0x1C1A,0x1C21,0x1CB2,0x1C3B,0x0C03,

// Case 7.4.1 (1056)
/* 0*037#257*/0x1A65,0x1419,0x11B2,0x17B4,0x04B1,
/* 1*088#134*/0x174B,0x1149,0x12B1,0x11B4,0x0830,
/* 2*026#346*/0x12A5,0x1485,0x1B76,0x1832,0x0582,
/* 3*091#025*/0x1A25,0x1238,0x1458,0x1528,0x0019,
/* 4*122#057*/0x1965,0x1390,0x1B63,0x1693,0x0784,
/* 5*082#136*/0x1695,0x112A,0x1093,0x1B36,0x0396,
/* 6*074#146*/0x1495,0x176A,0x110A,0x1870,0x07A0,
/* 7*094#027*/0x17A6,0x1780,0x11A0,0x1A70,0x03B2,

// Case 7.4.2 (1096)
/* 0*(06)*037#257*/0x1465,0x1459,0x121A,0x1AB2,0x1BA6,0x17B6,0x1476,0x15A1,0x0951,
/* 1*(06)*088#134*/0x1084,0x1748,0x18B7,0x1B83,0x12B3,0x1109,0x1130,0x1123,0x0904,
/* 2*(17)*026#346*/0x16A5,0x132B,0x1B83,0x1874,0x18B7,0x1576,0x1547,0x16B2,0x0A62,
/* 3*(17)*091#025*/0x1945,0x115A,0x1038,0x1023,0x1201,0x1A21,0x1519,0x1908,0x0498,
/* 4*(24)*122#057*/0x1465,0x1380,0x1890,0x1984,0x1594,0x1647,0x1B67,0x1783,0x0B73,
/* 5*(24)*082#136*/0x16A5,0x1A95,0x19A1,0x1091,0x1312,0x1301,0x1B32,0x12A6,0x0B26,
/* 6*(35)*074#146*/0x16A5,0x1109,0x19A1,0x1A95,0x1754,0x1765,0x1874,0x1490,0x0840,
/* 7*(35)*094#027*/0x1BA6,0x1380,0x1378,0x173B,0x167B,0x1AB2,0x11A2,0x1230,0x0120,

// Case 8 (1168)
/* 0*015#0123*/0x1B9A,0x09B8,
/* 1*102#0347*/0x1426,0x0024,
/* 2*051#0145*/0x1375,0x0135,

// Case 9 (1174)
/* 0*078#0237*/0x17A6,0x1180,0x1781,0x0A71,
/* 1*039#0134*/0x1B42,0x1129,0x1492,0x04B7,
/* 2*027#0125*/0x1A35,0x13A2,0x1853,0x0584,
/* 3*114#0457*/0x1965,0x13B0,0x10B6,0x0069,

// Case 10.1.1 (1190)
/* 0*105#0356*[-.-...]*/0x146A,0x194A,0x1028,0x082B,
/* 1*060#0167*[.-.-..]*/0x17A5,0x1189,0x1381,0x07BA,
/* 2*085#0246*[....--]*/0x1625,0x1340,0x1743,0x0521,
//(1202)
/* 0*105#0356*[+.+...]*/0x1846,0x190A,0x186B,0x002A,
/* 1*060#0167*[.+.+..]*/0x1795,0x11BA,0x1789,0x03B1,
/* 2*085#0246*[....++]*/0x1015,0x1236,0x1540,0x0763,

// Case 10.1.2 (1214)
/* 0*(06)(35)*105#0356*[-.-...]*/0x126A,0x1029,0x12A9,0x1B62,0x16B8,0x1468,0x1489,0x0809,
/* 1*(06)(17)*060#0167*[.-.-..]*/0x19A5,0x1789,0x1957,0x11A9,0x1A13,0x1BA3,0x1B37,0x0387,
/* 2*(06)(24)*085#0246*[....--]*/0x1405,0x1756,0x1015,0x1021,0x1320,0x1237,0x1627,0x0745,
//(1238)
/* 0*(17)(24)*105#0356*[+.+...]*/0x146A,0x1A94,0x1904,0x1408,0x1B80,0x1B02,0x1AB2,0x0A6B,
/* 1*(35)(24)*060#0167*[.+.+..]*/0x1BA5,0x1189,0x1138,0x13B8,0x18B7,0x157B,0x115A,0x0195,
/* 2*(17)(35)*085#0246*[....++]*/0x1465,0x1156,0x1621,0x1231,0x1130,0x1403,0x1437,0x0647,

// Case 10.2 (1262)
/* 0*105#0356*[+.-...]*/0x19CA,0x1AC6,0x190C,0x102C,0x12BC,0x18CB,0x184C,0x046C,
/* 1*060#0167*[.+.-..]*/0x1C95,0x1CBA,0x1C7B,0x1C57,0x19C8,0x1C38,0x1C13,0x01CA,
/* 2*085#0246*[....-+]*/0x1C15,0x12C6,0x1C21,0x14C5,0x1C76,0x17C3,0x1C03,0x0C40,
//(1286)
/* 0*105#0356*[-.+...]*/0x1C46,0x1C2A,0x1C94,0x1CA9,0x12C0,0x1C80,0x1CB8,0x0BC6,
/* 1*060#0167*[.-.+..]*/0x1CA5,0x17C5,0x1CBA,0x1C3B,0x11C9,0x13C1,0x1C89,0x08C7,
/* 2*085#0246*[....+-]*/0x16C5,0x12C6,0x123C,0x174C,0x137C,0x10C4,0x101C,0x015C,

// Case 11 (1310)
/* 0*077#0236*/0x16B5,0x1B80,0x15B0,0x0150,
/* 1*046#0137*/0x1786,0x1189,0x1126,0x0681,
/* 2*023#0124*/0x129A,0x1974,0x1792,0x0372,
/* 3*116#0467*/0x14A5,0x13BA,0x13A4,0x0340,
/* 4*099#0345*/0x1975,0x1902,0x1927,0x072B,
/* 5*057#0156*/0x116A,0x1846,0x1861,0x0813,

// Case 14 (1334)
/* 0*113#0456*/0x176A,0x1A90,0x17A0,0x0370,
/* 1*071#0234*/0x1B1A,0x1140,0x11B4,0x074B,
/* 2*043#0135*/0x1125,0x1285,0x12B8,0x0458,
/* 3*029#0126*/0x1695,0x1236,0x1396,0x0389,
/* 4*054#0147*/0x1496,0x1139,0x1369,0x063B,
/* 5*108#0367*/0x12A5,0x1785,0x1825,0x0802,

// Case 12.1.1 (1358)
/* 0*089#0256*[-...-.]*/0x1246,0x1192,0x1429,0x0038,
/* 1*075#0235*[--....]*/0x1945,0x181A,0x1180,0x0B8A,
/* 2*045#0136*[.--...]*/0x16A5,0x1129,0x1B92,0x089B,
/* 3*053#0146*[.-...-]*/0x16A5,0x1749,0x1179,0x0371,
/* 4*030#0127*[..--..]*/0x123A,0x17B6,0x18A3,0x09A8,
/* 5*101#0346*[..-..-]*/0x16A5,0x1274,0x172B,0x0024,
/* 6*092#0267*[...--.]*/0x1715,0x17B2,0x1172,0x0380,
/* 7*120#0567*[-..-..]*/0x19BA,0x1794,0x1B97,0x0803,
/* 8*086#0247*[..-.-.]*/0x1406,0x121A,0x103B,0x060B,
/* 9*083#0245*[.-..-.]*/0x1905,0x121A,0x1350,0x0753,
/*10*058#0157*[...-.-]*/0x1345,0x17B6,0x1384,0x0135,
/*11*106#0357*[-....-]*/0x1945,0x1786,0x1068,0x0260,
//(1406)
/* 0*089#0256*[+...+.]*/0x1246,0x1384,0x1342,0x0190,
/* 1*075#0235*[++....]*/0x1A45,0x14A8,0x18AB,0x0019,
/* 2*045#0136*[.++...]*/0x16B5,0x15B9,0x112A,0x09B8,
/* 3*053#0146*[.+...+]*/0x1495,0x116A,0x1617,0x0713,
/* 4*030#0127*[..++..]*/0x1786,0x168A,0x1A89,0x023B,
/* 5*101#0346*[..+..+]*/0x14A5,0x1B76,0x1A42,0x0240,
/* 6*092#0267*[...++.]*/0x1715,0x1180,0x1817,0x0B23,
/* 7*120#0567*[+..+..]*/0x19BA,0x13B0,0x10B9,0x0784,
/* 8*086#0247*[..+.+.]*/0x11A6,0x1160,0x1064,0x03B2,
/* 9*083#0245*[.+..+.]*/0x1A35,0x123A,0x1537,0x0019,
/*10*058#0157*[...+.+]*/0x1B65,0x1B53,0x1351,0x0784,
/*11*106#0357*[+....+]*/0x1065,0x1905,0x1602,0x0784,

// Case 12.1.2 (1454)
/* 0*(06)*089#0256*[-...-.]*/0x1846,0x1948,0x1980,0x1901,0x1863,0x1362,0x1132,0x0103,
/* 1*(35)*075#0235*[--....]*/0x1AB5,0x1159,0x11A5,0x1190,0x15B4,0x14B8,0x1048,0x0094,
/* 2*(06)*045#0136*[.--...]*/0x1685,0x126A,0x1589,0x11A5,0x12B6,0x16B8,0x12A1,0x0159,
/* 3*(06)*053#0146*[.-...-]*/0x19A5,0x1A36,0x191A,0x1A13,0x1954,0x1637,0x1467,0x0456,
/* 4*(17)*030#0127*[..--..]*/0x126A,0x1738,0x137B,0x1789,0x13B2,0x1796,0x169A,0x02B6,
/* 5*(06)*101#0346*[..-..-]*/0x10A5,0x1B6A,0x1745,0x1756,0x1540,0x176B,0x1A02,0x0BA2,
/* 6*(06)*092#0267*[...--.]*/0x1015,0x1102,0x1203,0x123B,0x1058,0x1857,0x1B87,0x0B38,
/* 7*(06)*120#0567*[-..-..]*/0x13BA,0x1784,0x137B,0x1738,0x13A0,0x10A9,0x1409,0x0480,
/* 8*(24)*086#0247*[..-.-.]*/0x1B6A,0x1BA2,0x1A64,0x1B23,0x1A41,0x1140,0x1310,0x0321,
/* 9*(24)*083#0245*[.-..-.]*/0x19A5,0x1320,0x1021,0x1237,0x1019,0x127A,0x1A75,0x091A,
/*10*(17)*058#0157*[...-.-]*/0x1165,0x1456,0x1467,0x1478,0x161B,0x1B13,0x18B3,0x087B,
/*11*(35)*106#0357*[-....-]*/0x1265,0x1809,0x1894,0x1902,0x1847,0x1925,0x1756,0x0745,
//(1550)
/* 0*(17)*089#0256*[+...+.]*/0x1946,0x1312,0x1013,0x1621,0x1803,0x1961,0x1498,0x0908,
/* 1*(17)*075#0235*[++....]*/0x1945,0x115A,0x1084,0x1904,0x1B80,0x11B0,0x1AB1,0x0195,
/* 2*(24)*045#0136*[.++...]*/0x16A5,0x1195,0x1A15,0x1891,0x1281,0x1B82,0x1B26,0x02A6,
/* 3*(35)*053#0146*[.+...+]*/0x16A5,0x1764,0x1546,0x1374,0x1934,0x1139,0x119A,0x095A,
/* 4*(35)*030#0127*[..++..]*/0x12A6,0x1B26,0x19A2,0x17B6,0x1392,0x1893,0x1837,0x03B7,
/* 5*(17)*101#0346*[..+..+]*/0x16A5,0x1B2A,0x16BA,0x102B,0x1754,0x170B,0x1407,0x0765,
/* 6*(35)*092#0267*[...++.]*/0x1B25,0x178B,0x13B8,0x157B,0x1038,0x1152,0x1120,0x0230,
/* 7*(24)*120#0567*[+..+..]*/0x194A,0x1490,0x1840,0x1380,0x17A4,0x1BA7,0x1B73,0x0783,
/* 8*(35)*086#0247*[..+.+.]*/0x1BA6,0x1130,0x1231,0x1403,0x1A21,0x1B43,0x164B,0x0B2A,
/* 9*(17)*083#0245*[.+..+.]*/0x1A95,0x119A,0x1759,0x121A,0x1079,0x1370,0x1302,0x0012,
/*10*(24)*058#0157*[...+.+]*/0x1465,0x13B8,0x178B,0x1138,0x167B,0x1418,0x1514,0x0476,
/*11*(24)*106#0357*[+....+]*/0x1945,0x1765,0x1754,0x1267,0x1827,0x1028,0x1089,0x0849,

// Case 12.2 (1646)
/* 0*089#0256*[-...+.]*/0x12C6,0x1C19,0x11C0,0x13C2,0x180C,0x18C3,0x194C,0x046C,
/* 1*075#0235*[-+....]*/0x1AC5,0x119C,0x14C9,0x145C,0x1ABC,0x10C8,0x1B8C,0x01C0,
/* 2*045#0136*[.-+...]*/0x1CA5,0x156C,0x12AC,0x16BC,0x1B8C,0x11C9,0x189C,0x02C1,
/* 3*053#0146*[.-...+]*/0x1CA5,0x1AC6,0x14C5,0x16C7,0x137C,0x11C9,0x113C,0x09C4,
/* 4*030#0127*[..-+..]*/0x12CA,0x1CB6,0x13BC,0x178C,0x167C,0x189C,0x19AC,0x03C2,
/* 5*101#0346*[..-..+]*/0x1CA5,0x154C,0x176C,0x1AC6,0x140C,0x1BC2,0x102C,0x07CB,
/* 6*092#0267*[...-+.]*/0x1C15,0x123C,0x101C,0x18C3,0x180C,0x1BC7,0x157C,0x02CB,
/* 7*120#0567*[+..-..]*/0x1CBA,0x1C38,0x17C4,0x1BC7,0x1C84,0x1C90,0x1CA9,0x03C0,
/* 8*086#0247*[..+.-.]*/0x16CA,0x11C2,0x1C03,0x12CB,0x1C3B,0x1C40,0x1C64,0x0AC1,
/* 9*083#0245*[.+..-.]*/0x1AC5,0x1C19,0x11C2,0x13C0,0x10C9,0x1C37,0x1C75,0x02CA,
/*10*058#0157*[...+.-]*/0x1C45,0x17C6,0x178C,0x1BC3,0x16CB,0x113C,0x151C,0x04C8,
/*11*106#0357*[+....-]*/0x1C45,0x1C26,0x184C,0x190C,0x159C,0x102C,0x17C6,0x08C7,
//(1742)
/* 0*089#0256*[+...-.]*/0x146C,0x190C,0x184C,0x13C0,0x138C,0x11C2,0x162C,0x09C1,
/* 1*075#0235*[+-....]*/0x1C45,0x10C9,0x14C8,0x159C,0x1B8C,0x11AC,0x1ABC,0x01C0,
/* 2*045#0136*[.+-...]*/0x16C5,0x16AC,0x15C9,0x11CA,0x189C,0x12BC,0x1B8C,0x02C1,
/* 3*053#0146*[.+...-]*/0x16C5,0x16AC,0x195C,0x1A1C,0x113C,0x1C74,0x137C,0x09C4,
/* 4*030#0127*[..+-..]*/0x16CA,0x12CB,0x17BC,0x17C6,0x19AC,0x138C,0x189C,0x03C2,
/* 5*101#0346*[..+..-]*/0x16C5,0x15CA,0x1BC6,0x1AC2,0x102C,0x174C,0x140C,0x07CB,
/* 6*092#0267*[...+-.]*/0x17C5,0x1C3B,0x18C7,0x103C,0x10C8,0x121C,0x115C,0x02CB,
/* 7*120#0567*[-..+..]*/0x19CA,0x1C80,0x1C94,0x18C7,0x1C47,0x1BC3,0x1CBA,0x03C0,
/* 8*086#0247*[..-.+.]*/0x12CA,0x1CB6,0x1C23,0x1BC3,0x1C64,0x1C01,0x1C40,0x0AC1,
/* 9*083#0245*[.-..+.]*/0x19C5,0x1C1A,0x11C0,0x1C90,0x1C75,0x13C2,0x1C37,0x02CA,
/*10*058#0157*[...-.+]*/0x1C65,0x17C4,0x1BC7,0x1B6C,0x151C,0x18C3,0x113C,0x04C8,
/*11*106#0357*[-....+]*/0x1C65,0x17C4,0x194C,0x19C5,0x126C,0x180C,0x102C,0x08C7,

// Case 13.1 (1838)
/* 0*090#0257*[------]*/0x1945,0x121A,0x17B6,0x0380,
/* 1*090#0257*[++++++]*/0x1A65,0x1190,0x1B23,0x0784,

// Case 13.2 (1846)
/* 0*090#0257*[-----+]*/0x1B65,0x1B59,0x121A,0x1380,0x1794,0x0B97,
/* 1*090#0257*[----+-]*/0x1945,0x17B6,0x181A,0x1801,0x1A38,0x0A23,
/* 2*090#0257*[---+--]*/0x1945,0x121A,0x10B6,0x103B,0x1680,0x0678,
/* 3*090#0257*[--+---]*/0x1945,0x11A6,0x1038,0x1172,0x17B2,0x0167,
/* 4*090#0257*[-+----]*/0x1A45,0x17B6,0x1380,0x1429,0x1219,0x04A2,
/* 5*090#0257*[-+++++]*/0x1A65,0x1794,0x1197,0x13B2,0x1817,0x0801,
/* 6*090#0257*[+-----]*/0x1905,0x121A,0x1345,0x17B6,0x1350,0x0384,
/* 7*090#0257*[+-++++]*/0x1065,0x1590,0x11A6,0x1784,0x1B23,0x0016,
/* 8*090#0257*[++-+++]*/0x1B65,0x135A,0x1784,0x1019,0x1B53,0x0A23,
/* 9*090#0257*[+++-++]*/0x1A65,0x1190,0x1342,0x1384,0x1472,0x07B2,
/*10*090#0257*[++++-+]*/0x1A65,0x1784,0x10B9,0x103B,0x1B29,0x0219,
/*11*090#0257*[+++++-]*/0x1A45,0x18A6,0x1190,0x13B2,0x14A8,0x0678,

// Case 13.3 (1918)
/* 0*090#0257*[---+-+]*/0x1C65,0x1C59,0x121A,0x14C9,0x10C8,0x1C78,0x1C47,0x1C3B,0x16CB,0x0C03,
/* 1*090#0257*[--++--]*/0x1945,0x1CA6,0x13C0,0x11C2,0x1CB2,0x1C3B,0x1C80,0x1C78,0x17C6,0x0C1A,
/* 2*090#0257*[--+--+]*/0x1C65,0x1CA6,0x19C5,0x11AC,0x1C21,0x1CB2,0x17C4,0x1BC7,0x1C94,0x0803,
/* 3*090#0257*[---++-]*/0x1945,0x12CA,0x1CB6,0x1C3B,0x1C23,0x1C1A,0x1C01,0x1C78,0x10C8,0x0C67,
/* 4*090#0257*[--++++]*/0x1C65,0x1A6C,0x17C4,0x159C,0x194C,0x178C,0x180C,0x11AC,0x11C0,0x03B2,
/* 5*090#0257*[--+-+-]*/0x1945,0x1CA6,0x17BC,0x18C3,0x1C23,0x1CB2,0x1C67,0x1C01,0x1AC1,0x0C80,
/* 6*090#0257*[-+---+]*/0x1C65,0x1C5A,0x1CB6,0x12CA,0x17C4,0x1C7B,0x1C19,0x14C9,0x1C21,0x0038,
/* 7*090#0257*[-++---]*/0x1AC5,0x1CA6,0x1C45,0x17C6,0x1C94,0x1C19,0x1CB2,0x11C2,0x1C7B,0x0380,
/* 8*090#0257*[-+++-+]*/0x1A65,0x1C3B,0x17C4,0x18C7,0x180C,0x103C,0x1B2C,0x1C19,0x121C,0x094C,
/* 9*090#0257*[-+--+-]*/0x1AC5,0x1C45,0x17B6,0x10C8,0x14C9,0x1C19,0x1C01,0x1C38,0x1C23,0x02CA,
/*10*090#0257*[-++-++]*/0x1A65,0x119C,0x11C0,0x13C2,0x138C,0x180C,0x194C,0x17BC,0x17C4,0x0B2C,
/*11*090#0257*[-++++-]*/0x1AC5,0x1A6C,0x1C19,0x194C,0x145C,0x167C,0x180C,0x18C7,0x101C,0x0B23,
/*12*090#0257*[+-++-+]*/0x1C65,0x16CA,0x12CB,0x159C,0x121C,0x11AC,0x103C,0x10C9,0x13BC,0x0784,
/*13*090#0257*[+--+--]*/0x1C45,0x17C6,0x1C59,0x121A,0x1C84,0x1C78,0x1CB6,0x1C3B,0x1C90,0x03C0,
/*14*090#0257*[+----+]*/0x1C65,0x19C5,0x121A,0x1C38,0x17C4,0x1BC7,0x1C84,0x1C03,0x1C90,0x0CB6,
/*15*090#0257*[+-+++-]*/0x1C45,0x16CA,0x10C9,0x14C8,0x159C,0x101C,0x11AC,0x167C,0x178C,0x023B,
/*16*090#0257*[+--+++]*/0x1C65,0x12CA,0x13C2,0x159C,0x11C0,0x11AC,0x13BC,0x1B6C,0x190C,0x0784,
/*17*090#0257*[+---+-]*/0x19C5,0x12CA,0x1C45,0x17B6,0x1AC1,0x1C01,0x1C90,0x1C84,0x1C23,0x08C3,
/*28*090#0257*[+++--+]*/0x1A65,0x14C8,0x10C9,0x103C,0x138C,0x147C,0x17BC,0x121C,0x12CB,0x019C,
/*19*090#0257*[++----]*/0x1AC5,0x15C4,0x17B6,0x1C19,0x11C2,0x13C0,0x1C90,0x1CA2,0x1C84,0x0C38,
/*20*090#0257*[++-+-+]*/0x1AC5,0x1B6C,0x1C19,0x1A2C,0x121C,0x190C,0x103C,0x1BC3,0x165C,0x0784,
/*21*090#0257*[+++-+-]*/0x1AC5,0x16CA,0x12CB,0x145C,0x167C,0x17BC,0x123C,0x138C,0x14C8,0x0019,
/*22*090#0257*[++--++]*/0x1C65,0x15AC,0x17C4,0x17BC,0x1B6C,0x1A2C,0x138C,0x13C2,0x184C,0x0190,
/*23*090#0257*[++-++-]*/0x1AC5,0x1B6C,0x145C,0x1C78,0x1BC3,0x167C,0x184C,0x1A2C,0x123C,0x0019,

// Case 13.4 (2158)
/* 0*090#0257*[++---+]*/0x1C65,0x1C5A,0x1C90,0x1C03,0x1C38,0x1C84,0x1C47,0x1C7B,0x1CB6,0x1CA2,0x1C21,0x0C19,
/* 1*090#0257*[-++-+-]*/0x1AC5,0x1A6C,0x138C,0x123C,0x1B2C,0x145C,0x17BC,0x167C,0x194C,0x119C,0x101C,0x080C,
/* 2*090#0257*[--++-+]*/0x1C65,0x1CA6,0x1CB2,0x1C59,0x1C21,0x1C1A,0x1C94,0x1C47,0x1C78,0x1C80,0x1C03,0x0C3B,
/* 3*090#0257*[+--++-]*/0x1C45,0x1B6C,0x159C,0x11AC,0x101C,0x190C,0x184C,0x178C,0x167C,0x13BC,0x123C,0x0A2C,

// Case 13.5.2 (2206)
/* 0*(06)*090#0257*[-++--+]*/0x1A65,0x1784,0x1B87,0x18B3,0x1804,0x1094,0x1190,0x1310,0x1213,0x0B23,
/* 1*(17)*090#0257*[++--+-]*/0x1945,0x115A,0x17B6,0x1380,0x1302,0x1012,0x1908,0x1498,0x1951,0x0A21,
/* 2*(24)*090#0257*[+--+-+]*/0x1A65,0x19A5,0x1A91,0x1A26,0x12B6,0x13B2,0x1132,0x1031,0x1901,0x0784,
/* 3*(35)*090#0257*[--+++-]*/0x1945,0x1BA6,0x1380,0x1837,0x13B7,0x1230,0x1120,0x121A,0x12AB,0x067B,
//(2246)
/* 0*(06)*090#0257*[-++--+]*/0x1465,0x1BA6,0x1519,0x121A,0x12AB,0x15A1,0x1594,0x1476,0x17B6,0x0803,
/* 1*(17)*090#0257*[++--+-]*/0x1A65,0x13B2,0x18B3,0x1574,0x1B87,0x1B62,0x16A2,0x1756,0x1847,0x0190,
/* 2*(24)*090#0257*[+--+-+]*/0x1465,0x1459,0x121A,0x1380,0x1089,0x1849,0x1783,0x1B73,0x17B6,0x0764,
/* 3*(35)*090#0257*[--+++-]*/0x1A65,0x1190,0x1A91,0x19A5,0x1940,0x1480,0x1784,0x1574,0x1675,0x03B2,

// Case 13.5.1 (2286)
/* 0*090#0257*[-++--+]*/0x1A65,0x1380,0x1219,0x1942,0x14B2,0x07B4,
/* 1*090#0257*[++--+-]*/0x1A35,0x1584,0x17B6,0x1190,0x123A,0x0385,
/* 2*090#0257*[+--+-+]*/0x1965,0x121A,0x1B03,0x1B60,0x1690,0x0784,
/* 3*090#0257*[--+++-]*/0x1945,0x17A6,0x13B2,0x1018,0x1178,0x01A7
}
{}

MC33::~MC33() {
	clear_temp_isosurface();
}

int MC33::DefaultColor = DEFAULT_SURFACE_COLOR;

void MC33::set_default_surface_color(unsigned char *color) {
	DefaultColor = *(reinterpret_cast<int*>(color));
}

/******************************************************************
Vertices:           Faces:
    3 __________2        ___________
   /|          /|      /|          /|
  / |         / |     / |   2     / |
7/__________6/  |    /  |     4  /  |
|   |       |   |   |¯¯¯¯¯¯¯¯¯¯¯| 1 |     z
|   0_______|___1   | 3 |_______|___|     |
|  /        |  /    |  /  5     |  /      |____y
| /         | /     | /     0   | /      /
4/__________5/      |/__________|/      x


This function returns a vector with all six test face results (face[6]). Each
result value is 1 if the positive face vertices are joined, -1 if the negative
vertices are joined, and 0 (unchanged) if the test should not be applied. The
return value of this function is the the sum of all six results.*/
int MC33::face_tests(int *face, int i) const {
	if (i&0x80) { //vertex 0
		face[0] = ((i&0xCC) == 0x84? (v[0]*v[5] < v[1]*v[4]? -1: 1): 0);//0x84 = 10000100, vertices 0 and 5
		face[3] = ((i&0x99) == 0x81? (v[0]*v[7] < v[3]*v[4]? -1: 1): 0);//0x81 = 10000001, vertices 0 and 7
		face[4] = ((i&0xF0) == 0xA0? (v[0]*v[2] < v[1]*v[3]? -1: 1): 0);//0xA0 = 10100000, vertices 0 and 2
	} else {
		face[0] = ((i&0xCC) == 0x48? (v[0]*v[5] < v[1]*v[4]? 1: -1): 0);//0x48 = 01001000, vertices 1 and 4
		face[3] = ((i&0x99) == 0x18? (v[0]*v[7] < v[3]*v[4]? 1: -1): 0);//0x18 = 00011000, vertices 3 and 4
		face[4] = ((i&0xF0) == 0x50? (v[0]*v[2] < v[1]*v[3]? 1: -1): 0);//0x50 = 01010000, vertices 1 and 3
	}
	if (i&0x02) { //vertex 6
		face[1] = ((i&0x66) == 0x42? (v[1]*v[6] < v[2]*v[5]? -1: 1): 0);//0x42 = 01000010, vertices 1 and 6
		face[2] = ((i&0x33) == 0x12? (v[3]*v[6] < v[2]*v[7]? -1: 1): 0);//0x12 = 00010010, vertices 3 and 6
		face[5] = ((i&0x0F) == 0x0A? (v[4]*v[6] < v[5]*v[7]? -1: 1): 0);//0x0A = 00001010, vertices 4 and 6
	} else {
		face[1] = ((i&0x66) == 0x24? (v[1]*v[6] < v[2]*v[5]? 1: -1): 0);//0x24 = 00100100, vertices 2 and 5
		face[2] = ((i&0x33) == 0x21? (v[3]*v[6] < v[2]*v[7]? 1: -1): 0);//0x21 = 00100001, vertices 2 and 7
		face[5] = ((i&0x0F) == 0x05? (v[4]*v[6] < v[5]*v[7]? 1: -1): 0);//0x05 = 00000101, vertices 5 and 7
	}
	return face[0] + face[1] + face[2] + face[3] + face[4] + face[5];
}

/* Faster function for the face test, the test is applied to only one face
(int face). This function is only used for the cases 3 and 6 of MC33*/
int MC33::face_test1(int face) const {
	switch (face) {
	case 0:
		return (v[0]*v[5] < v[1]*v[4]? 0x48: 0x84);
	case 1:
		return (v[1]*v[6] < v[2]*v[5]? 0x24: 0x42);
	case 2:
		return (v[3]*v[6] < v[2]*v[7]? 0x21: 0x12);
	case 3:
		return (v[0]*v[7] < v[3]*v[4]? 0x18: 0x81);
	case 4:
		return (v[0]*v[2] < v[1]*v[3]? 0x50: 0xA0);
	default:
		return (v[4]*v[6] < v[5]*v[7]? 0x05: 0x0A);
	}
}

// Silence dereferencing type-punned pointer warning in GCC
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif
// an ugly signbit:
#if MC33_DOUBLE_PRECISION
inline unsigned int signbf(double x) {
	return ((*(reinterpret_cast<unsigned long long int*>(&x))>>32)&0x80000000);
}
#else
inline unsigned int signbf(float x) {
	return (*(reinterpret_cast<unsigned int*>(&x))&0x80000000);
}
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

/******************************************************************
Interior test function. If the test is positive, the function returns a value
different from 0. The integer i must be 0 to test if the vertices 0 and 6 are
joined. 1 for vertices 1 and 7, 2 for vertices 2 and 4, and 3 for 3 and 5.
For case 13, the integer flag13 must be 1, and the function returns 2 if one
of the vertices 0, 1, 2 or 3 is joined to the center point of the cube (case
13.5.2), returns 1 if one of the vertices 4, 5, 6 or 7 is joined to the
center point of the cube (case 13.5.2 too), and it returns 0 if the vertices
are not joined (case 13.5.1)*/
int MC33::interior_test(int i, int flag13) const {
	//Signs of cube vertices were changed to use signbit function in calc_isosurface
	//A0 = -v[0], B0 = -v[1], C0 = -v[2], D0 = -v[3]
	//A1 = -v[4], B1 = -v[5], C1 = -v[6], D1 = -v[7]
	//But the function still works
	MC33_real At = v[4] - v[0], Bt = v[5] - v[1],
				Ct = v[6] - v[2], Dt = v[7] - v[3];
	MC33_real t = At*Ct - Bt*Dt;//the "a" value.
	if (signbf(t)) {
		if (i&0x01) return 0;
	} else {
		if (!(i&0x01) || t == 0) return 0;
	}
	t = 0.5f*(v[3]*Bt + v[1]*Dt - v[2]*At - v[0]*Ct)/t;//t = -b/2a

	if (t > 0 && t < 1) {
		At = v[0] + At*t;
		Bt = v[1] + Bt*t;
		Ct = v[2] + Ct*t;
		Dt = v[3] + Dt*t;
		Ct *= At;
		Dt *= Bt;
		if (i&0x01) {
			if (Ct < Dt && signbf(Dt) == 0)
				return (signbf(Bt) == signbf(v[i])) + flag13;
		} else {
			if (Ct > Dt && signbf(Ct) == 0)
				return (signbf(At) == signbf(v[i])) + flag13;
		}
	}
	return 0;
}

/******************************************************************
Auxiliary function that calculates the normal if a vertex
of the cube lies on the isosurface.
*/
unsigned int MC33::surfint(unsigned int x, unsigned int y, unsigned int z, MC33_real *r) {
	r[0] = x; r[1] = y; r[2] = z;
	if (x == 0)
		r[3] = F[z][y][0] - F[z][y][di];
	else if (x == nx) {
		x *= di;
		r[3] = F[z][y][x - di] - F[z][y][x];
	} else {
		x *= di;
		r[3] = 0.5f*(F[z][y][x - di] - F[z][y][x + di]);
	}
	if (y == 0)
		r[4] = F[z][0][x] - F[z][1][x];
	else if (y == ny)
		r[4] = F[z][y - 1][x] - F[z][y][x];
	else
		r[4] = 0.5f*(F[z][y - 1][x] - F[z][y + 1][x]);
	if (z == 0)
		r[5] = F[0][y][x] - F[1][y][x];
	else if (z == nz)
		r[5] = F[z - 1][y][x] - F[z][y][x];
	else
		r[5] = 0.5f*(F[z - 1][y][x] - F[z + 1][y][x]);
	return store_point(r);
}

/******************************************************************
This function find the MC33 case (using the index i, and the face and interior
tests). The correct triangle pattern is obtained from the arrays contained in
the MC33_LookUpTable.h file. The necessary vertices (intersection points) are
also calculated here (using trilinear interpolation).
       _____2_____
     /|          /|
   11 |<-3     10 |
   /____6_____ /  1     z
  |   |       |   |     |
  |   |_____0_|___|     |____y
  7  /        5  /     /
  | 8         | 9     x
  |/____4_____|/

The temporary matrices: Dx, Dy, Ux, Uy, and Lz are filled
and used here.*/
#define FF 0xFFFFFFFF
void MC33::find_case(unsigned int x, unsigned int y, unsigned int z, unsigned int i) {
	unsigned int p[13] = {FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF};
	unsigned int ti[3];//for vertex indices of a triangle
	union { // memory saving
		int f[6];//for the face tests
		MC33_real r[6];//for intercept and normal coordinates
	};
	const unsigned short int *pcase = table;
	unsigned int k, c, m, n;
	MC33_real t;
	if (i&0x80) {
		c = pcase[i^0xFF];
		m = (c&0x800) == 0;
		n = !m;
	} else {
		c = pcase[i];
		n = (c&0x800) == 0;
		m = !n;
	}
	k = c&0x7FF;
	switch (c>>12) { //find the MC33 case
		case 0: // case 1, 2, 5, 8, 9, 11 and 14
			pcase += k;
			break;
		case 1: // case 3
			pcase += ((m? i: i^0xFF)&face_test1(k>>2)? 183 + (k<<1): 159 + k);
			break;
		case 2: // case 4
			pcase += (interior_test(k,0)? 239 + 6*k: 231 + (k<<1));
			break;
		case 3: // case 6
			if ((m? i: i^0xFF)&face_test1(k%6))
				pcase += 575 + 5*k; //6.2
			else
				pcase += (interior_test(k/6,0)? 407 + 7*k: 335 + 3*k); //6.1
			break;
		case 4: // case 7
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -3:
				pcase += 695 + 3*k; //7.1
				break;
			case -1: //7.2
				pcase += (f[4] + f[5] < 0? (f[0] + f[2] < 0? 759: 799): 719) + 5*k;
				break;
			case 1: //7.3
				pcase += (f[4] + f[5] < 0? 983: (f[0] + f[2] < 0? 839: 911)) + 9*k;
				break;
			default: //7.4
				pcase += (interior_test(k>>1,0)? 1095 + 9*k: 1055 + 5*k);
			}
			break;
		case 5: // case 10
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -2:
				if (k&2? interior_test(0,0): interior_test(0,0)||interior_test(k? 1: 3,0))
					pcase += 1213 + (k<<3); //10.1.2
				else
					pcase += 1189 + (k<<2); //10.1.1
				break;
			case 0: //10.2
				pcase += (f[2 + k] < 0? 1261: 1285) + (k<<3);
				break;
			default:
				if (k&2? interior_test(1,0): interior_test(2,0)||interior_test(k? 3: 1,0))
					pcase += 1237 + (k<<3); //10.1.2
				else
					pcase += 1201 + (k<<2); //10.1.1
			}
			break;
		case 6: // case 12
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -2: //12.1
				pcase += (interior_test((0xDA010C>>(k<<1))&3,0)? 1453 + (k<<3): 1357 + (k<<2));
				break;
			case 0: //12.2
				pcase += (f[k>>1] < 0? 1645: 1741) + (k<<3);
				break;
			default: //12.1
				pcase += (interior_test((0xA7B7E5>>(k<<1))&3,0)? 1549 + (k<<3): 1405 + (k<<2));
			}
			break;
		default: // case 13
			switch (abs(face_tests(f, 165))) {
			case 0:
				k = ((f[1] < 0)<<1)|(f[5] < 0);
				if (f[0]*f[1] == f[5]) //13.4
					pcase += 2157 + 12*k;
				else {
					c = interior_test(k, 1); // 13.5.1 if c == 0 else 13.5.2
					pcase += 2285 + (c? 10*k - 40*c: 6*k);
				}
				break;
			case 2: //13.3
				pcase += 1917 + 10*((f[0] < 0? f[2] > 0: 12 + (f[2] < 0)) + (f[1] < 0? f[3] < 0: 6 + (f[3] > 0)));
				if (f[4] > 0)
					pcase += 30;
				break;
			case 4: //13.2
				k = 21 + 11*f[0] + 4*f[1] + 3*f[2] + 2*f[3] + f[4];
				if (k >> 4)
					k -= (k&32? 20: 10);
				pcase += 1845 + 3*k;
				break;
			default: //13.1
				pcase += 1839 + 2*f[0];
			}
	}
	while (i) {
		i = *(++pcase);
		for (k = 3; k;) {
			c = i&0x0F;
			i >>= 4;
			if (p[c] == FF) {
				switch (c) { // the vertices r[3] and normals (r + 3)[3] are calculated here
				case 0:
					if (z || x)
						ti[--k] = p[0] = Dy[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[0] = p[3];
							else if (p[8] != FF)
								p[0] = p[8];
							else if (y && signbf(v[3]))
								p[0] = Lz[y][0];
							else if (y && signbf(v[4]))
								p[0] = Dx[y][0];
							else if (y? signbf(S->iso - F[0][y - 1][0]): 0)
								p[0] = Dy[y - 1][0];
							else
								p[0] = surfint(0,y,0,r);
						} else if (v[1] == 0) {
							if (p[9] != FF)
								p[0] = p[9];
							else
								p[0] = (p[1] != FF? p[1]: surfint(0,y + 1,0,r));
						} else {
							t = v[0]/(v[0] - v[1]);
							r[0] = r[2] = 0;
							r[1] = y + t;
							r[3] = (v[4] - v[0])*(1 - t) + (v[5] - v[1])*t;
							r[4] = v[1] - v[0];
							r[5] = (v[3] - v[0])*(1 - t) + (v[2] - v[1])*t;
							p[0] = store_point(r);
						}
						Dy[y][0] = ti[--k] = p[0];
					}
					break;
				case 1:
					if (x)
						ti[--k] = p[1] = Lz[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[1] = p[0];
							else if (p[9] != FF)
								p[1] = p[9];
							else if (z && signbf(v[0]))
								p[1] = Dy[y][0];
							//else if (z && signbf(v[5]))
							//	p[1] = Dx[y + 1][0];
							else if (z && y + 1 < ny? signbf(S->iso - F[z][y + 2][0]): 0)
								p[1] = Dy[y + 1][0];
							else if (z? signbf(S->iso - F[z - 1][y + 1][0]): 0) {
								ti[--k] = p[1] = Lz[y + 1][0]; // value of previous slice
								break;
							} else
								p[1] = surfint(0,y + 1,z,r);
						} else if (v[2] == 0) {
							if (p[10] != FF)
								p[1] = p[10];
							else
								p[1] = (p[2] != FF? p[2]: surfint(0,y + 1,z + 1,r));
						} else {
							t = v[1]/(v[1] - v[2]);
							r[0] = 0; r[1] = y + 1;
							r[2] = z + t;
							r[3] = (v[5] - v[1])*(1 - t) + (v[6] - v[2])*t;
							r[4] = (y + 1 < ny? 0.5f*((F[z][y][0] - F[z][y + 2][0])*(1 - t)
										+ (F[z + 1][y][0] - F[z + 1][y + 2][0])*t):
										(v[1] - v[0])*(1 - t) + (v[2] - v[3])*t);
							r[5] = v[2] - v[1];
							p[1] = store_point(r);
						}
						Lz[y + 1][0] = ti[--k] = p[1];
					}
					break;
				case 2:
					if (x)
						ti[--k] = p[2] = Uy[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[2] = p[3];
							else if (p[11] != FF)
								p[2] = p[11];
							else if (y && signbf(v[0]))
								p[2] = Lz[y][0];
							else if (y && signbf(v[7]))
								p[2] = Ux[y][0];
							else if (y? signbf(S->iso - F[z + 1][y - 1][0]): 0)
								p[2] = Uy[y - 1][0];
							else
								p[2] = surfint(0,y,z + 1,r);
						} else if (v[2] == 0) {
							if (p[10] != FF)
								p[2] = p[10];
							else
								p[2] = (p[1] != FF? p[1]: surfint(0,y + 1,z + 1,r));
						} else {
							t = v[3]/(v[3] - v[2]);
							r[0] = 0; r[2] = z + 1;
							r[1] = y + t;
							r[3] = (v[7] - v[3])*(1 - t) + (v[6] - v[2])*t;
							r[4] = v[2] - v[3];
							r[5] = (z + 1 < nz? 0.5f*((F[z][y][0] - F[z + 2][y][0])*(1 - t)
										+ (F[z][y + 1][0] - F[z + 2][y + 1][0])*t):
										(v[3] - v[0])*(1 - t) + (v[2] - v[1])*t);
							p[2] = store_point(r);
						}
						Uy[y][0] = ti[--k] = p[2];
					}
					break;
				case 3:
					if (y || x)
						ti[--k] = p[3] = Lz[y][x];
					else {
						if (v[0] == 0) {
							if (p[0] != FF)
								p[3] = p[0];
							else if (p[8] != FF)
								p[3] = p[8];
							else if (z && signbf(v[1]))
								p[3] = Dy[0][0];
							else if (z && signbf(v[4]))
								p[3] = Dx[0][0];
							else if (z? signbf(S->iso - F[z - 1][0][0]): 0) {
								ti[--k] = p[3] = Lz[0][0]; // value of previous slice
								break;
							} else
								p[3] = surfint(0,0,z,r);
						} else if (v[3] == 0) {
							if (p[2] != FF)
								p[3] = p[2];
							else
								p[3] = (p[11] != FF? p[11]: surfint(0,0,z + 1,r));
						} else {
							t = v[0]/(v[0] - v[3]);
							r[0] = r[1] = 0;
							r[2] = z + t;
							r[3] = (v[4] - v[0])*(1 - t) + (v[7] - v[3])*t;
							r[4] = (v[1] - v[0])*(1 - t) + (v[2] - v[3])*t;
							r[5] = v[3] - v[0];
							p[3] = store_point(r);
						}
						Lz[0][0] = ti[--k] = p[3];
					}
					break;
				case 4:
					if (z)
						ti[--k] = p[4] = Dy[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[4] = p[8];
							//else if (p[7] != FF)
							//	p[4] = p[7];
							else if (y && signbf(v[7]))
								p[4] = Lz[y][x + 1];
							else if (y && signbf(v[0]))
								p[4] = Dx[y][x];
							else if (y? signbf(S->iso - F[0][y - 1][di*(x + 1)]): 0)
								p[4] = Dy[y - 1][x + 1];
							else if (y && x + 1 < nx? signbf(S->iso - F[0][y][di*(x + 2)]): 0)
								p[4] = Dx[y][x + 1];
							else
								p[4] = surfint(x + 1,y,0,r);
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[4] = p[5];
							else
								p[4] = (p[9] != FF? p[9]: surfint(x + 1,y + 1,0,r));
						} else {
							t = v[4]/(v[4] - v[5]);
							r[0] = x + 1; r[2] = 0;
							r[1] = y + t;
							r[3] = (x + 1 < nx? 0.5f*((F[0][y][di*x] - F[0][y][di*(x + 2)])*(1 - t)
										+ (F[0][y + 1][di*x] - F[0][y + 1][di*(x + 2)])*t):
										(v[4] - v[0])*(1 - t) + (v[5] - v[1])*t);
							r[4] = v[5] - v[4];
							r[5] = (v[7] - v[4])*(1 - t) + (v[6] - v[5])*t;
							p[4] = store_point(r);
						}
						Dy[y][x + 1] = ti[--k] = p[4];
					}
					break;
				case 5:
					if (v[5] == 0) {
						if (z) {
							if (signbf(v[4]))
								p[5] = p[4] = Dy[y][x + 1];
							else if (signbf(v[1]))
								p[5] = p[9] = Dx[y + 1][x];
							else if (x + 1 < nx? signbf(S->iso - F[z][y + 1][di*(x + 2)]): 0)
								p[5] = Dx[y + 1][x + 1];
							else if (y + 1 < ny? signbf(S->iso - F[z][y + 2][di*(x + 1)]): 0)
								p[5] = Dy[y + 1][x + 1];
							else if (signbf(S->iso - F[z - 1][y + 1][di*(x + 1)])) {
								ti[--k] = p[5] = Lz[y + 1][x + 1]; // value of previous slice
								break;
							} else
								p[5] = surfint(x + 1,y + 1,z,r);
						} else
							p[5] = surfint(x + 1,y + 1,0,r);
					} else if (v[6] == 0)
						p[5] = surfint(x + 1,y + 1,z + 1,r);
					else {
						t = v[5]/(v[5] - v[6]);
						r[0] = x + 1; r[1] = y + 1;
						r[2] = z + t;
						c = di*x;
						r[3] = (x + 1 < nx? 0.5f*((F[z][y + 1][c] - F[z][y + 1][c + 2*di])*(1 - t)
									+ (F[z + 1][y + 1][c] - F[z + 1][y + 1][c + 2*di])*t):
									(v[5] - v[1])*(1 - t) + (v[6] - v[2])*t);
						r[4] = (y + 1 < ny? 0.5f*((F[z][y][c + di] - F[z][y + 2][c + di])*(1 - t)
									+ (F[z + 1][y][c + di] - F[z + 1][y + 2][c + di])*t):
									(v[5] - v[4])*(1 - t) + (v[6] - v[7])*t);
						r[5] = v[6] - v[5];
						p[5] = store_point(r);
					}
					Lz[y + 1][x + 1] = ti[--k] = p[5];
					break;
				case 6:
					if (v[7] == 0) {
						if (y) {
							if (signbf(v[3]))
								p[6] = p[11] = Ux[y][x];
							else if (signbf(v[4]))
								p[6] = p[7] = Lz[y][x + 1];
							else if (signbf(S->iso - F[z + 1][y - 1][di*(x + 1)]))
								p[6] = Uy[y - 1][x + 1];
							else if (x + 1 < nx? signbf(S->iso - F[z + 1][y][di*(x + 2)]): 0)
								p[6] = Ux[y][x + 1];
							else
								p[6] = surfint(x + 1,y,z + 1,r);
						} else if (p[11] != FF)
								p[6] = p[11];
							//else if (p[7] != FF)
							//	p[6] = p[7];
							else
								p[6] = surfint(x + 1,0,z + 1,r);
					} else if (v[6] == 0) {
						if (p[5] != FF)
							p[6] = p[5];
						else
							p[6] = (p[10] != FF? p[10]: surfint(x + 1,y + 1,z + 1,r));
					} else {
						t = v[7]/(v[7] - v[6]);
						r[0] = x + 1;
						r[1] = y + t; r[2] = z + 1;
						c = di*x;
						r[3] = (x + 1 < nx? 0.5f*((F[z + 1][y][c] - F[z + 1][y][c + 2*di])*(1 - t)
									+ (F[z + 1][y + 1][c] - F[z + 1][y + 1][c + 2*di])*t):
									(v[7] - v[3])*(1 - t) + (v[6] - v[2])*t);
						r[4] = v[6] - v[7];
						r[5] = (z + 1 < nz? 0.5f*((F[z][y][c + di] - F[z + 2][y][c + di])*(1 - t)
										+ (F[z][y + 1][c + di] - F[z + 2][y + 1][c + di])*t):
									(v[7] - v[4])*(1 - t) + (v[6] - v[5])*t);
						p[6] = store_point(r);
					}
					Uy[y][x + 1] = ti[--k] = p[6];
					break;
				case 7:
					if (y)
						ti[--k] = p[7] = Lz[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[7] = p[8];
							else if (p[4] != FF)
								p[7] = p[4];
							else if (z && signbf(v[0]))
								p[7] = Dx[0][x];
							//else if (z && signbf(v[5]))
							//	p[7] = Dy[0][x + 1];
							else if (z && x + 1 < nx? signbf(S->iso - F[z][0][di*(x + 2)]): 0)
								p[7] = Dx[0][x + 1];
							else if (z? signbf(S->iso - F[z - 1][0][di*(x + 1)]): 0) {
								ti[--k] = p[7] = Lz[0][x + 1]; // value of previous slice
								break;
							} else
								p[7] = surfint(x + 1,0,z,r);
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[7] = p[6];
							else
								p[7] = (p[11] != FF? p[11]: surfint(x + 1,0,z + 1,r));
						} else {
							t = v[4]/(v[4] - v[7]);
							r[0] = x + 1; r[1] = 0;
							r[2] = z + t;
							r[3] = (x + 1 < nx? 0.5f*((F[z][0][di*x] - F[z][0][di*(x + 2)])*(1 - t)
										+ (F[z + 1][0][di*x] - F[z + 1][0][di*(x + 2)])*t):
										(v[4] - v[0])*(1 - t) + (v[7] - v[3])*t);
							r[4] = (v[5] - v[4])*(1 - t) + (v[6] - v[7])*t;
							r[5] = v[7] - v[4];
							p[7] = store_point(r);
						}
						Lz[0][x + 1] = ti[--k] = p[7];
					}
					break;
				case 8:
					if (z || y)
						ti[--k] = p[8] = Dx[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[8] = p[3];
							else if (p[0] != FF)
								p[8] = p[0];
							else if (x && signbf(v[3]))
								p[8] = Lz[0][x];
							else if (x && signbf(v[1]))
								p[8] = Dy[0][x];
							else if (x? signbf(S->iso - F[0][0][di*(x - 1)]): 0)
								p[8] = Dx[0][x - 1];
							else
								p[8] = surfint(x,0,0,r);
						} else if (v[4] == 0) {
							if (p[4] != FF)
								p[8] = p[4];
							else
								p[8] = (p[7] != FF? p[7]: surfint(x + 1,0,0,r));
						} else {
							t = v[0]/(v[0] - v[4]);
							r[1] = r[2] = 0;
							r[0] = x + t;
							r[3] = v[4] - v[0];
							r[4] = (v[1] - v[0])*(1 - t) + (v[5] - v[4])*t;
							r[5] = (v[3] - v[0])*(1 - t) + (v[7] - v[4])*t;
							p[8] = store_point(r);
						}
						Dx[0][x] = ti[--k] = p[8];
					}
					break;
				case 9:
					if (z)
						ti[--k] = p[9] = Dx[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[9] = p[0];
							//else if (p[1] != FF)
							//	p[9] = p[1];
							else if (x && signbf(v[0]))
								p[9] = Dy[y][x];
							else if (x && signbf(v[2]))
								p[9] = Lz[y + 1][x];
							else if (x? signbf(S->iso - F[0][y + 1][di*(x - 1)]): 0)
								p[9] = Dx[y + 1][x - 1];
							else
								p[9] = surfint(x,y + 1,0,r);
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[9] = p[5];
							else
								p[9] = (p[4] != FF? p[4]: surfint(x + 1,y + 1,0,r));
						} else {
							t = v[1]/(v[1] - v[5]);
							r[1] = y + 1; r[2] = 0;
							r[0] = x + t;
							r[3] = v[5] - v[1];
							r[4] = (y + 1 < ny? 0.5f*((F[0][y][di*x] - F[0][y + 2][di*x])*(1 - t)
										+ (F[0][y][di*(x + 1)] - F[0][y + 2][di*(x + 1)])*t):
										(v[1] - v[0])*(1 - t) + (v[5] - v[4])*t);
							r[5] = (v[2] - v[1])*(1 - t) + (v[6] - v[5])*t;
							p[9] = store_point(r);
						}
						Dx[y + 1][x] = ti[--k] = p[9];
					}
					break;
				case 10:
					if (v[2] == 0) {
						if (x) {
							if (signbf(v[1]))
								p[10] = p[1] = Lz[y + 1][x];
							else if (signbf(v[3]))
								p[10] = p[2] = Uy[y][x];
							else if (signbf(S->iso - F[z + 1][y + 1][di*(x - 1)]))
								p[10] = Ux[y + 1][x - 1];
							else
								p[10] = surfint(x,y + 1,z + 1,r);
						} else if (p[2] != FF)
								p[10] = p[2];
							//else if (p[1] != FF)
							//	p[10] = p[1];
							else
								p[10] = surfint(0,y + 1,z + 1,r);
					} else if (v[6] == 0) {
						if (p[5] != FF)
							p[10] = p[5];
						else
							p[10] = (p[6] != FF? p[6]: surfint(x + 1,y + 1,z + 1,r));
					} else {
						t = v[2]/(v[2] - v[6]);
						r[0] = x + t;
						r[1] = y + 1; r[2] = z + 1;
						r[3] = v[6] - v[2];
						c = di*x;
						r[4] = (y + 1 < ny? 0.5f*((F[z + 1][y][c] - F[z + 1][y + 2][c])*(1 - t)
									+ (F[z + 1][y][c + di] - F[z + 1][y + 2][c + di])*t):
									(v[2] - v[3])*(1 - t) + (v[6] - v[7])*t);
						r[5] = (z + 1 < nz? 0.5f*((F[z][y + 1][c] - F[z + 2][y + 1][c])*(1 - t)
									+ (F[z][y + 1][c + di] - F[z + 2][y + 1][c + di])*t):
									(v[2] - v[1])*(1 - t) + (v[6] - v[5])*t);
						p[10] = store_point(r);
					}
					Ux[y + 1][x] = ti[--k] = p[10];
					break;
				case 11:
					if (y)
						ti[--k] = p[11] = Ux[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[11] = p[3];
							else if (p[2] != FF)
								p[11] = p[2];
							else if (x && signbf(v[0]))
								p[11] = Lz[0][x];
							else if (x && signbf(v[2]))
								p[11] = Uy[0][x];
							else if (x? signbf(S->iso - F[z + 1][0][di*(x - 1)]): 0)
								p[11] = Ux[0][x - 1];
							else
								p[11] = surfint(x,0,z + 1,r);
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[11] = p[6];
							else
								p[11] = (p[7] != FF? p[7]: surfint(x + 1,0,z + 1,r));
						} else {
							t = v[3]/(v[3] - v[7]);
							r[1] = 0; r[2] = z + 1;
							r[0] = x + t;
							r[3] = v[7] - v[3];
							r[4] = (v[2] - v[3])*(1 - t) + (v[6] - v[7])*t;
							r[5] = (z + 1 < nz? 0.5f*((F[z][0][di*x] - F[z + 2][0][di*x])*(1 - t)
										+ (F[z][0][di*(x + 1)] - F[z + 2][0][di*(x + 1)])*t):
										(v[3] - v[0])*(1 - t) + (v[7] - v[4])*t);
							p[11] = store_point(r);
						}
						Ux[0][x] = ti[--k] = p[11];
					}
				break;
				default:
					r[0] = x + 0.5f; r[1] = y + 0.5f; r[2] = z + 0.5f;
					r[3] = v[4] + v[5] + v[6] + v[7] - v[0] - v[1] - v[2] - v[3];
					r[4] = v[1] + v[2] + v[5] + v[6] - v[0] - v[3] - v[4] - v[7];
					r[5] = v[2] + v[3] + v[6] + v[7] - v[0] - v[1] - v[4] - v[5];
					ti[--k] = p[12] = store_point(r);
				}
			} else
				ti[--k] = p[c];//now ti contains the vertex indices of the triangle
		}
		if (ti[0] != ti[1] && ti[0] != ti[2] && ti[1] != ti[2]) { //to avoid zero area triangles
			if (!(S->nT&0x0FFF)) {
				try {
					S->T.resize(S->nT + 0x1000);
				}
				catch (...) {
					memoryfault = 1;
					if (S->nT)
						S->nT = 1;
					return;
				}
			}
			unsigned int *vp = S->T[S->nT++].v;
#ifndef MC33_NORMAL_NEG
			*vp = ti[n]; *(++vp) = ti[m]; *(++vp) = ti[2];
#else
			*vp = ti[m]; *(++vp) = ti[n]; *(++vp) = ti[2];
#endif
		}
	}
}

void MC33::case_count(unsigned int x, unsigned int y, unsigned int z, unsigned int i) {
	unsigned int p[13] = {FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF};
	union {
		int f[6];
		unsigned int ti[3];
	};
	const unsigned short int *pcase = table;
	unsigned int c, k, m = i>>7;
	c = pcase[m? i^0xFF: i];
	m = (m^((c&0xFFF)>>11));
	k = c&0x7FF;
	switch (c>>12) { //find the MC33 case
		case 0: // cases 1, 2, 5, 8, 9, 11 and 14
			pcase += k;
			break;
		case 1: // case 3
			pcase += ((m? i: i^0xFF)&face_test1(k>>2)? 183 + (k<<1): 159 + k);
			break;
		case 2: // case 4
			pcase += (interior_test(k,0)? 239 + 6*k: 231 + (k<<1));
			break;
		case 3: // case 6
			if ((m? i: i^0xFF)&face_test1(k%6))
				pcase += 575 + 5*k; //6.2
			else
				pcase += (interior_test(k/6,0)? 407 + 7*k: 335 + 3*k); //6.1
			break;
		case 4: // case 7
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -3:
				pcase += 695 + 3*k; //7.1
				break;
			case -1: //7.2
				pcase += (f[4] + f[5] < 0? (f[0] + f[2] < 0? 759: 799): 719) + 5*k;
				break;
			case 1: //7.3
				pcase += (f[4] + f[5] < 0? 983: (f[0] + f[2] < 0? 839: 911)) + 9*k;
				break;
			default: //7.4
				pcase += (interior_test(k>>1,0)? 1095 + 9*k: 1055 + 5*k);
			}
			break;
		case 5: // case 10
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -2:
				if (k&2? interior_test(0,0): interior_test(0,0)||interior_test(k? 1: 3,0))
					pcase += 1213 + (k<<3); //10.1.2
				else
					pcase += 1189 + (k<<2); //10.1.1
				break;
			case 0: //10.2
				pcase += (f[2 + k] < 0? 1261: 1285) + (k<<3);
				break;
			default:
				if (k&2? interior_test(1,0): interior_test(2,0)||interior_test(k? 3: 1,0))
					pcase += 1237 + (k<<3); //10.1.2
				else
					pcase += 1201 + (k<<2); //10.1.1
			}
			break;
		case 6: // case 12
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -2: //12.1
				pcase += (interior_test((0xDA010C>>(k<<1))&3,0)? 1453 + (k<<3): 1357 + (k<<2));
				break;
			case 0: //12.2
				pcase += (f[k>>1] < 0? 1645: 1741) + (k<<3);
				break;
			default: //12.1
				pcase += (interior_test((0xA7B7E5>>(k<<1))&3,0)? 1549 + (k<<3): 1405 + (k<<2));
			}
			break;
		default: // case 13
			switch (abs(face_tests(f,(m? 90: 165)))) {
			case 0:
				k = ((f[1] < 0)<<1)|(f[5] < 0);
				if (f[0]*f[1] == f[5]) //13.4
					pcase += 2157 + 12*k;
				else {
					c = interior_test(k, 1); // 13.5.1 if c == 0 else 13.5.2
					pcase += 2285 + (c? 10*k - 40*c: 6*k);
				}
				break;
			case 2: //13.3
				pcase += 1917 + 10*((f[0] < 0? f[2] > 0: 12 + (f[2] < 0)) + (f[1] < 0? f[3] < 0: 6 + (f[3] > 0)));
				if (f[4] > 0)
					pcase += 30;
				break;
			case 4: //13.2
				k = 21 + 11*f[0] + 4*f[1] + 3*f[2] + 2*f[3] + f[4];
				if (k >> 4)
					k -= (k&32? 20: 10);
				pcase += 1845 + 3*k;
				break;
			default: //13.1
				pcase += 1839 + 2*f[0];
			}
	}
	while (i) {
		i = *(++pcase);
		for (k = 3; k;) {
			c = i&0x0F;
			i >>= 4;
			if (p[c] == FF) {
				switch (c) {
				case 0:
					if (z || x)
						p[0] = Dy[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[0] = p[3];
							else if (p[8] != FF)
								p[0] = p[8];
							else if (y && signbf(v[3]))
								p[0] = Lz[y][0];
							else if (y && signbf(v[4]))
								p[0] = Dx[y][0];
							else if (y? signbf(S->iso - F[0][y - 1][0]): 0)
								p[0] = Dy[y - 1][0];
							else
								p[0] = S->nV++;;
						} else if (v[1] == 0) {
							if (p[1] != FF)
								p[0] = p[1];
							else if (p[9] != FF)
								p[0] = p[9];
							else
								p[0] = S->nV++;;
						} else
							p[0] = S->nV++;;
						Dy[y][0] = p[0];
					}
					break;
				case 1:
					if (x)
						p[1] = Lz[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[1] = p[0];
							else if (p[9] != FF)
								p[1] = p[9];
							else if (z && signbf(v[0]))
								p[1] = Dy[y][0];
							else if (z && signbf(v[5]))
								p[1] = Dx[y + 1][0];
							else if (z && y + 1 < ny? signbf(S->iso - F[z][y + 2][0]): 0)
								p[1] = Dy[y + 1][0];
							else if (z? signbf(S->iso - F[z - 1][y + 1][0]): 0)
								p[1] = Lz[y + 1][0];
							else
								p[1] = S->nV++;
						} else if (v[2] == 0) {
							if (p[2] != FF)
								p[1] = p[2];
							else if (p[10] != FF)
								p[1] = p[10];
							else
								p[1] = S->nV++;
						} else
							p[1] = S->nV++;
						Lz[y + 1][0] = p[1];
					}
					break;
				case 2:
					if (x)
						p[2] = Uy[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[2] = p[3];
							else if (p[11] != FF)
								p[2] = p[11];
							else if (y && signbf(v[0]))
								p[2] = Lz[y][0];
							else if (y && signbf(v[7]))
								p[2] = Ux[y][0];
							else if (y? signbf(S->iso - F[z + 1][y - 1][0]): 0)
								p[2] = Uy[y - 1][0];
							else
								p[2] = S->nV++;
						} else if (v[2] == 0) {
							if (p[1] != FF)
								p[2] = p[1];
							else if (p[10] != FF)
								p[2] = p[10];
							else
								p[2] = S->nV++;
						} else
							p[2] = S->nV++;
						Uy[y][0] = p[2];
					}
					break;
				case 3:
					if (y || x)
						p[3] = Lz[y][x];
					else {
						if (v[0] == 0) {
							if (p[0] != FF)
								p[3] = p[0];
							else if (p[8] != FF)
								p[3] = p[8];
							else if (z && signbf(v[1]))
								p[3] = Dy[0][0];
							else if (z && signbf(v[4]))
								p[3] = Dx[0][0];
							else if (z? signbf(S->iso - F[z - 1][0][0]): 0)
								p[3] = Lz[0][0];
							else
								p[3] = S->nV++;
						} else if (v[3] == 0) {
							if (p[2] != FF)
								p[3] = p[2];
							else if (p[11] != FF)
								p[3] = p[11];
							else
								p[3] = S->nV++;
						} else
							p[3] = S->nV++;
						Lz[0][0] = p[3];
					}
					break;
				case 4:
					if (z)
						p[4] = Dy[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[4] = p[8];
							//else if (p[7] != FF)
							//	p[4] = p[7];
							else if (y && signbf(v[0]))
								p[4] = Dx[y][x];
							else if (y && signbf(v[7]))
								p[4] = Lz[y][x + 1];
							else if (y? signbf(S->iso - F[0][y - 1][di*(x + 1)]): 0)
								p[4] = Dy[y - 1][x + 1];
							else if (y && x + 1 < nx? signbf(S->iso - F[0][y][di*(x + 2)]): 0)
								p[4] = Dx[y][x + 1];
							else
								p[4] = S->nV++;
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[4] = p[5];
							else if (p[9] != FF)
								p[4] = p[9];
							else
								p[4] = S->nV++;
						} else
							p[4] = S->nV++;
						Dy[y][x + 1] = p[4];
					}
					break;
				case 5:
					if (v[5] == 0) {
						if (signbf(v[4]))
							p[5] = p[4] = (z? Dy[y][x + 1]: S->nV++);
						else if (signbf(v[1]))
							p[5] = p[9] = (z? Dx[y + 1][x]: S->nV++);
						else {
							if (z && x + 1 < nx? signbf(S->iso - F[z][y + 1][di*(x + 2)]): 0)
								p[5] = Dx[y + 1][x + 1];
							else if (z && y + 1 < ny? signbf(S->iso - F[z][y + 2][di*(x + 1)]): 0)
								p[5] = Dy[y + 1][x + 1];
							else if (z? signbf(S->iso - F[z - 1][y + 1][di*(x + 1)]): 0)
								p[5] = Lz[y + 1][x + 1]; // value of previous slice
							else
								p[5] = S->nV++;
						}
					} else
						p[5] = S->nV++;
					Lz[y + 1][x + 1] = p[5];
					break;
				case 6:
					if (v[7] == 0) {
						if (y) {
							if (signbf(v[3]))
								p[6] = p[11] = Ux[y][x];
							else if (signbf(v[4]))
								p[6] = p[7] = Lz[y][x + 1];
							else if (signbf(S->iso - F[z + 1][y - 1][di*(x + 1)]))
								p[6] = Uy[y - 1][x + 1];
							else if (x + 1 < nx? signbf(S->iso - F[z + 1][y][di*(x + 2)]): 0)
								p[6] = Ux[y][x + 1];
							else
								p[6] = S->nV++;
						} else if (p[11] != FF)
								p[6] = p[11];
							//else if (p[7] != FF)
							//	p[6] = p[7];
							else
								p[6] = S->nV++;
					} else if (v[6] == 0) {
						if (p[5] == FF)
							p[6] = (p[10] == FF? S->nV++: p[10]);
						else
							p[6] = p[5];
					} else
						p[6] = S->nV++;
					Uy[y][x + 1] = p[6];
					break;
				case 7:
					if (y)
						p[7] = Lz[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[7] = p[8];
							else if (p[4] != FF)
								p[7] = p[4];
							else if (z && signbf(v[0]))
								p[7] = Dx[0][x];
							else if (z && signbf(v[5]))
								p[7] = Dy[0][x + 1];
							else if (z && x + 1 < nx? signbf(S->iso - F[z][0][di*(x + 2)]): 0)
								p[7] = Dx[0][x + 1];
							else if (z? signbf(S->iso - F[z - 1][0][di*(x + 1)]): 0)
								p[7] = Lz[0][x + 1];
							else
								p[7] = S->nV++;
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[7] = p[6];
							else if (p[11] != FF)
								p[7] = p[11];
							else
								p[7] = S->nV++;
						} else
							p[7] = S->nV++;
						Lz[0][x + 1] = p[7];
					}
					break;
				case 8:
					if (z || y)
						p[8] = Dx[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[8] = p[3];
							else if (p[0] != FF)
								p[8] = p[0];
							else if (x && signbf(v[3]))
								p[8] = Lz[0][x];
							else if (x && signbf(v[1]))
								p[8] = Dy[0][x];
							else if (x? signbf(S->iso - F[0][0][di*(x - 1)]): 0)
								p[8] = Dx[0][x - 1];
							else
								p[8] = S->nV++;
						} else if (v[4] == 0) {
							if (p[4] != FF)
								p[8] = p[4];
							else if (p[7] != FF)
								p[8] = p[7];
							else
								p[8] = S->nV++;
						} else
							p[8] = S->nV++;
						Dx[0][x] = p[8];
					}
					break;
				case 9:
					if (z)
						p[9] = Dx[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[9] = p[0];
							//else if (p[1] != FF)
							//	p[9] = p[1];
							else if (x && signbf(v[0]))
								p[9] = Dy[y][x];
							else if (x && signbf(v[2]))
								p[9] = Lz[y + 1][x];
							else if (x? signbf(S->iso - F[0][y + 1][di*(x - 1)]): 0)
								p[9] = Dx[y + 1][x - 1];
							else
								p[9] = S->nV++;
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[9] = p[5];
							else if (p[4] != FF)
								p[9] = p[4];
							else
								p[9] = S->nV++;
						} else
							p[9] = S->nV++;
						Dx[y + 1][x] = p[9];
					}
					break;
				case 10:
					if (v[2] == 0) {
						if (x) {
							if (signbf(v[1]))
								p[10] = p[1] = Lz[y + 1][x];
							else if (signbf(v[3]))
								p[10] = p[2] = Uy[y][x];
							else if (signbf(S->iso - F[z + 1][y + 1][di*(x - 1)]))
								p[10] = Ux[y + 1][x - 1];
							else
								p[10] = S->nV++;
						} else if (p[2] != FF)
								p[10] = p[2];
							//else if (p[1] != FF)
							//	p[10] = p[1];
							else
								p[10] = S->nV++;
					} else if (v[6] == 0) {
						if (p[5] == FF)
							p[10] = (p[6] == FF? S->nV++: p[6]);
						else
							p[10] = p[5];
					} else
						p[10] = S->nV++;
					Ux[y + 1][x] = p[10];
					break;
				case 11:
					if (y)
						p[11] = Ux[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[11] = p[3];
							else if (p[2] != FF)
								p[11] = p[2];
							else if (x && signbf(v[0]))
								p[11] = Lz[0][x];
							else if (x && signbf(v[2]))
								p[11] = Uy[0][x];
							else if (x? signbf(S->iso - F[z + 1][0][di*(x - 1)]): 0)
								p[11] = Ux[0][x - 1];
							else
								p[11] = S->nV++;
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[11] = p[6];
							else if (p[7] != FF)
								p[11] = p[7];
							else
								p[11] = S->nV++;
						} else
							p[11] = S->nV++;
						Ux[0][x] = p[11];
					}
				break;
				default:
					p[12] = S->nV++;
				}
			}
			ti[--k] = p[c];
		}
		if (ti[0] != ti[1] && ti[0] != ti[2] && ti[1] != ti[2])
			++S->nT;
	}
}
#undef FF

void MC33::free_temp_D_U() {
	delete[] Dx; delete[] Ux; delete[] Dy; delete[] Uy;
	delete[] Lz;
}

void MC33::clear_temp_isosurface() {
	if (F) {
		for (unsigned int y = 0; y != ny; y++) {
			delete[] Dx[y]; delete[] Ux[y]; delete[] Dy[y]; delete[] Uy[y];
			delete[] Lz[y];
		}
		delete[] Dx[ny]; delete[] Ux[ny];
		delete[] Lz[ny];
		free_temp_D_U();
		F = 0;
	}
}

#define code_in_define_part01\
			unsigned int nv = S->nV++;\
			if (!(nv&0x0FFF)) {\
				try {\
					S->V.resize(nv + 0x1000);\
					S->N.resize(nv + 0x1000);\
				}\
				catch (...) {\
					memoryfault = 1;\
					S->nV = (nv? 1: 0);\
					return 0u;\
				}\
			}

#ifndef MC33_NORMAL_NEG
#define code_in_define_part02\
			MC33_real t = invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);\
			float *q = S->N[nv].v;\
			*q = t * *r; *(++q) = t * *(++r); *(++q) = t * *(++r);\
			return nv;
#else //MC33_NORMAL_NEG reverses the direction of the normal
#define code_in_define_part02\
			MC33_real t = -invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);\
			float *q = S->N[nv].v;\
			*q = t * *r; *(++q) = t * *(++r); *(++q) = t * *(++r);\
			return nv;
#endif

int MC33::set_grid3d(grid3d &G) {
	clear_temp_isosurface();
	nx = G.N[0];
	ny = G.N[1];
	nz = G.N[2];
#ifndef GRD_ORTHOGONAL
	if (G.nonortho) {
		store_point = [this] (MC33_real* r) {
			code_in_define_part01
			multAbf(_A,r,r,0);
			for (int i = 0; i != 3; i++)
				S->V[nv].v[i] = *(r++) + MC_O[i];
			multAbf(A_,r,r,1);
			code_in_define_part02
		};
		for (int j = 0; j != 3; j++)
			for (int i = 0; i != 3; i++) {
				_A[j][i] = G._A[j][i]*G.d[i]; // true transformation matrices
				A_[j][i] = G.A_[j][i]/G.d[j];
			}
	} else
#endif
	if (G.d[0] != G.d[1] || G.d[1] != G.d[2]) {
		ca = G.d[2]/G.d[0];
		cb = G.d[2]/G.d[1];
		store_point = [this] (MC33_real* r) {
			code_in_define_part01
			for (int i = 0; i != 3; i++)
				S->V[nv].v[i] = *(r++)*MC_D[i] + MC_O[i];
			r[0] *= ca; // normal[0]
			r[1] *= cb; // normal[1]
			code_in_define_part02
		};
	} else if (G.d[0] == 1 && G.r0[0] == 0 && G.r0[1] == 0 && G.r0[2] == 0)
		store_point = [this] (MC33_real* r) {
			code_in_define_part01
			for (int i = 0; i != 3; i++)
				S->V[nv].v[i] = *(r++);
			code_in_define_part02
		};
	else
		store_point = [this] (MC33_real* r) {
			code_in_define_part01
			for (int i = 0; i != 3; i++)
				S->V[nv].v[i] = *(r++)*MC_D[i] + MC_O[i];
			code_in_define_part02
		};

	for (int j = 0; j != 3; j++) {
		MC_O[j] = G.r0[j];
		MC_D[j] = G.d[j];
	}
	F = const_cast<const GRD_data_type***>(G.F);
	if (!F)
		return 0;
	Lz = new (std::nothrow) unsigned int*[ny + 1]; //edges 1, 3, 5 (only write) and 7
	Dy = new (std::nothrow) unsigned int*[ny]; //edges 0 and 4
	Uy = new (std::nothrow) unsigned int*[ny]; //edges 2 and 6 (only write)
	Dx = new (std::nothrow) unsigned int*[ny + 1]; //edges 8 and 9
	Ux = new (std::nothrow) unsigned int*[ny + 1]; //edges 10 (only write) and 11
	if (!Ux) {
		free_temp_D_U();
		return -1;
	}
	for (unsigned int j = 0; j != ny; j++) {
		Dx[j] = new (std::nothrow) unsigned int[nx];
		Ux[j] = new (std::nothrow) unsigned int[nx];
		Lz[j] = new (std::nothrow) unsigned int[nx + 1];
		Dy[j] = new (std::nothrow) unsigned int[nx + 1];
		Uy[j] = new (std::nothrow) unsigned int[nx + 1];
	}
	if (Uy[ny - 1]) {
		Dx[ny] = new (std::nothrow) unsigned int[nx];
		Ux[ny] = new (std::nothrow) unsigned int[nx];
		Lz[ny] = new (std::nothrow) unsigned int[nx + 1];
		di = G.x_data > 1? G.x_data: 1;
		if (Lz[ny])
			return 0;
	} else
		Dx[ny] = Ux[ny] = Lz[ny] = 0;
	clear_temp_isosurface();
	return -1;
}
#undef code_in_define_part01
#undef code_in_define_part02

int MC33::set_grid3d(grid3d *G) {
	if (!G)
		return -1;
	return set_grid3d(*G);
}

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-pointer" 
#endif

int MC33::calculate_isosurface(surface &Sf, MC33_real iso) {
	const GRD_data_type ***FG = F;
	if (!FG)//The set_grid3d function was not executed
		return -2;
	if (Sf.nV)
		Sf.clear();
	S = &Sf;
	Sf.iso = iso;
	memoryfault = 0;
	unsigned int d = di, Nx = nx;
	MC33_real Vt[12];
	MC33_real *v2 = Vt;
	v = Vt + 4;
	for (unsigned int z = 0; z != nz; z++) {
		const GRD_data_type **F0 = *FG;
		const GRD_data_type **F1 = *(++FG);
		for (unsigned int y = 0; y != ny; y++) {
			const GRD_data_type *V00 = *F0;
			const GRD_data_type *V01 = *(++F0);
			const GRD_data_type *V10 = *F1;
			const GRD_data_type *V11 = *(++F1);
			v2[0] = iso - *V00;//the difference was inverted to use the signbit function
			v2[1] = iso - *V01;
			v2[2] = iso - *V11;
			v2[3] = iso - *V10;
			//the eight least significant bits of i correspond to vertex indices. (x...x01234567)
			//If the bit is 1 then the vertex value is greater than zero.
			unsigned int i = signbf(v2[3]) != 0;
			if (signbf(v2[0])) i |= 8;
			if (signbf(v2[1])) i |= 4;
			if (signbf(v2[2])) i |= 2;
			for (unsigned int x = 0; x != Nx; x++) {
				std::swap(v2, v);
				V00 += d;
				v2[0] = iso - *V00;//*(++V00);
				V01 += d;
				v2[1] = iso - *V01;//*(++V01);
				V11 += d;
				v2[2] = iso - *V11;//*(++V11);
				V10 += d;
				v2[3] = iso - *V10;//*(++V10);
				i = ((i&0x0F)<<4)|(signbf(v2[3]) != 0);
				if (signbf(v2[0])) i |= 8;
				if (signbf(v2[1])) i |= 4;
				if (signbf(v2[2])) i |= 2;
				if (i && i^0xFF) { //i is different from 0 and 0xFF
					if (v2 < v) { MC33_real *t = v2; MC33_real *s = t + 8; *s = *t;
						*(++s) = *(++t); *(++s) = *(++t); *(++s) = *(++t); }
					find_case(x, y, z, i);
				}
			}
		}
		std::swap(Dx, Ux);
		std::swap(Dy, Uy);
	}
	try {
		Sf.color = std::vector<int>(Sf.nV, DefaultColor);
	}
	catch (...) {
		memoryfault = 1;
	}
	if (memoryfault) {
		Sf.clear();
		return -1;
	}
	return 0;
}

size_t MC33::size_of_isosurface(MC33_real iso, unsigned int &nV, unsigned int &nT) {
	const GRD_data_type ***FG = F;
	if (!FG)//The set_grid3d function was not executed
		return -2;
	surface Sf;
	S = &Sf;
	Sf.iso = iso;
	memoryfault = 0;
	unsigned int d = di, Nx = nx;
	MC33_real Vt[12];
	MC33_real *v2 = Vt;
	v = Vt + 4;
	for (unsigned int z = 0; z != nz; z++) {
		const GRD_data_type **F0 = *FG;
		const GRD_data_type **F1 = *(++FG);
		for (unsigned int y = 0; y != ny; y++) {
			const GRD_data_type *V00 = *F0;
			const GRD_data_type *V01 = *(++F0);
			const GRD_data_type *V10 = *F1;
			const GRD_data_type *V11 = *(++F1);
			v2[0] = iso - *V00;
			v2[1] = iso - *V01;
			v2[2] = iso - *V11;
			v2[3] = iso - *V10;
			unsigned int i = signbf(v2[3]) != 0;
			if (signbf(v2[2])) i |= 2;
			if (signbf(v2[1])) i |= 4;
			if (signbf(v2[0])) i |= 8;
			for (unsigned int x = 0; x != Nx; x++) {
				std::swap(v, v2);
				V00 += d;
				v2[0] = iso - *V00;//*(++V00);
				V01 += d;
				v2[1] = iso - *V01;//*(++V01);
				V11 += d;
				v2[2] = iso - *V11;//*(++V11);
				V10 += d;
				v2[3] = iso - *V10;//*(++V10);
				i = ((i&0x0F)<<4)|(signbf(v2[3]) != 0);
				if (signbf(v2[0])) i |= 8;
				if (signbf(v2[1])) i |= 4;
				if (signbf(v2[2])) i |= 2;
				if (i && i^0xFF) {
					if (v > v2) { MC33_real *t = v2; MC33_real *s = t + 8; *s = *t;
						*(++s) = *(++t); *(++s) = *(++t); *(++s) = *(++t); }
					case_count(x, y, z, i);
				}
			}
		}
		std::swap(Dx, Ux);
		std::swap(Dy, Uy);
	}
	nV = Sf.nV;
	nT = Sf.nT;
	// number of vertices * (size of vertex and normal + size of color ) +
	// number of triangle * size of triangle + size of class surface
	return nV*(3*(sizeof(MC33_real) + sizeof(float)) + sizeof(int)) + nT*(3*sizeof(int)) + sizeof(surface);
}
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

size_t MC33::size_of_isosurface(MC33_real iso) {
	unsigned int nV, nT;
	return size_of_isosurface(iso, nV, nT);
}

//******************************************************************
const unsigned int* grid3d::get_N() {return N;}
const float* grid3d::get_L() {return L;}
const double* grid3d::get_r0() {return r0;}
const double* grid3d::get_d() {return d;}
const char* grid3d::get_title() {return title;}

GRD_data_type grid3d::get_grid_value(unsigned int i, unsigned int j, unsigned int k) {
	if (F && i <= N[0] && j <= N[1] && k <= N[2])
		return F[k][j][x_data > 1? i*x_data: i];
	return sqrt(-1); // NaN value
}

GRD_data_type grid3d::interpolated_value(double x, double y, double z) {
	double r[3] = {x - r0[0], y - r0[1], z - r0[2]};
#ifndef GRD_ORTHOGONAL
	if(nonortho)
		multAb(A_,r,r,0);
#endif
	return (this->*interpolation)(r);
}

GRD_data_type grid3d::bad_value(double*) const {
	return sqrt(-1); // NaN value
}

int grid3d::set_interpolation(int i) {
	if (!F || x_data > 0)
		i = 0;
	switch (i) {
		case 1:
			interpolation = &grid3d::trilinear;
			break;
		case 3:
			interpolation = &grid3d::tricubic;
			break;
		default:
			i = 0;
			interpolation = &grid3d::bad_value;
			break;
	}
	return (i? 0: -1);
}


GRD_data_type grid3d::trilinear(double *r) const {
	unsigned int i[3];
	double t[3], s[3];
	for (int j = 0; j < 3; j++) {
		if ((periodic>>j)&0x01) {
			t[j] = (r[j] - L[j]*floor(r[j]/L[j]))/d[j];
			i[j] = (int)t[j];
			t[j] -= i[j];
			i[j] = (i[j]? i[j] - 1: N[j] - 1);
		} else {
			t[j] = r[j]/d[j];
			if (t[j] < 0 || t[j] > N[j])
				return sqrt(-1);
			i[j] = (int)t[j];
			t[j] -= i[j];
			if (i[j] == N[j]) {
				t[j] = 1.0;
				--i[j];
			}
		}
		s[j] = 1.0 - t[j];
	}
	return s[0]*s[1]*s[2]*F[i[2]][i[1]][i[0]] + t[0]*s[1]*s[2]*F[i[2]][i[1]][i[0] + 1] +
		s[0]*t[1]*s[2]*F[i[2]][i[1] + 1][i[0]] + s[0]*s[1]*t[2]*F[i[2] + 1][i[1]][i[0]] +
		t[0]*t[1]*s[2]*F[i[2]][i[1] + 1][i[0] + 1] + s[0]*t[1]*t[2]*F[i[2] + 1][i[1] + 1][i[0]] +
		t[0]*s[1]*t[2]*F[i[2] + 1][i[1]][i[0] + 1] + t[0]*t[1]*t[2]*F[i[2] + 1][i[1] + 1][i[0] + 1];
}

//from https://facyt-quimicomp.neocities.org/Vega_en.html#c_library
GRD_data_type grid3d::tricubic(double *r) const {
	unsigned int iR[3], x, y, z, i, j, k;
	double f = 0, t, w[3][4];
	for (i = 0; i < 3; i++)
	{
		if ((periodic>>i)&0x01)
		{
			t = (r[i] - L[i]*floor(r[i]/L[i]))/d[i];
			iR[i] = (int)t;
			t -= iR[i];
			iR[i] = (iR[i]? iR[i] - 1: N[i] - 1);
		}
		else
		{
			t = r[i]/d[i];
			if (t < 0 || t > N[i])
				return sqrt(-1);
			iR[i] = (int)t;
			t -= iR[i];
			if (!iR[i])
				t--;
			else if (iR[i] > N[i] - 2)
			{
				z = iR[i] - N[i] + 2;
				iR[i] -= z + 1;
				t += z;
			}
			else
				iR[i]--;
		}
		w[i][0] = (-2.0 + (3.0 - t)*t)*t*(1.0/6.0);
		w[i][1] = (2.0 + (-1.0 + (-2.0 + t)*t)*t)*0.5;
		w[i][2] = (2.0 + (1.0 - t)*t)*t*0.5;
		w[i][3] = (t*t - 1.0)*t*(1.0/6.0);
	}
	z = iR[2];
	for (k = 0; k < 4; k++)
	{
		y = iR[1];
		for (j = 0; j < 4; j++)
		{
			x = iR[0];
			for (i = 0; i < 4; i++)
			{
				f += w[0][i]*w[1][j]*w[2][k]*F[z][y][x];
				if ((++x) > N[0]) x = 1;
			}
			if ((++y) > N[1]) y = 1;
		}
		if ((++z) > N[2]) z = 1;
	}
	return f;
}


#ifndef GRD_ORTHOGONAL
const float* grid3d::get_Ang() {return Ang;}
const double (*grid3d::get__A())[3] {return _A;}
const double (*grid3d::get_A_())[3] {return A_;}
int grid3d::isnotorthogonal() {return nonortho;}

void grid3d::update_matrices() {
	if (Ang[0] != 90 || Ang[1] != 90 || Ang[2] != 90) {
		nonortho = 1;
		double ca = cos(Ang[0]*(M_PI/180.0));
		double cb = cos(Ang[1]*(M_PI/180.0));
		double aux1 = Ang[2]*(M_PI/180.0);
		double sg = sin(aux1);
		double cg = cos(aux1);
		aux1 = ca - cb*cg;
		double aux2 = sqrt(sg*sg + 2*ca*cb*cg - ca*ca - cb*cb);
		_A[0][0] = A_[0][0] = 1.0;
		_A[0][1] = cg;
		_A[0][2] = cb;
		_A[1][1] = sg;
		A_[1][1] = cb = 1.0/sg;
		A_[0][1] = -cg*cb;
		_A[1][2] = aux1*cb;
		_A[2][2] = aux2*cb;
		aux2 = 1.0/aux2;
		A_[0][2] = (cg*aux1 - ca*sg*sg)*cb*aux2;
		A_[1][2] = -aux1*cb*aux2;
		A_[2][2] = sg*aux2;
		_A[1][0] = _A[2][0] = _A[2][1] = 0.0;
		A_[1][0] = A_[2][0] = A_[2][1] = 0.0;
		multAbf = mult_TSAbf;
	} else {
		nonortho = 0;
		setIdentMat3x3d(_A);
		setIdentMat3x3d(A_);
	}
}

void grid3d::set_Ang(float angle_bc, float angle_ca, float angle_ab) {
	Ang[0] = angle_bc; Ang[1] = angle_ca; Ang[2] = angle_ab;
	update_matrices();
}
#endif

void grid3d::set_ratio_aspect(double rx, double ry, double rz) {
	if (F) {
		d[0] = rx; d[1] = ry; d[2] = rz;
		for (int i = 0; i != 3; i++)
			L[i] = N[i]*d[i];
	}
}

void grid3d::set_r0(double x, double y, double z) {
	r0[0] = x; r0[1] = y; r0[2] = z;
}

void grid3d::set_title(const char *s) {
	char *t = title, *te = title + sizeof(title) - 1;
	while (*s && t < te)
		*(t++) = *(s++);
	*t = (char)0;
}

void grid3d::set_grid_value(unsigned int i, unsigned int j, unsigned int k, GRD_data_type value) {
	if (!x_data && F && i <= N[0] && j <= N[1] && k <= N[2])
		F[k][j][i] = value;
}

//******************************************************************
void grid3d::free_F() {
	if (F) {
		clear_subgrid();
		if (x_data)
			for (unsigned int k = 0; k <= N[2]; k++)
				delete[] F[k];
		else
			for (unsigned int k = 0; k <= N[2]; k++) {
				if (F[k]) {
					for (unsigned int j = 0; j <= N[1]; j++)
						delete[] F[k][j];
				} else 
					break;
				delete[] F[k];
			}
		delete[] F;
		F = 0;
		interpolation = &grid3d::bad_value;
	}
}

grid3d::grid3d() : F(0), subgrid(0), nsg(0), maxnsg(0) {}

grid3d::grid3d(const grid3d &G) {
	if (!G.F)
		return;
	memcpy(this, &G, sizeof (grid3d));
	F = 0;
	subgrid = 0;
	nsg = maxnsg = 0;
	if (alloc_F())
		return;
	for (unsigned int k = 0; k <= N[2]; k++)
		for (unsigned int j = 0; j <= N[1]; j++) {
			if (G.x_data > 1)
				for (unsigned int i = 0; i <= N[0]; i++)
					F[k][j][i] = G.F[k][j][i*G.x_data];
			else
				memcpy(F[k][j], G.F[k][j], (N[0] + 1)*sizeof(GRD_data_type));
		}
}

grid3d::~grid3d() {
	free_F();
}

int grid3d::alloc_F() {
	unsigned int j, k;
	F = new (std::nothrow) GRD_data_type**[N[2] + 1];
	if (!F)
		return -1;
	for (k = 0; k <= N[2]; k++) {
		F[k] = new (std::nothrow) GRD_data_type*[N[1] + 1];
		if (!F[k])
			return -1;
		for (j = 0; j <= N[1]; j++) {
			F[k][j] = new (std::nothrow) GRD_data_type[N[0] + 1];
			if (!F[k][j]) {
				while (j)
					delete[] F[k][--j];
				delete[] F[k];
				F[k] = 0;
				return -1;
			}
		}
	}
	x_data = 0;
	interpolation = &grid3d::trilinear;
	return 0;
}

void grid3d::delete_grid_data() {
	free_F();
	N[0] = N[1] = N[2] = 0;
}

int grid3d::set_grid_dimensions(unsigned int Nx, unsigned int Ny, unsigned int Nz) {
	free_F();
	if (Nx == 0 || Ny == 0 || Nz == 0) {
		N[0] = N[1] = N[2] = 0;
		return -1;
	}
	N[0] = Nx - 1;
	N[1] = Ny - 1;
	N[2] = Nz - 1;
	if (alloc_F()) {
		free_F();
		return -2;
	}
	for (int i = 0; i != 3; i++) {
		L[i] = N[i];
		d[i] = 1.0;
		r0[i] = 0.0;
#ifndef GRD_ORTHOGONAL
		Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_ORTHOGONAL
	nonortho = 0;
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return 0;
}

int grid3d::set_data_pointer(unsigned int Nx, unsigned int Ny, unsigned int Nz, GRD_data_type* data) {
	free_F();
	if (!data || Nx == 0 || Ny == 0 || Nz == 0) {
		N[0] = N[1] = N[2] = 0;
		return -1;
	}
	N[0] = Nx - 1;
	N[1] = Ny - 1;
	N[2] = Nz - 1;
	F = new (std::nothrow) GRD_data_type**[Nz];
	if (!F)
		return -2;
	for (unsigned int k = 0; k < Nz; k++) {
		F[k] = new (std::nothrow) GRD_data_type*[Ny];
		if (!F[k]) {
			while (k)
				delete[] F[--k];
			delete[] F;
			F = 0;
			return -2;
		}
		for (unsigned int j = 0; j < Ny; j++)
			F[k][j] = data + j*Nx;
		data += Ny*Nx;
	}
	for (int i = 0; i != 3; i++) {
		L[i] = N[i];
		d[i] = 1.0;
		r0[i] = 0.0;
#ifndef GRD_ORTHOGONAL
		Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_ORTHOGONAL
	nonortho = 0;
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	x_data = -1;
	return 0;
}

void grid3d::clear_subgrid() {
	if (subgrid) {
		for (unsigned int i = 0; i != nsg; i++)
			delete subgrid[i];
		nsg = maxnsg = 0;
		delete[] subgrid;
		subgrid = 0;
	}
}

grid3d *grid3d::get_subgrid(unsigned int i) {
	return (i < nsg? subgrid[i]: 0);
}

void grid3d::del_subgrid(unsigned int i) {
	if (i < nsg) {
		delete subgrid[i];
		while (++i < nsg)
			subgrid[i - 1] = subgrid[i];
		--nsg;
	}
}

unsigned int grid3d::subgrid_size() {
	return nsg;
}

int grid3d::add_subgrid(unsigned int Oi, unsigned int Oj, unsigned int Ok,
					unsigned int Ni, unsigned int Nj, unsigned int Nk,
					unsigned int Si, unsigned int Sj, unsigned int Sk) {
	if (!F || x_data > 0 || !Ni || !Nj || !Nk || !Si || !Sj || !Sk)
		return -1;
	if (Oi + Si*(Ni - 1) > N[0] || Oj + Sj*(Nj - 1) > N[1] || Ok + Sk*(Nk - 1) > N[2])
		return -1;
	if (nsg == maxnsg) {
			grid3d **t = new (std::nothrow) grid3d*[maxnsg + 8];
			if (!t)
				return -2;
			if (subgrid) {
				memcpy(t, subgrid, nsg*sizeof(void*));
				delete[] subgrid;
		}
		subgrid = t;
		maxnsg += 8;
	}
	grid3d *G = new (std::nothrow) grid3d();
	if (!G)
		return -2;

	G->N[0] = Ni - 1;
	G->N[1] = Nj - 1;
	G->N[2] = Nk - 1;
	G->F = new (std::nothrow) GRD_data_type**[Nk];
	if (!G->F) {
		delete G;
		return -2;
	}
	for (unsigned int k = 0; k < Nk; k++) {
		G->F[k] = new (std::nothrow) GRD_data_type*[Nj];
		if (!G->F[k]) {
			while (++k < Nk)
				G->F[k] = 0;
			delete G;
			return -2;
		}
		for (unsigned int j = 0; j < Nj; j++)
			G->F[k][j] = F[Ok + k*Sk][Oj + j*Sj] + Oi;
	}
	G->d[0] = Si*d[0];
	G->d[1] = Sj*d[1];
	G->d[2] = Sk*d[2];
	G->title[0] = 0;
	double d0[3] = {Oi*d[0], Oj*d[1], Ok*d[2]};
#ifndef GRD_ORTHOGONAL
	memcpy(G->A_, A_, sizeof A_);
	memcpy(G->_A, _A, sizeof _A);
	memcpy(G->Ang, Ang, sizeof Ang);
	G->nonortho = nonortho;
	if (nonortho)
		multAb(_A, d0, d0, 0);
#endif
	for (int i = 0; i != 3; i++) {
		G->r0[i] = r0[i] + d0[i];
		G->L[i] = G->N[i]*G->d[i];
	}

	G->x_data = Si;
	G->periodic = 0;
	if (Si*Ni == N[0] + 1)
		G->periodic = periodic&1;
	if (Sj*Nj == N[1] + 1)
		G->periodic |= periodic&2;
	if (Sk*Nk == N[2] + 1)
		G->periodic |= periodic&4;
	subgrid[nsg++] = G;
	interpolation = &grid3d::bad_value;
	return 0;
}

int grid3d::generate_grid_from_fn(double xi, double yi, double zi, double xf, double yf, double zf,
	double dx, double dy, double dz, double (*fn)(double x, double y, double z)) {
	free_F();
	if (dx <= 0 || dy <= 0 || dz <= 0 || xi == xf || yi == yf || zi == zf) {
		N[0] = N[1] = N[2] = 0;
		return -1;
	}
	if (xi > xf)
		std::swap(xi, xf);
	if (xf - xi < dx)
		dx = xf - xi;
	if (yi > yf)
		std::swap(yi, yf);
	if (yf - yi < dy)
		dy = yf - yi;
	if (zi > zf)
		std::swap(zi, zf);
	if (zf - zi < dz)
		dz = zf - zi;
	N[0] = int((xf - xi)/dx + 0.5);
	N[1] = int((yf - yi)/dy + 0.5);
	N[2] = int((zf - zi)/dz + 0.5);
	if (alloc_F()) {
		free_F();
		return -2;
	}
	d[0] = dx; d[1] = dy; d[2] = dz;
	r0[0] = xi; r0[1] = yi; r0[2] = zi;
	if (fn) {
		double x, y, z = zi;
		for (size_t k = 0; k <= N[2]; k++) {
			y = yi;
			for (size_t j = 0; j <= N[1]; j++) {
				x = xi;
				for (size_t i = 0; i <= N[0]; i++) {
					F[k][j][i] = (GRD_data_type)fn(x,y,z);
					x += dx;
				}
				y += dy;
			}
			z += dz;
		}
	}
	for (int i = 0; i != 3; i++) {
		L[i] = N[i]*d[i];
#ifndef GRD_ORTHOGONAL
		Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_ORTHOGONAL
	nonortho = 0;
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return 0;
}


//******************************************************************
int grid3d::read_grd(const char *filename) {
	char cs[32];
	unsigned int i, j, k;
	int xi[3], order;

	free_F();
	std::ifstream in(filename);
	if (!in)
		return -1;
	nonortho = 0;
	periodic = 0;
	in.getline(title, 159);
	in.ignore(60, '\n');
#ifdef GRD_ORTHOGONAL
	float Ang[3];
#endif
	in >> L[0] >> L[1] >> L[2] >> Ang[0] >> Ang[1] >> Ang[2] >> N[0] >> N[1] >> N[2];
	in >> order >> xi[0] >> cs >> xi[1] >> cs >> xi[2];
	in.ignore(20, '\n');
	if (N[0] < 1 || N[1] < 1 || N[2] < 1) return -1;
	if (order != 1 && order != 3) return -1;
	for (i = 0; i != 3; i++) {
		d[i] = L[i]/N[i];
		r0[i] = xi[i]*d[i];
	}
	periodic = (xi[0] == 0)|((xi[1] == 0)<<1)|((xi[2] == 0)<<2);

#ifndef GRD_ORTHOGONAL
	update_matrices();
#endif

	if (alloc_F()) {
		free_F();
		return -2;
	}
	for (k = 0; k <= N[2]; k++)
		if (order == 1) {
			for (j = 0; j <= N[1]; j++)
				for (i = 0; i <= N[0]; i++)
				{
					in.getline(cs, 31);
					F[k][j][i] = std::stof(cs);
				}
		} else {
			for (i = 0; i <= N[0]; i++)
				for (j = 0; j <= N[1]; j++)
				{
					in.getline(cs, 31);
					F[k][j][i] = std::stof(cs);
				}
		}
	return (in.good()? 0: -4);
}

int grid3d::read_grd_binary(const char* filename) {
	unsigned int i;
	free_F();
	std::ifstream in(filename, std::ios::binary);
	if (!in)
		return -1;
	in.read(reinterpret_cast<char*>(&i), sizeof(int));
	if (i != 0x4452475f) // _GRD
		return -1;

	in.read(reinterpret_cast<char*>(&i), sizeof(int));
	if (i > 159)
		return -1;
	in.read(title, i);
	in.read(reinterpret_cast<char*>(N), sizeof N);
	in.read(reinterpret_cast<char*>(L), sizeof L);
	in.read(reinterpret_cast<char*>(r0), sizeof r0);
	in.read(reinterpret_cast<char*>(d), sizeof d);
#ifndef GRD_ORTHOGONAL
	in.read(reinterpret_cast<char*>(&nonortho), sizeof(int));
	if (nonortho) {
		in.read(reinterpret_cast<char*>(Ang),3*sizeof(float));
		in.read(reinterpret_cast<char*>(_A), sizeof _A);
		in.read(reinterpret_cast<char*>(A_), sizeof A_);
		multAbf = mult_Abf;
	} else {
		setIdentMat3x3d(_A);
		setIdentMat3x3d(A_);
	}
#else
	in.read(reinterpret_cast<char*>(&i), sizeof(int));
	if (i)
		in.ignore(3*sizeof(float) + 18*sizeof(double));
#endif
	periodic = (r0[0] == 0 && r0[1] == 0 && r0[2] == 0? 7: 0);
	if (alloc_F()) {
		free_F();
		return -2;
	}
	for (unsigned int k = 0; k <= N[2]; k++)
		for (unsigned int j = 0; j <= N[1]; j++)
#if GRD_TYPE_SIZE == 8
			for (unsigned int i = 0; i <= N[0]; i++) {
				float f;
				in.read(reinterpret_cast<char*>(&f), sizeof(float));
				F[k][j][i] = f;
			}
#else
			in.read(reinterpret_cast<char*>(F[k][j]),(N[0] + 1)*sizeof(float));
#endif
	return (in.good()? 0: -4);
}
//******************************************************************
int grid3d::read_scanfiles(const char *filename, unsigned int res, int order) {
	std::string nm(filename);
	std::ifstream in;
	unsigned int i, j, l;
	int m, k = -1;
	unsigned short int n;
	free_F();
#ifndef GRD_ORTHOGONAL
	nonortho = 0;
#endif
	periodic = 0;
	r0[0] = r0[1] = r0[2] = 0.0;
	for (i = 0; i != 2; i++) {
		d[i] = 1.0;
		N[i] = res - 1;
		L[i] = float(res - 1);
	}
	d[2] = 1.0;
	l = unsigned(nm.size() - 1);
	while (nm[l] >= '0' && nm[l] <= '9')
		l--;
	m = stoi(nm.substr(++l));
	while (1) {
		nm.replace(l,6,std::to_string(m++));
		in.open(nm, std::ios::binary|std::ifstream::in);
		if (!in) {
			m = 0;
			break;
		}
		if (!((++k)&63)) {
			GRD_data_type ***pt = new (std::nothrow) GRD_data_type**[k + 64];
			if (!pt) {
				m = -2;
				break;
			}
			if (F) {
				memcpy(pt, F, k*sizeof(void*));
				delete[] F;
			}
			F = pt;
		}
		F[k] = new (std::nothrow) GRD_data_type*[res];
		if (!F[k]) {
			m = -2;
			break;
		}
		for (j = 0; j != res; j++)
			F[k][j] = new (std::nothrow) GRD_data_type[res];
		if (!F[k][N[1]]) {
			m = -2;
			break;
		}

		if (order)
			for (j = 0; j != res; j++)
				for (i = 0; i != res; i++) {
					in.read(reinterpret_cast<char*>(&n), sizeof(short int));
					F[k][j][i] = static_cast<unsigned short int>((n>>8)|(n<<8));
				}
		else
			for (j = 0; j != res; j++)
				for (i = 0; i != res; i++) {
					in.read(reinterpret_cast<char*>(&n), sizeof(short int));
					F[k][j][i] = n;
				}
		if (in.fail()) {
			for (j = 0; j != res; j++)
				delete[] F[k][j];
			delete[] F[k];
			m = (--k > 0? 0: -4);//m = -4;
			break;
		}
		in.close();
	}
	N[2] = k;
	L[2] = float(k);
	x_data = 0;
	interpolation = &grid3d::trilinear;
	if (m == -2)
		free_F();
	else {
		j = (k + 1)/2;
		for (i = 0; i != j; i++)
			std::swap(F[i], F[k - i]);
	}
#ifndef GRD_ORTHOGONAL
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return m;
}

//******************************************************************
int grid3d::read_raw_file(const char *filename, unsigned int *n, int byte, int isfloat) {
	unsigned int i, j, k;
	free_F();
	std::ifstream in(filename, std::ios::binary);
	if (!in)
		return -1;
	unsigned int ui = 0;
	if (isfloat) {
		if (byte != sizeof(float) && byte != sizeof(double))
			return -1;
	} else if (abs(byte) > 4 || abs(byte) == 3 || !byte)
		return -1;
	if (byte == -1) byte = 1;
	if (n[0] == 0 || n[1] == 0 || n[2] == 0)
		return -1; // bad input parameters
#ifndef GRD_ORTHOGONAL
	nonortho = 0;
#endif
	periodic = 0;
	r0[0] = r0[1] = r0[2] = 0.0;
	for (int h = 0; h != 3; ++h) {
		d[h] = 1.0;
		N[h] = n[h] - 1;
		L[h] = float(N[h]);
	}
	if (alloc_F()) {
		free_F();
		return -2;
	}
#ifdef GRD_INTEGER
	if (!isfloat && GRD_TYPE_SIZE == byte)
#else
	if (isfloat && GRD_TYPE_SIZE == byte)
#endif
	{
		byte *= n[0];
		for (k = 0; k != n[2]; k++)
			for (j = 0; j != n[1]; j++)
				in.read(reinterpret_cast<char*>(F[k][j]), byte);
	} else if (isfloat) {
		if (byte == 8) {
#if defined (GRD_INTEGER) || GRD_TYPE_SIZE == 4
			double df;
			for (k = 0; k != n[2]; k++)
				for (j = 0; j != n[1]; j++)
					for (i = 0; i != n[0]; i++)
					{
						in.read(reinterpret_cast<char*>(&df), byte);
						F[k][j][i] = GRD_data_type(df);
					}
#endif
		} else {
#if defined (GRD_INTEGER) || GRD_TYPE_SIZE == 8
			float f;
			for (k = 0; k != n[2]; k++)
				for (j = 0; j != n[1]; j++)
					for (i = 0; i != n[0]; i++) {
						in.read(reinterpret_cast<char*>(&f), byte);
						F[k][j][i] = GRD_data_type(f);
					}
#endif
		}
	} else if (byte < 0) {
		for (k = 0; k != n[2]; k++)
			for (j = 0; j != n[1]; j++)
				for (i = 0; i != n[0]; i++) {
					in.read(reinterpret_cast<char*>(&ui), -byte);
					if (byte == -2)
						F[k][j][i] = GRD_data_type((ui>>8)|(ui<<8));
					else if (byte == -4)
						F[k][j][i] = GRD_data_type((ui>>24)|((ui>>8)&0x00f0)|((ui<<8)&0x0f00)|(ui<<24));
					else
						F[k][j][i] = GRD_data_type((ui>>16)|(ui&0x00f0)|(ui<<16));
				}
	} else {
		for (k = 0; k != n[2]; k++)
			for (j = 0; j != n[1]; j++)
				for (i = 0; i != n[0]; i++) {
					in.read(reinterpret_cast<char*>(&ui), byte);
					F[k][j][i] = GRD_data_type(ui);
				}
	}
#ifndef GRD_ORTHOGONAL
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return (in.good()? 0: -4);
}

//******************************************************************
int grid3d::read_dat_file(const char *filename) {
	unsigned short int n, nx, ny, nz;
	free_F();
	std::ifstream in(filename, std::ios::binary);
	if (!in)
		return -1;
#ifndef GRD_ORTHOGONAL
	nonortho = 0;
#endif
	periodic = 0;
	r0[0] = r0[1] = r0[2] = 0.0;
	for (int h = 0; h != 3; ++h)
		d[h] = 1.0;
	in.read(reinterpret_cast<char*>(&nx), sizeof(short int));
	in.read(reinterpret_cast<char*>(&ny), sizeof(short int));
	in.read(reinterpret_cast<char*>(&nz), sizeof(short int));
	N[0] = nx - 1;
	N[1] = ny - 1;
	N[2] = nz - 1;
	for (int h = 0; h != 3; ++h)
		L[h] = float(N[h]);
	if (alloc_F()) {
		free_F();
		return -2;
	}
	while (nz--)
		for (unsigned int j = 0; j < ny; j++)
			for (unsigned int i = 0; i < nx; i++) {
				in.read(reinterpret_cast<char*>(&n), sizeof(short int));
				F[nz][j][i] = n;
			}
#ifndef GRD_ORTHOGONAL
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return (in.good()? 0: -4);
}

//******************************************************************
int grid3d::save_raw_file(const char *filename) {
	if (!F)
		return -1;
	std::ofstream out(filename, std::ios::binary);
	if (!out)
		return -1;
	for (unsigned int k = 0; k <= N[2]; k++)
		for (unsigned int j = 0; j <= N[1]; j++) {
			if (x_data > 1)
				for (unsigned int i = 0; i <= N[0]; i++)
					out.write(reinterpret_cast<char*>(F[k][j] + i*x_data), sizeof(GRD_data_type));
			else
				out.write(reinterpret_cast<char*>(F[k][j]), (N[0] + 1)*sizeof(GRD_data_type));
		}
	return (out.good()? 0: -4);
}

MC33_real surface::get_isovalue() { return iso; }

unsigned int surface::get_num_vertices() { return nV; }

unsigned int surface::get_num_triangles() { return nT; }

surface::surface() : nV(0), nT(0), sflag(8) {}

void surface::clear() {
	T.clear();
	V.clear();
	N.clear();
	color.clear();
	nV = nT = 0;
	sflag = 8;
}

void surface::adjustvectorlenght() {
	T.resize(nT);
	T.shrink_to_fit();
	V.resize(nV);
	V.shrink_to_fit();
	N.resize(nV);
	N.shrink_to_fit();
	color.shrink_to_fit();
	sflag &= ~8;
}


const unsigned int *surface::getTriangle(unsigned int n) { return (n < nT? T[n].v: 0); }

const MC33_real *surface::getVertex(unsigned int n) { return (n < nV? V[n].v: 0); }

const float *surface::getNormal(unsigned int n) { return (n < nV? N[n].v: 0); }

void surface::flipNormals() {
	unsigned int n = 3*nV;
	float *v = N[0].v;
	for (unsigned int i = 0; i != n; i++)
		v[i] = -v[i];
	sflag ^= 4;
}

void surface::flipTriangles() {
	for (auto &t : T)
		std::swap(t.v[0], t.v[1]);
	sflag ^= 2;
}

void surface::setColor(unsigned int n, unsigned char *pcolor) {
	if (n < nV)
		color[n] = *(reinterpret_cast<int*>(pcolor));
}

const unsigned char* surface::getColor(unsigned int n) {
	return (n < nV? reinterpret_cast<unsigned char*>(&color[n]): 0);
}


#if MC33_DOUBLE_PRECISION
#define MC33_surf_magic_num 0x6575732e //".sud"
#define MC33_surf_magic_nu2 0x7075732e
#define MC33_real2 float
#else
#define MC33_surf_magic_num 0x7075732e //".sup"
#define MC33_surf_magic_nu2 0x6575732e
#define MC33_real2 double
#endif

int surface::save_bin(const char *filename) {
	if (sflag&8)
		adjustvectorlenght();
	std::ofstream out(filename, std::ios::binary);
	if (!out)
		return -1;

	int n = MC33_surf_magic_num;
	out.write(reinterpret_cast<char*>(&n),sizeof(int));
	out.write(reinterpret_cast<char*>(&iso),sizeof(MC33_real));
	out.write(reinterpret_cast<char*>(&nV),sizeof(int));
	out.write(reinterpret_cast<char*>(&nT),sizeof(int));

	out.write(reinterpret_cast<char*>(&T[0]),3*sizeof(int)*nT);
	out.write(reinterpret_cast<char*>(&V[0]),3*sizeof(MC33_real)*nV);
	out.write(reinterpret_cast<char*>(&N[0]),3*sizeof(float)*nV);
	out.write(reinterpret_cast<char*>(&color[0]),sizeof(int)*nV);
	return (out.good()? 0: -1);
}

int surface::read_bin(const char *filename) {
	std::ifstream in(filename, std::ios::binary);
	if (!in)
		return -1;
	int n;
	in.read(reinterpret_cast<char*>(&n),sizeof(int));
	if (n == MC33_surf_magic_num || n == MC33_surf_magic_nu2)
	{
		if (n == MC33_surf_magic_nu2) {
			MC33_real2 t;
			in.read(reinterpret_cast<char*>(&t),sizeof(MC33_real2));
			iso = t;
		} else
			in.read(reinterpret_cast<char*>(&iso),sizeof(MC33_real));
		in.read(reinterpret_cast<char*>(&nV),sizeof(int));
		in.read(reinterpret_cast<char*>(&nT),sizeof(int));
		adjustvectorlenght();
		in.read(reinterpret_cast<char*>(&T[0]),3*sizeof(int)*nT);
		if (n == MC33_surf_magic_nu2) {
			for (unsigned int j = 0; j != nV; j++) {
				MC33_real2 t[3];
				in.read(reinterpret_cast<char*>(t),3*sizeof(MC33_real2));
				for (int i = 0; i != 3; i++)
					V[j].v[i] = t[i];
			}
		} else
			in.read(reinterpret_cast<char*>(V[0].v),3*sizeof(MC33_real)*nV);
		in.read(reinterpret_cast<char*>(&N[0]),3*sizeof(int)*nV);
		in.read(reinterpret_cast<char*>(&color[0]),sizeof(int)*nV);
	} else
		return -1;
	return (in.good()? 0: -1);
}
#undef MC33_surf_magic_num
#undef MC33_surf_magic_nu2
#undef MC33_real2

#if MC33_DOUBLE_PRECISION
#define MC33_prec 11
#else
#define MC33_prec 6
#endif
int surface::save_txt(const char* filename) {
	std::ofstream out(filename);
	if (sflag&8)
		adjustvectorlenght();
	if (!out)
		return -1;

	out << "isovalue: ";
	out.setf(std::ios_base::scientific, std::ios_base::floatfield);
	out.precision(MC33_prec);
	out << iso << "\n\nVERTICES:\n" << nV << "\n\n";
	out.setf(std::ios_base::fixed, std::ios_base::basefield);
	out.precision(MC33_prec);
	for (const auto &r: V)
		out << std::setw(7 + MC33_prec) << r.v[0] << std::setw(8 + MC33_prec) << r.v[1] <<
			std::setw(8 + MC33_prec) << r.v[2] << std::endl;

	out << "\n\nTRIANGLES:\n" << nT << "\n\n";
	for (const auto &t: T)
		out << std::setw(8) << t.v[0] << " " << std::setw(8) << t.v[1] << " " << std::setw(8) << t.v[2] << std::endl;

	out << "\n\nNORMALS:\n";
	out.precision(5);
	for (const auto &r: N)
		out << std::setw(12) << r.v[0] << std::setw(13) << r.v[1] << std::setw(13) << r.v[2] << std::endl;

	out << "\n\nCOLORS:\n";
	for (const auto &c: color)
		out << c << std::endl;
	out << "\nEND\n";
	return (out.good()? 0: -1);
}

int surface::save_obj(const char* filename) {
	std::ofstream out(filename);
	if (sflag&8)
		adjustvectorlenght();
	if (!out)
		return -1;

	out << "# isovalue: ";
	out.precision(MC33_prec);
	out << iso << "\n# VERTICES " << nV << ":\n";
	out.precision(MC33_prec);
	for (const auto &r: V)
		out << "v " << r.v[0] << " " << r.v[1] << " " << r.v[2] << std::endl;

	out << "# NORMALS:\n";
	out.precision(5);
	for (const auto &r: N)
		out << "vn " << r.v[0] << " " << r.v[1] << " " << r.v[2] << std::endl;

	out << "# TRIANGLES " << nT << ":\n";
	for (const auto &t: T) {
		std::string s0 = std::to_string(t.v[0] + 1);
		std::string s1 = std::to_string(t.v[1] + 1);
		std::string s2 = std::to_string(t.v[2] + 1);
		out << "f " << s0 << "//" << s0 << " " << s1 << "//" << s1 << " " << s2 << "//" << s2 << std::endl;
	}
	out << "# END";
	return (out.good()? 0: -1);
}

int surface::save_ply(const char* filename, const char* author, const char* object) {
	std::ofstream out(filename);
	if (sflag&8)
		adjustvectorlenght();
	if (!out)
		return -1;
	char empty = 0;
	if (!author)
		author = &empty;
	if (!object)
		object = &empty;
	out << "ply\nformat ascii 1.0\ncomment author: " << author << "\ncomment object: " << object;
	out << "\nelement vertex " << nV << "\nproperty float x\nproperty float y\nproperty float z";
	out << "\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement face " << nT;
	out << "\nproperty list uchar int vertex_index\nend_header";
	for (unsigned int i = 0; i < nV; i++) {
		auto v = V[i].v;
		out << "\n" << v[0] << " " << v[1] << " " << v[2] << " ";
		unsigned char* c = reinterpret_cast<unsigned char*>(&color[i]);
		out << (int)c[0] << " " << (int)c[1] << " " << (int)c[2];
	}
	for (unsigned int i = 0; i < nT; i++) {
		auto v = T[i].v;
		out << "\n3 " << v[0] << " " << v[1] << " " << v[2];
	}
	out << "\n";
	return (out.good()? 0: -1);
}

#undef MC33_prec

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#endif // mc33cpp_implementation

#endif // MC33cpp_h_
