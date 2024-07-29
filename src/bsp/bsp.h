#pragma once
#include <cstdint>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "mesh.h"
#include "util.h"

#include <functional>


#include <slimcpplib/long_int.h>




typedef Eigen::Matrix<int32_t, 2, 1>                            Vector2i;
typedef Eigen::Matrix<int32_t, 3, 1>                            Vector3i;
typedef Eigen::Matrix<double, 3, 1>                             vector3d;


using namespace slim::literals;
using uint128_t = slim::uint128_t;
using int128_t = slim::int128_t;
using uint256_t = slim::uint256_t;
using int256_t = slim::int256_t;


// https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

// we can use the double to validate our exatc arithmetic
template<typename T>
int min_axis(T& d)
{
	int minj = 0;
	double minval = std::abs(d[0]);
	for (int j = 1; j < 3; j++)
	{
		double absval = std::abs(d[j]);
		if (absval < minval) {
			minval = absval;
			minj = j;
		}
	}

	return minj;
}

template<typename T>
int max_axis(T& d)
{
	int maxj = 0;
	double maxval = std::abs(d[0]);
	for (int j = 1; j < 3; j++)
	{
		double absval = std::abs(d[j]);
		if (absval > maxval) {
			maxval = absval;
			maxj = j;
		}
	}

	return maxj;
}


extern int cubeface[6][4];
extern int childoctantcenteroffset[8][3];

// 128 bit
struct BspPlane
{
	int32_t a, b, c;
	int64_t d;
};


struct RawQuad
{
	Eigen::Vector3d p0, p1, p2, p3;
};

inline void buildframe(Eigen::Vector3d& n, Eigen::Vector3d& x, Eigen::Vector3d& y)
{
	Eigen::Vector3d refd;
	refd << 0, 0, 0;
	refd[min_axis(n)] = 1.0;

	x = refd.cross(n);
	y = n.cross(x);
	x.normalize();
	y.normalize();
}



inline RawQuad plane2quad(const BspPlane& plane, double r)
{
	Eigen::Vector3d n;
	n << plane.a, plane.b, plane.c;

	const double len = n.norm();
	n.normalize();
	Eigen::Vector3d axisx, axisy;
	buildframe(n, axisx, axisy);

	const double d = -plane.d / len;
	Eigen::Vector3d p0 = n * d;
	p0 -= axisx * r / 2.0;
	p0 -= axisy * r / 2.0;

	Eigen::Vector3d p1 = p0 + axisx * r;
	Eigen::Vector3d p2 = p1 + axisy * r;
	Eigen::Vector3d p3 = p0 + axisy * r;

	return { p0, p1, p2, p3 };
}

inline RawQuad plane2quad(const BspPlane& plane, const Eigen::Vector3d& p, double r)
{
	Eigen::Vector3d n;
	n << plane.a, plane.b, plane.c;

	const double len = n.norm();
	n.normalize();
	Eigen::Vector3d axisx, axisy;
	buildframe(n, axisx, axisy);

	const double d = -plane.d / len;
	Eigen::Vector3d p0 = n * d;
	p0 -= axisx * r / 2.0;
	p0 -= axisy * r / 2.0;

	auto dir = p - p0;
	auto x = axisx.dot(dir);
	auto y = axisy.dot(dir);

	p0 += x * axisx + y * axisy;

	Eigen::Vector3d p1 = p0 + axisx * r;
	Eigen::Vector3d p2 = p1 + axisy * r;
	Eigen::Vector3d p3 = p0 + axisy * r;

	return { p0, p1, p2, p3 };
}

template<typename T>
pmp::Vertex addvertex(Mesh& sh, T* p) {
	pmp::Point point;
	point[0] = p[0];
	point[1] = p[1];
	point[2] = p[2];
	return sh.msh.add_vertex(point);
};

template<typename Mesh, typename T>
pmp::Vertex addvertex0(Mesh& sh, T* p) {
	pmp::Point point;
	point[0] = p[0];
	point[1] = p[1];
	point[2] = p[2];
	return sh.add_vertex(point);
};


void dumpplane(const BspPlane& pl, const Eigen::Vector3d& refpt, const double r);



inline void addrawquad(Mesh& sh, const BspPlane& plane, std::vector<pmp::Vertex>& vertices, double r=5) {
	auto rq = plane2quad(plane, r);

	vertices.clear();
	vertices.push_back(addvertex(sh, rq.p0.data()));
	vertices.push_back(addvertex(sh, rq.p1.data()));
	vertices.push_back(addvertex(sh, rq.p2.data()));
	vertices.push_back(addvertex(sh, rq.p3.data()));

	sh.msh.add_face(vertices);
};


struct Homogeneous4DPoint
{
	int128_t x1, x2, x3, x4;
};

inline Homogeneous4DPoint tohd4(int *p)
{
	Homogeneous4DPoint hd;
	hd.x1 = p[0];
	hd.x2 = p[1];
	hd.x3 = p[2];
	hd.x4 = 1;

	return hd;
}

enum PolyCutSplit
{
	PC_LEFT,
	PC_RIGHT,
	PC_CUT,
	PC_ON,
	PC_NON,
};

struct TripleTriple
{
	int p0[3];
	int p1[3];
	int p2[3];
};

struct PolyCut
{
	// 
	std::vector<Homogeneous4DPoint> pts;
	std::vector<vector3d> p3d;
	std::vector<int> hnxt;
	std::vector<int> hpre;
	std::vector<int> htov;
	std::vector< BspPlane > planes;
	std::vector< std::pair<int, int> > hplanes;
	vector3d center;

	// reuse during split
	std::vector<int> sh;
	std::vector< std::pair<int, int> > hssign;
	std::vector<int> edges;

	// journal
	bool isjournal = false;
	std::string fn;
	FILE* fp = nullptr;
	~PolyCut() {
		closejournal();
	}

	void init(int* p0, int* p1, int* p2);
	int split(const BspPlane& pl, int starth, int newhs[2]);

	void closejournal() const{
		if (nullptr != fp) {
			fclose(fp);
		}
	}

	void dumpoly(const std::string &fn) const;
};



struct CutEntity
{
	int entity;
	int type;
};


struct Node
{
	// uint32_t f = -1;
	int left = -1;
	int right = -1;
	int plane = -1;
	bool isleftout = true;// normal direction side is left
	bool isleft = false;
	bool isright = false;
};


struct EdgePlane
{
	int p0, p1;
};


struct Bsp
{
	PolyCut pc;

	std::vector<Node> nodes;
	std::vector<BspPlane> planes;

	pmp::FaceProperty<int> supportplaneprop;
	pmp::VertexProperty<Homogeneous4DPoint> hdpts;
	pmp::EdgeProperty< EdgePlane> edgeplanepair;
	pmp::EdgeProperty< bool > edgeiscut;
	pmp::VertexProperty<char> vertexflag;
	pmp::EdgeProperty<bool> edgeflag;
	pmp::FaceProperty<bool> faceflag;
};


void plane_from_points(int* p0, int* p1, int* p3, BspPlane& plane);
int planelinesegmentintersection(int* p0, int* p1, const BspPlane& pl, Homogeneous4DPoint& hd4);
int virtualtriangleplanes(int* p0, int* p1, int* p2,
	BspPlane& pl0, BspPlane& pl1, BspPlane& pl2);


void intersect_3_planes(
	const BspPlane& p0,
	const BspPlane& p1, 
	const BspPlane& p2, 
	Homogeneous4DPoint &p);

int classify_vertex(const BspPlane &s, const Homogeneous4DPoint& p);
int classify_vertex(const BspPlane& s, int* p);


// arithmatic
int64_t determinant(int64_t* c0, int64_t* c1, int64_t* c2);

int128_t determinant128(int64_t* c0, int64_t* c1, int64_t* c2);

void build_cube(Bsp &bsp, Mesh& sh, const double l);

int convexmesh_cut(Bsp &bsp, Mesh& sh, BspPlane &pl);


int mesh2bsp(Mesh& sh, bool cutter=false, bool iscuttercheckpoint=false, const char* fn=nullptr);

struct Octant
{
	float x, y, z;// center
	float extent;// half size

	uint32_t childid = std::numeric_limits<uint32_t>::max();
	std::vector<uint32_t> childtris;
};


typedef std::function<void(const int triid, float* p0, float* p1, float* p2)> TriFun;

constexpr int Octant_MAX_TRI = 20;
struct Octree
{
	std::vector<Octant> octantlist;

	void build(uint32_t nodeid, std::vector<uint32_t> &trilist, TriFun fun);
};





// debug
void printedgeimpl(pmp::SurfaceMesh& msh, pmp::Halfedge h);
std::vector<pmp::Face> collectfacefromseed(pmp::SurfaceMesh& msh, pmp::Face f);
std::vector<pmp::Face> collectfacefromseed(Mesh& sh, pmp::Face f);
std::vector<pmp::Face> collectfacefromseedv(Mesh& sh, pmp::Vertex v3);
void dump(std::string fn, pmp::SurfaceMesh& msh, std::vector<pmp::Face>& faces, bool isdump=false);
void dumpmesh(Mesh& sh, const std::string& fn__, const int i, bool isdump);
void dumpmesh(const int i, Mesh& sh);

