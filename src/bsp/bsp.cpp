#include "bsp.h"

#include <assert.h>
#include <array>

#include "pmp/bounding_box.h"
#include "pmp/io/write_obj.h"
#include "pmp/algorithms/normals.h"
#include "bvh.h"


int triBoxOverlap(float boxcenter[3],float boxhalfsize[3],float triverts[3][3]);

int cubeface[6][4] = {
	{1, 0, 3, 2},// btm
	{4, 5, 6, 7},// top
	{3, 0, 4, 7},// left
	{1, 2, 6, 5},// right
	{0, 1, 5, 4},// front 
	{2, 3, 7, 6},// back

};


Vector3i cross(Vector3i& a, Vector3i& b)
{
	Vector3i n;
	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = a[2] * b[0] - a[0] * b[2];
	n[2] = a[0] * b[1] - a[1] * b[0];

	return n;
}

int64_t determinant(int64_t* c0, int64_t* c1, int64_t* c2)
{
	auto d0 = c0[0] * (c1[1] * c2[2] - c2[1] * c1[2]);
	auto d1 = c1[0] * (c0[1] * c2[2] - c2[1] * c0[2]);
	auto d2 = c2[0] * (c0[1] * c1[2] - c1[1] * c0[2]);

#if 0
	Eigen::Matrix3d m;
	m << c0[0], c1[0], c2[0],
		c0[1], c1[1], c2[1],
		c0[2], c1[2], c2[2];

	const double det = m.determinant();
#endif

	const int64_t det0 =  d0 - d1 + d2;

	return det0;
}

int128_t determinant128(int64_t* c0_, int64_t* c1_, int64_t* c2_)
{
	int128_t c0[3] = { c0_[0], c0_[1], c0_[2] };
	int128_t c1[3] = { c1_[0], c1_[1], c1_[2] };
	int128_t c2[3] = { c2_[0], c2_[1], c2_[2] };

	int128_t d0 = c0[0] * (c1[1] * c2[2] - c2[1] * c1[2]);
	int128_t d1 = c1[0] * (c0[1] * c2[2] - c2[1] * c0[2]);
	int128_t d2 = c2[0] * (c0[1] * c1[2] - c1[1] * c0[2]);

#if 0
	Eigen::Matrix3d m;
	m << c0[0], c1[0], c2[0],
		c0[1], c1[1], c2[1],
		c0[2], c1[2], c2[2];

	const double det = m.determinant();
#endif

	const int128_t det0 = d0 - d1 + d2;

	return det0;
}

void plane_from_points(int* p0, int* p1, int* p2, BspPlane& plane)
{
	Vector3i v0(p0[0], p0[1], p0[2]);
	Vector3i v1(p1[0], p1[1], p1[2]);
	Vector3i v2(p2[0], p2[1], p2[2]);

	Vector3i d01 = v1 - v0;
	Vector3i d02 = v2 - v0;

	Vector3i n = cross(d01, d02);

	int d = v0.dot(n);
	int d1 = v1.dot(n);
	int d2 = v2.dot(n);
	assert(d == v1.dot(n) && d == v2.dot(n));

	plane.a = n[0];
	plane.b = n[1];
	plane.c = n[2];
	plane.d = -d;
}

void intersect_3_planes(
	const BspPlane& p0,
	const BspPlane& p1,
	const BspPlane& p2,
	Homogeneous4DPoint& p)
{
	int64_t ca[3] = { p0.a, p1.a, p2.a };
	int64_t cb[3] = { p0.b, p1.b, p2.b };
	int64_t cc[3] = { p0.c, p1.c, p2.c };
	int64_t cd[3] = { p0.d, p1.d, p2.d };

	p.x1 = determinant128(cd, cb, cc);
	p.x2 = determinant128(ca, cd, cc);
	p.x3 = determinant128(ca, cb, cd);
	p.x4 = -determinant128(ca, cb, cc);
}





struct CubePlane
{
	BspPlane btmplane ;
	BspPlane topplane ;
	BspPlane leftplane ;
	BspPlane rightplane ;
	BspPlane frontplane ;
	BspPlane backplane ;

	void initcubeplane(const double l)
	{
		int il = int(l);
		btmplane = { 0, 0, -1, 0 };
		topplane = { 0, 0, 1, -il };
		leftplane = { 0, -1, 0, 0 };
		rightplane = { 0, 1, 0, -il };
		frontplane = { 1, 0, 0, -il };
		backplane = { -1, 0, 0, 0 };
	}
};

struct CubePlaneIds
{
	int btm;
	int top ;
	int left ;
	int right ;
	int front ;
	int back ;

	void init(const int offset)
	{
		btm = offset;
		top = btm + 1;
		left = top + 1;
		right = left + 1;
		front = right + 1;
		back = front + 1;
	}

	std::vector< std::array<int, 3> > faceindex()
	{
		std::vector< std::array<int, 3> > cubept = {
			{front, btm, left},
			{right, btm, front},
			{right, btm, back},
			{back, btm, left},

			{front, top, left},
			{right, top, front},
			{right, top, back},
			{back, top, left},
		};

		return cubept;
	}

};


void buildcubebsp(Bsp& bsp, CubePlaneIds& cpids)
{
	bsp.nodes.resize(6);

	bsp.nodes[0].right = 1;
	bsp.nodes[0].plane = cpids.left;

	bsp.nodes[1].right = 2;
	bsp.nodes[1].plane = cpids.front;

	bsp.nodes[2].right = 3;
	bsp.nodes[2].plane = cpids.right;

	bsp.nodes[3].right = 4;
	bsp.nodes[3].plane = cpids.top;

	bsp.nodes[4].right = 5;
	bsp.nodes[4].plane = cpids.back;

	bsp.nodes[5].plane = cpids.btm;

	bsp.nodes[5].right = 3;
	bsp.nodes[5].isright = true;
}

void buildcubemesh(Bsp& bsp, Mesh &sh, std::vector< std::array<int, 3> > &cubept)
{

	bsp.supportplaneprop = sh.msh.add_face_property<int>("facesupportplane");
	bsp.hdpts = sh.msh.add_vertex_property< Homogeneous4DPoint>("vertexhd");
	bsp.edgeplanepair = sh.msh.add_edge_property<EdgePlane>("edgeplanepair");
	bsp.edgeiscut = sh.msh.add_edge_property<bool>("edgeiscut");
	bsp.vertexflag = sh.msh.add_vertex_property<char>("vertexflag");
	bsp.edgeflag = sh.msh.add_edge_property<bool>("edgeflag");
	bsp.faceflag = sh.msh.add_face_property<bool>("faceflag");



	std::vector<pmp::Vertex> vertices;
	for (int i = 0; i < 8; i++)
	{
		auto& pl0 = bsp.planes[cubept[i][0]];
		auto& pl1 = bsp.planes[cubept[i][1]];
		auto& pl2 = bsp.planes[cubept[i][2]];

		auto rwpl0 = bspplane2plane(pl0);
		auto rwpl1 = bspplane2plane(pl1);
		auto rwpl2 = bspplane2plane(pl2);

		auto ipt0 = IntersectPlanesLU(rwpl0, rwpl1, rwpl2);
		pmp::Point point(ipt0[0], ipt0[1], ipt0[2]);
		auto v = sh.msh.add_vertex(point);
		vertices.push_back(v);

		intersect_3_planes(pl0, pl1, pl2, bsp.hdpts[v]);


#if 0
		auto hd4 = bsp.hdpts[v];
		double x4 = double(hd4.x4);
		Eigen::Vector3d pt;
		pt << hd4.x1 / x4, hd4.x2 / x4, hd4.x3 / x4;

		double d0 = ipt0[0] - pt[0];
		double d1 = ipt0[1] - pt[1];
		double d2 = ipt0[2] - pt[2];
		d0 = std::abs(d0);
		d1 = std::abs(d1);
		d2 = std::abs(d2);

		assert(d0 < PLANEEPSILON);
		assert(d1 < PLANEEPSILON);
		assert(d2 < PLANEEPSILON);
#endif
	}

	std::vector<pmp::Vertex> facevertices;
	for (int i = 0; i < 6; i++)
	{
		facevertices.clear();
		for (int j = 0; j < 4; j++)
		{
			facevertices.push_back(vertices[cubeface[i][j]]);
		}

		auto newface = sh.msh.add_face(facevertices);
		bsp.supportplaneprop[newface] = i;
	}

	for (auto e : sh.msh.edges())
	{
		auto h0 = sh.msh.halfedge(e, 0);
		auto h1 = sh.msh.opposite_halfedge(h0);
		auto f0 = sh.msh.face(h0);
		auto f1 = sh.msh.face(h1);
		bsp.edgeplanepair[e].p0 = bsp.supportplaneprop[f0];
		bsp.edgeplanepair[e].p1 = bsp.supportplaneprop[f1];
	}

	if (0)
		sh.write("cube.obj");
}


void build_cube(Bsp& bsp, Mesh& sh, const double l)
{
	CubePlane cp;
	cp.initcubeplane(l);

	const int offset = int(bsp.planes.size());
	bsp.planes.push_back(cp.btmplane);
	bsp.planes.push_back(cp.topplane);
	bsp.planes.push_back(cp.leftplane);
	bsp.planes.push_back(cp.rightplane);
	bsp.planes.push_back(cp.frontplane);
	bsp.planes.push_back(cp.backplane);

	CubePlaneIds cpids;
	cpids.init(offset);
	auto cubept = cpids.faceindex();

	buildcubemesh(bsp, sh, cubept);

	buildcubebsp(bsp, cpids);

	
	// const int nodeoffset = int(bsp.nodes.size());
	// bsp.nodes.emplace_back();
	// bsp.nodes[nodeoffset].left = 0;
}


inline double planedist(const pmp::Point& p_, Plane& rawpl)
{
	Eigen::Vector3d p;
	p << p_[0], p_[1], p_[2];

	return rawpl.Normal.dot(p) - rawpl.Distance;
}


union CutEntityUnion
{
	pmp::Halfedge h;
	pmp::Vertex v;
};

struct CutEntity0
{
	CutEntityUnion entity;
	int type;
};


std::vector< CutEntity> convexmeshcut_traversal(Bsp& bsp, Mesh& sh, BspPlane& pl, pmp::Vertex v, pmp::Halfedge nxth, int type);

std::vector< CutEntity> convexmeshcut_traversal_(Bsp& bsp, Mesh& sh, BspPlane& pl, pmp::Vertex v, pmp::Halfedge nxth, int type)
{
	std::vector< CutEntity> entities;

	
	if (type == 0)
	{
		auto range = sh.msh.halfedges(v);
		int pretype = -1;
		bool initpretype = false;
		for (auto h : range) {
			auto v1 = sh.msh.to_vertex(h);

			const int type1 = classify_vertex(pl, bsp.hdpts[v1]);

			if (!initpretype) {
				initpretype = true;
				pretype = type1;
			}
			else {
				if (type1 != pretype&&0!=type1&&0!=pretype)
				{
					nxth = h;
					break;
				}
			}
		}

		if (!nxth.is_valid())
			return entities;// convext cut this plane coninside exiting face

		auto h0 = nxth;
		auto opnxth = sh.msh.opposite_halfedge(nxth);
		nxth = sh.msh.next_halfedge(opnxth);// this face will be cut
		// nxth = opnxth;

		// assert(sh.msh.to_vertex(h0) == sh.msh.to_vertex(nxth));

		CutEntity ent = { int(opnxth.idx()), 0 };
		entities.push_back(ent);
	}
	else {
		assert(nxth.is_valid());
		auto oph = sh.msh.opposite_halfedge(nxth);
		auto f0 = sh.msh.face(nxth);
		auto f1 = sh.msh.face(oph);

		CutEntity ent = { int(nxth.idx()), 1};
		entities.push_back(ent);
	}


	auto poppointcut = [&](pmp::Halfedge h1, pmp::Halfedge crrh) {
		auto opnxtcrrh = sh.msh.prev_halfedge(h1);
		// h1 = sh.msh.next_halfedge(h1);
		CutEntity ent = { int(opnxtcrrh.idx()), 0 };

		auto preh = pmp::Halfedge(entities.back().entity);
		assert(sh.msh.face(preh) == sh.msh.face(crrh));
		entities.push_back(ent);

		auto h1f = sh.msh.face(h1);
		auto opnxtcrrhf = sh.msh.face(opnxtcrrh);
		assert(h1f == opnxtcrrhf);
	};

	auto crrh = nxth;
	bool initface = false;
	int pretype = -1;
	while (true)
	{
		auto h1 = sh.msh.next_halfedge(crrh);
		auto v = sh.msh.to_vertex(crrh);
		int type1 = classify_vertex(pl, bsp.hdpts[v]);
		
		if (!initface) {
			pretype = type1;
			initface = true;
		}
		else {
			if (pretype == 0) {
				if (0 == type1) {// similar pre post status update as the normal case
					h1 = sh.msh.opposite_halfedge(crrh);
					
					poppointcut(h1, crrh);

					initface = false;
					pretype = -1;
				}
				else {
					pretype = type1;// update 
				}
			}
			else if (pretype != type1) {
				// jump to neighbor face
				h1 = sh.msh.opposite_halfedge(crrh);

				if (0 == type1) {// point cut
					poppointcut(h1, crrh);
				}
				else if (pretype != type1) {// edge cut
					CutEntity ent = { int(h1.idx()), 1 };
					entities.push_back(ent);
				}

				initface = false;
				pretype = -1;
			}
		}

		crrh = h1;
		pretype = type1;

		if (crrh == nxth) {
			break;// reach to the edge
		}
	}

	return entities;
}


struct InitCfg
{
	pmp::Halfedge nxth;
	pmp::Vertex v;
	int type = -3;
};

struct ScopedCleanVertexFalg
{
	std::vector<int>& taggedvertex;
	Bsp& bsp;
	~ScopedCleanVertexFalg()
	{
		for (auto v : taggedvertex)
		{
			bsp.vertexflag[pmp::Vertex(v)] = 0;
		}
	}
};
InitCfg vertex_traversal(Bsp& bsp, Mesh& sh, BspPlane& pl, Plane &rawpl,
	pmp::Vertex v, const int type)
{
	if (0)
	{
		int types[3] = { 0, 0, 0 };
		const size_t nv = sh.msh.n_vertices();
		for (int i = 0; i < nv; i++)
		{
			int type = classify_vertex(pl, bsp.hdpts[pmp::Vertex(i)]);
			types[type + 1]++;
		}

		printf("types, %d %d %d\n", types[0], types[1], types[2]);
	}

	std::vector<int> taggedvertex;
	ScopedCleanVertexFalg scopped = { taggedvertex, bsp };

	InitCfg rtn;

	while (true)
	{
		if (0 == type) {
			assert(false);
			break;// found point cutopposite_halfedge
		}

		bsp.vertexflag[v] = 1;
		taggedvertex.push_back(v.idx());

		auto range = sh.msh.halfedges(v);
		pmp::Vertex nxtv;
		double mind = std::numeric_limits<double>::max();
		for (auto h : range) {
			auto tov = sh.msh.to_vertex(h);

			if (1 == bsp.vertexflag[tov])
				continue;

			auto p0 = sh.msh.position(tov);
			int type0 = classify_vertex(pl, bsp.hdpts[tov]);
			if (0 == type0) {
				rtn.v = tov;
				rtn.type = 0;
				return rtn;
				// break;// found point cut;
			}
			if (type != type0) {
				rtn.nxth = h;
				rtn.v = tov;	
				rtn.type = 1;
				return rtn;
				// break;// found edge cut
			}

			const double d0 = planedist(p0, rawpl)*type0;
			assert(d0 > PLANEEPSILON);

			if (mind > d0) {
				mind = d0;
				nxtv = tov;
			}
		}

		if (!nxtv.is_valid()) {
			break;// corner case for linear search
		}

		assert(nxtv != v);
		v = nxtv;
	}

	return rtn;
}

struct ScopedCleanEdgeFalg
{
	std::vector<pmp::Edge> &edges;
	Bsp& bsp;
	~ScopedCleanEdgeFalg()
	{
		for (auto val : edges)
		{
			bsp.edgeiscut[val] = false;
		}
	}
};


static int id = 0;
std::pair<int ,int> convexmesh_cut(Bsp& bsp, Mesh& sh, const int newplid, pmp::Vertex v3= pmp::Vertex(3))
{
	auto testfun = [&]() {
		auto faces = collectfacefromseedv(sh, v3);

		std::string fn = "";
		fn += std::to_string(newplid);
		fn += ".obj";
		dump(fn, sh.msh, faces);

		return faces.size();
	};

	id++;

	auto sz0 = testfun();
	if (0)
	{



		{
			Mesh plsh;
			std::vector<pmp::Vertex> vertices;
			addrawquad(plsh, bsp.planes[newplid], vertices, 10000);
			std::string fn = std::to_string(newplid);
			fn += "plane.obj";
			plsh.write(fn);
		}



	}


	std::pair<int, int> rntcode;
	BspPlane& pl = bsp.planes[newplid];

	auto rawpl = bspplane2plane(pl);

	// auto v3 = pmp::Vertex(3);
	auto p3 = sh.msh.position(v3);
	const double d0 = planedist(p3, rawpl);

	const int type = classify_vertex(pl, bsp.hdpts[v3]);

	InitCfg initcfg;
	if (0 == type) {
		initcfg.v = v3;
		initcfg.type = 0;
	}
	else {
		initcfg = vertex_traversal(bsp, sh, pl, rawpl, v3, type);
		if (0 == initcfg.type) {
			initcfg.nxth = pmp::Halfedge();
		}
		else if (1 == initcfg.type)
		{
			initcfg.v = pmp::Vertex();// make it invalid for edge cut
			assert(1 == initcfg.type);
		}
	}

	if (-3 == initcfg.type)
		return { -1, 0 };//


	std::vector<pmp::Edge> cutedges;
	std::vector<pmp::Halfedge> cuthalfedges;
	ScopedCleanEdgeFalg scoped = { cutedges, bsp };
	auto popcutedge = [&](pmp::Halfedge h) {
		cuthalfedges.push_back(h);
		auto e = sh.msh.edge(h);
		cutedges.push_back(e);
		bsp.edgeiscut[e] = true;
	};

	auto entities = convexmeshcut_traversal(bsp, sh, pl,
		initcfg.v, initcfg.nxth, initcfg.type);

	int total = 0;
	for (auto ent : entities)
		total += ent.entity;



	if (entities.empty())
		return { -2, 0 };


	for (int i = 0; i < entities.size() - 1; i++)
	{
		if (0 == entities[i].type)
			continue;// point cut

		assert(entities[i].type == 1);
		auto h = pmp::Halfedge(entities[i].entity);
		auto oph = sh.msh.opposite_halfedge(h);
		auto v0 = sh.msh.to_vertex(oph);
		auto v1 = sh.msh.to_vertex(h);

		auto newv = sh.msh.new_vertex();
		auto newh = sh.msh.insert_vertex(h, newv);

		assert(sh.msh.to_vertex(h) == newv);

		auto f0 = sh.msh.face(h);
		auto f1 = sh.msh.face(sh.msh.opposite_halfedge(h));
		auto e = sh.msh.edge(h);
		int id0 = bsp.edgeplanepair[e].p0;
		int id1 = bsp.edgeplanepair[e].p1;
		auto newe = sh.msh.edge(newh);
		bsp.edgeplanepair[newe] = bsp.edgeplanepair[e];

		auto pl0 = bsp.planes[id0];// bsp.supportplaneprop[f0]];
		auto pl1 = bsp.planes[id1];// bsp.supportplaneprop[f1]];

		auto ipt = IntersectPlanesLU(pl0, pl1, pl);

		Homogeneous4DPoint hd4;
		intersect_3_planes(pl0, pl1, pl, hd4);
		bsp.hdpts[newv] = hd4;

		auto& p = sh.msh.position(newv);
		p[0] = ipt[0];
		p[1] = ipt[1];
		p[2] = ipt[2];
	}

	auto gettailh = [&](pmp::Halfedge h, pmp::Vertex v) {

		while (true)
		{
			auto tov = sh.msh.to_vertex(h);
			if (v == tov)
				return h;

			h = sh.msh.next_halfedge(h);
		}

		return pmp::Halfedge();
	};

	auto getheadh = [&](pmp::Vertex v, pmp::Halfedge h) {
		auto oph = sh.msh.opposite_halfedge(h);
		auto preh = sh.msh.prev_halfedge(oph);

		return gettailh(oph, v);
	};

	auto issharededge = [&](pmp::Vertex v0, pmp::Vertex v1) {
		auto range = sh.msh.halfedges(v0);
		for (auto h : range) {
			auto tov = sh.msh.to_vertex(h);
			if (v1 == tov) {
				return h;
			}
		}

		return pmp::Halfedge();
	};

	auto vv = [&](pmp::Vertex v0, pmp::Vertex v1) {

		auto r0 = sh.msh.halfedges(v0);
		auto r1 = sh.msh.halfedges(v1);

		for (auto h0 : r0)
		{
			auto f0 = sh.msh.face(h0);
			for (auto h1 : r1)
			{
				auto f1 = sh.msh.face(h1);
				if (f0 == f1)
				{
					return std::make_pair(sh.msh.prev_halfedge(h0), sh.msh.prev_halfedge(h1));
				}
			}
		}

		return std::make_pair(pmp::Halfedge(), pmp::Halfedge());
	};

	// special case the head and tail split
	// tail h will change
	pmp::Halfedge headpreoph1;
	if (0 == entities[0].type) {
		auto h1 = pmp::Halfedge(entities[0].entity);
		auto h1nxt = sh.msh.next_halfedge(h1);
		auto oph1nxt = sh.msh.opposite_halfedge(h1nxt);// inverse op of type==0 entit pop
		headpreoph1 = oph1nxt;
	}
	for (int i = 0; i < entities.size() - 1; i++)
	{
		// assert(entities[i].type == 1);
		// assert(entities[i + 1].type == 1);
		if (0 == entities[i].type && 0 == entities[i + 1].type) {
			auto sv0 = pmp::Vertex(entities[i].entity);
			auto sv1 = pmp::Vertex(entities[i + 1].entity);
			auto sharedh = issharededge(sv0, sv1);
			if (sharedh.is_valid()) {
				popcutedge(sharedh);
				continue;
			}
		}
		auto h0 = pmp::Halfedge(entities[i].entity);
		auto h1 = pmp::Halfedge(entities[i + 1].entity);
		auto oph1 = sh.msh.opposite_halfedge(h1);
		auto preoph1 = sh.msh.prev_halfedge(oph1);

		if (0 == entities[i].type && 0 == entities[i + 1].type) {
			auto v0 = pmp::Vertex(entities[i].entity);
			auto v1 = pmp::Vertex(entities[i + 1].entity);

			auto p = vv(v0, v1);
			h0 = p.first;
			preoph1 = p.second;
			assert(h0.is_valid());
			assert(sh.msh.to_vertex(h0) == v0);
			assert(sh.msh.to_vertex(preoph1) == v1);
		}
		else {
			if (0 == entities[i + 1].type) {
				// assert(0 != entities[i].type);
				auto tailv = pmp::Vertex(entities[i + 1].entity);
				preoph1 = gettailh(h0, tailv);
			}

			if (0 == entities[i].type) {
				if (0 == entities[i + 1].type) {
					assert(false);// should captured
				}
				else {
					auto headv = pmp::Vertex(entities[i].entity);
					auto tailoph = pmp::Halfedge(entities[i + 1].entity);
					h0 = getheadh(headv, tailoph);
				}
			}
		}


		auto f = sh.msh.face(h0);
		auto newh = sh.msh.insert_edge(h0, preoph1);
		auto f0 = sh.msh.face(newh);
		auto f1 = sh.msh.face(sh.msh.opposite_halfedge(newh));
		bsp.supportplaneprop[f0] = bsp.supportplaneprop[f];
		bsp.supportplaneprop[f1] = bsp.supportplaneprop[f];

		auto newe = sh.msh.edge(newh);
		bsp.edgeplanepair[newe].p0 = bsp.supportplaneprop[f];
		bsp.edgeplanepair[newe].p1 = newplid;

		popcutedge(newh);
	}

	const auto sz1 = testfun();

	auto fid0 = int(sh.msh.n_faces());
	{
		if (0)// got disabled
		{
			auto h = sh.msh.halfedge(cutedges.front(), 0);// this zero might be wrong!!!
			auto tov = sh.msh.to_vertex(h);

			std::vector<pmp::Halfedge> hs;
			hs.push_back(h);
			for (int i = 1; i < cutedges.size(); i++)// collect hs
			{
				auto e = cutedges[i];
				auto h0 = sh.msh.halfedge(e, 0);
				auto h1 = sh.msh.halfedge(e, 1);

				auto v1 = sh.msh.to_vertex(h0);
				auto v0 = sh.msh.to_vertex(h1);

				if (v0 == tov) {
					hs.push_back(h0);
					tov = v1;
				}
				else {
					assert(v1 == tov);
					hs.push_back(h1);
					tov = v0;
				}
			}

			if (0) {
				for (auto h : hs) {
					printedgeimpl(sh.msh, h);
				}
			}
		}

		auto hs = cuthalfedges;

		int offsetvid = int(sh.msh.n_vertices());
		int offsetfaceid = int(sh.msh.n_faces());
		int offsetedgeid = int(sh.msh.n_edges());
		sh.msh.cut_mesh(hs);

		const auto sz2 = testfun();// might not help issue is edge splitting???
		// assert(sz1 + 2 == sz2);

		//update vertex
		for (int i = 0; i < hs.size(); i++) {
			auto h = hs[i];

			auto i0 = (hs.size() - 1 + i) % hs.size();
			auto v0 = sh.msh.to_vertex(hs[i0]);

			bsp.hdpts[pmp::Vertex(offsetvid + i)] = bsp.hdpts[v0];

			auto e = pmp::Edge(offsetedgeid + i);
			auto newe = sh.msh.edge(h);
			bsp.edgeplanepair[e].p0 = bsp.edgeplanepair[newe].p0;
			bsp.edgeplanepair[e].p1 = bsp.edgeplanepair[newe].p1;
		}
	}
	bsp.supportplaneprop[pmp::Face(fid0)] = newplid;
	bsp.supportplaneprop[pmp::Face(fid0+1)] = newplid;

	return { 0, fid0 };
}






void PolyCut::init(int* p0, int* p1, int* p2)
{
	if (isjournal) {

		if (nullptr == fp) {
			fn = "";
			static int fncnt = 0;
			fn += std::to_string(fncnt++);
			fn += ".txt";

			fn = "pc.txt";
			fp = fopen(fn.c_str(), "wb");
		}

		// fprintf(fp, "%d %d %d %d %d %d %d %d %d\n",
		// 	p0[0], p0[1], p0[2],
		// 	p1[0], p1[1], p1[2],
		// 	p2[0], p2[1], p2[2]);

		TripleTriple t;
		std::copy(p0, p0 + 3, t.p0);
		std::copy(p1, p1 + 3, t.p1);
		std::copy(p2, p2 + 3, t.p2);
		fwrite(&t, sizeof(t), 1, fp);
	}

	{
		pts.clear();
		hnxt.clear();
		hpre.clear();
		htov.clear();
		planes.clear();
		hplanes.clear();
	}

	BspPlane pl;
	plane_from_points(p0, p1, p2, pl);
	planes.push_back(pl);

	BspPlane pl0, pl1, pl2;
	virtualtriangleplanes(p0, p1, p2, pl0, pl1, pl2);
	planes.push_back(pl0);
	planes.push_back(pl1);
	planes.push_back(pl2);

	hplanes.push_back(std::make_pair(0, 1));
	hplanes.push_back(std::make_pair(0, 2));
	hplanes.push_back(std::make_pair(0, 3));


	pts.push_back(tohd4(p0));
	pts.push_back(tohd4(p1));
	pts.push_back(tohd4(p2));
	p3d.push_back(vector3d(p0[0], p0[1], p0[2]));
	p3d.push_back(vector3d(p1[0], p1[1], p1[2]));
	p3d.push_back(vector3d(p2[0], p2[1], p2[2]));

	center = (p3d[0] + p3d[1] + p3d[2]) / 3.0;


	htov.push_back(1);
	htov.push_back(2);
	htov.push_back(0);

	hnxt.push_back(1);
	hnxt.push_back(2);
	hnxt.push_back(0);

	hpre.push_back(2);
	hpre.push_back(0);
	hpre.push_back(1);
}

void PolyCut::dumpoly(const std::string &fn) const
{
	pmp::SurfaceMesh msh;
	std::vector<double> crds3;
	std::vector<pmp::Vertex> vertlist;
	for (size_t i = 0; i < p3d.size(); i++)
	{
		crds3.push_back(p3d[i].x());
		crds3.push_back(p3d[i].y());
		crds3.push_back(p3d[i].z());

		pmp::Point p(p3d[i].x(), p3d[i].y(), p3d[i].z());
		vertlist.push_back(msh.add_vertex(p));
	}


	std::vector<pmp::Vertex> facevertlist;
	std::vector<bool> hflag;
	hflag.resize(hnxt.size());
	for (size_t i = 0; i < hnxt.size(); i++)
	{
		if (hflag[i])
			continue;

		facevertlist.clear();
		// hflag[i] = true;
		auto h = i;
		while (true)
		{
			if (hflag[h])
				break;

			facevertlist.push_back(vertlist[htov[h]]);
			hflag[h] = true;
			h = hnxt[h];
		}

		msh.add_face(facevertlist);
	}

	pmp::write_obj(msh, fn.c_str(), {});
}

int PolyCut::split(const BspPlane& pl, int starth, int newhs[2])
{
	if (isjournal) {
		fwrite(&pl, sizeof(pl), 1, fp);
		fwrite(&starth, sizeof(starth), 1, fp);
	}

	sh.clear();// collect sign changed edge;
	hssign.clear();
	edges.clear();

	int tov0 = htov[starth];
	int presign = classify_vertex(pl, pts[tov0]);

	int counter[3] = { 0, 0, 0 };
	int ntotal = 0;
	int crrh = hnxt[starth];
	for (;;)
	{
		int tov = htov[crrh];
		int sign = classify_vertex(pl, pts[tov]);
		counter[sign + 1]++;
		ntotal++;

		if (presign != sign) {
			sh.push_back(crrh);
			hssign.push_back(std::make_pair(presign, sign));
		}

		presign = sign;

		if (crrh == starth)
			break;
		crrh = hnxt[crrh];
	}

	auto setnxt = [&](int h0, int h1) {
		hnxt[h0] = h1;
		hpre[h1] = h0;
	};

	auto addnewedge = [&]() {
		int offset = int(hnxt.size());

		hnxt.resize(offset + 1);
		hpre.resize(offset + 1);
		htov.resize(offset + 1);
		hplanes.resize(offset + 1);
		return offset;
	};

	auto addnewvertex = [&]() {
		int offset = int(pts.size());
		pts.resize(offset + 1);
		p3d.resize(offset + 1);
		return offset;
	};

	auto setv = [&](int h, int tov) {
		htov[h] = tov;
	};

	if (counter[0] == 0&&counter[2]>0) { // 0 = -1+1; right
		return PC_LEFT;
	}
	if (counter[2] == 0&&counter[0]>0) {// 1 = 1 +1; left
		return PC_RIGHT;
	}
	if (counter[1] == ntotal) {
		return PC_ON;
	}

	if (!sh.empty())
	{
		// todo!!!
		// 3 mean one signe is zero and counter[0]>0 && counter[2]>0  1 points
		// 4 mean two sign are zeor and counter[0]>0 && counter[2]>0  zero points
		// 2 couter[0]==0 return left || counter[2]==0 return right
		assert(counter[0] > 0 || counter[1] > 0);
		// 2 counter[1]==0 and counter[0]>0 && counter[2]>0 cut generate two points
		if (3 == sh.size())
		{
			int i = 0;
			for (; i < sh.size(); i++)
			{
				const int i0 = i;
				const int i1 = (i + 1) % sh.size();

				if (hnxt[sh[i0]] == sh[i1] &&
					hssign[i0].second == 0 &&
					hssign[i1].first == 0)
				{
					
					break;
				}
			}
			assert(i < sh.size());
			for (int j = i; j < sh.size()-1; j++)
			{
				hssign[j] = hssign[j + 1];
				sh[j] = sh[j + 1];
			}

			sh.resize(2);
			hssign.resize(2);
		}

		
		if (4 == sh.size())
		{
			assert(counter[1] == 2 && counter[0] > 0 && counter[2] > 0);
			for (int i = 0; i < sh.size(); i++)
			{
				const int i0 = i;
				const int i1 = (i + 1) % sh.size();

				if (hssign[i0].second == 0 && hssign[i1].first == 0) {
					edges.push_back(sh[i0]);
					assert(hssign[i0].first + hssign[i1].second == 0);// opposite sign
				}
			}

			assert(edges.size() == 2);
		}
		else {
			if (2 != sh.size()) {
				closejournal();
			}
			assert(2 == sh.size());

			// 
			for (int i = 0; i < sh.size(); i++)
			{
				if (0 == hssign[i].first) {
					edges.push_back(hpre[sh[i]]);
				}
				else if (0 == hssign[i].second) {
					edges.push_back(sh[i]);
				}
				else {
					int newedge = addnewedge();
					int newv = addnewvertex();

					hplanes[newedge] = hplanes[sh[i]];

					auto pl0 = planes[hplanes[sh[i]].first];
					auto pl1 = planes[hplanes[sh[i]].second];
					intersect_3_planes(pl, pl0, pl1, pts[newv]);

					auto rwpl0 = bspplane2plane(pl);
					auto rwpl1 = bspplane2plane(pl0);
					auto rwpl2 = bspplane2plane(pl1);

					auto ipt0 = IntersectPlanesLU(rwpl0, rwpl1, rwpl2);
					p3d[newv] = ipt0;

					int n0 = hnxt[sh[i]];

					setnxt(sh[i], newedge);
					setnxt(newedge, n0);

					setv(newedge, htov[sh[i]]);
					setv(sh[i], newv);

					edges.push_back(sh[i]);
				}
			}
		}

		auto e0 = addnewedge();
		auto e1 = addnewedge();

		const int newplane = int(planes.size());
		planes.push_back(pl);
		{
			hplanes[e0] = (std::make_pair(newplane, 0));
			hplanes[e1] = hplanes[e0];
		}

		auto h0 = hnxt[edges[0]];
		auto h1 = hnxt[edges[1]];

		setnxt(edges[0], e0);
		setnxt(e0, h1);
		setv(e0, htov[edges[1]]);

		setnxt(e1, h0);
		setnxt(edges[1], e1);
		setv(e1, htov[edges[0]]);

		// assert(hssign[0].first == hssign[1].second);
		// assert(hssign[0].second == hssign[1].first);

		if (0 == hssign[0].first) {
			assert(0 != hssign[0].second);
			if (1 == hssign[0].second) {
				newhs[0] = e1;
				newhs[1] = e0;
			}
			else {
				newhs[0] = e0;
				newhs[1] = e1;
			}
		}
		else {
			if (1 == hssign[0].first) {
				newhs[0] = e0;
				newhs[1] = e1;
			}
			else {
				newhs[0] = e1;
				newhs[1] = e0;
			}
		}

		return PC_CUT;
	}

	return PC_NON;
}


struct BspSpliter
{
	Mesh& cubesh;
	Bsp& bsp;

	BspPlane& pl;
	int plid;
	int* p0;
	int* p1;
	int* p2;


	int leafsplit(int seedvid);

	bool isleftleaf(int nodeid) {
		return bsp.nodes[nodeid].isleft;
	}

	bool isrightleaf(int nodeid) {
		return bsp.nodes[nodeid].isright;
	}

	void split(int nodeid, int polyhid)
	{
		if (-1 == nodeid) {
			return;
			bsp.pc.closejournal();
		}
		const auto node = bsp.nodes[nodeid];
		auto pl = bsp.planes[node.plane];


		auto popleftchild = [&]() {
			if (node.left > -1) {
				int childid = leafsplit(node.left);
				if (childid > -1) {
					bsp.nodes[nodeid].left = childid;
					bsp.nodes[nodeid].isleft = false;
				}
			}
		};

		auto poprightchild = [&]() {
			if (node.right > -1) {
				int childid = leafsplit(node.right);
				if (childid > -1) {
					bsp.nodes[nodeid].right = childid;
					bsp.nodes[nodeid].isright = false;
				}
			}
		};

		int s0 = classify_vertex(pl, p0);
		int s1 = classify_vertex(pl, p1);
		int s2 = classify_vertex(pl, p2);

		int newpolyhs[2] = { 0, 0 };
		static int pcsplitcnt = 0;
		pcsplitcnt++;
		int pcstatus = bsp.pc.split(pl, polyhid, newpolyhs);// more test case one side of the polygon
		if (PC_ON == pcstatus)
			return;

		if (PC_RIGHT==pcstatus) {
			if (isrightleaf(nodeid)) {
				poprightchild();
			}
			else
				split(node.right, polyhid);
		}
		else if (PC_LEFT==pcstatus) {
			if (isleftleaf(nodeid)) {
				popleftchild();
			}
			else
				split(node.left, polyhid);
		}
		else
		{
			assert(PC_CUT == pcstatus);

			if (isleftleaf(nodeid)) {
				popleftchild();
			} else {
				split(node.left, newpolyhs[0]);
			}

			if (isrightleaf(nodeid)) {
				poprightchild();
			}
			else {
				split(node.right, newpolyhs[1]);
			}
		}
	}
};



int BspSpliter::leafsplit(int seedvid)
{
	auto splitfacetype = [&](const int f) {
		auto h = cubesh.msh.halfedge(pmp::Face(f));
		auto oph = cubesh.msh.opposite_halfedge(h);
		auto nxth = cubesh.msh.next_halfedge(oph);
		auto v0 = cubesh.msh.to_vertex(nxth);

		int type0 = classify_vertex(pl, bsp.hdpts[v0]);
		return type0;
	};

	auto getseedvid = [&](const int f) {
		auto h = cubesh.msh.halfedge(pmp::Face(f));
		auto v = cubesh.msh.to_vertex(h);

		return v.idx();
	};


	auto code = convexmesh_cut(bsp, cubesh, plid, pmp::Vertex(seedvid));
	if (0 == code.first) {
		int f0 = code.second;
		int f1 = f0 + 1;

		auto t0 = splitfacetype(f0);
		auto t1 = splitfacetype(f1);


		try
		{
			if (1) {
				Mesh sh0, sh1;
				auto f0list = collectfacefromseed(cubesh, pmp::Face(f0));
				auto f1list = collectfacefromseed(cubesh, pmp::Face(f1));

				std::string fn = std::to_string(plid);
				fn += "_f0.obj";
				dump(fn, cubesh.msh, f0list);
				fn = std::to_string(plid);
				fn += "_f1.obj";
				dump(fn, cubesh.msh, f1list);
			}
		}
		catch (const std::exception&)
		{
			printf("connectivity is not right!\n");
		}


		Node childnode;
		childnode.left = getseedvid(f0);
		childnode.right = getseedvid(f1);
		
		childnode.plane = plid;
		childnode.isleft = true;
		childnode.isright = true;

		if (t0 == -1) {
			assert(1 == t1);
			std::swap(childnode.left, childnode.right);
		}
		else {
			assert(-1 == t1);

		}

		int offset = int(bsp.nodes.size());
		bsp.nodes.push_back(childnode);
		return offset;
	}

	return -1;
}


struct TriVert {
	pmp::Vertex v0, v1, v2;
};

TriVert gettriverts(pmp::SurfaceMesh& msh, pmp::Face f)
{
	auto h0 = msh.halfedge(f);
	auto h1 = msh.next_halfedge(h0);
	auto h2 = msh.next_halfedge(h1);

	auto v0 = msh.to_vertex(h2);
	auto v1 = msh.to_vertex(h0);
	auto v2 = msh.to_vertex(h1);

	TriVert tv;
	tv.v0 = v0;
	tv.v1 = v1;
	tv.v2 = v2;

	return tv;
};

TriVert gettriverts(Mesh& sh, pmp::Face f) 
{
	return gettriverts(sh.msh, f);
};


void vertex_normals__(pmp::SurfaceMesh& mesh)
{
	auto vnormal = mesh.vertex_property<pmp::Normal>("v:normal");
	for (auto v : mesh.vertices())
		vnormal[v] = vertex_normal(mesh, v);
}

auto top__(pmp::Point& p)
{
	return vector3d(p[0], p[1], p[2]);
}

int cutterfaces[5][4] = {
	{0, 1, 2, 3},
	{1, 0, 4, 5},
	{2, 1, 5, 6},
	{3, 2, 6, 7},
	{0, 3, 7, 4}
};

void mesh2millbin2d(const Mesh& sh, const std::vector<Eigen::Vector3i>& intposlist,
	pmp::SurfaceMesh& outputmsh, std::vector<Eigen::Vector3i>& intcuterposlist)
{
	// for bunny 60 
	const double iscale = 1000;
	const double avglen = 60;
	const double avglenhalf = avglen / 2.0;
	auto msh = sh.msh;

	auto popintpos = [&](const vector3d& p_) {
		pmp::Point p(p_[0], p_[1], p_[2]);
		auto v = outputmsh.add_vertex(p);
		intcuterposlist.push_back({ int(p_[0]), int(p_[1]), int(p_[2]) });
	};

	// intcuterposlist.resize(msh.n_vertices() * 8);

	// update to the scale position
	for (auto v : msh.vertices())
	{
		msh.position(v)[0] = intposlist[v.idx()][0];
		msh.position(v)[1] = intposlist[v.idx()][1];
		msh.position(v)[2] = intposlist[v.idx()][2];
	}

	double rect = 17;
	const double recthalf = rect / 2.0;
	const double bboxscale = 1000;
	const int ucnt = 2 + int(bboxscale / rect);
	const int vcnt = 2 + int(bboxscale / rect);

	struct BinStatus
	{
		double mind = std::numeric_limits<double>::max();
		vector3d p;
		int triid = -1;
	};

	std::vector< BinStatus > bins;
	bins.resize(ucnt* vcnt);
	auto bindid = [&](const int i, const int j) {
		return i + j * ucnt;
	};
	auto udpatebin = [&](const int triid, const int i, const int j,
		const vector3d& p0, const vector3d& p1, const vector3d& p2)
	{
		auto& bin = bins[bindid(i, j)];

		double maxz = p0.z();
		maxz = std::max(maxz, p1.z());
		maxz = std::max(maxz, p2.z());
		const double depth = bboxscale - maxz;
		if (depth < bin.mind) {
			bin.mind = depth;
			bin.triid = triid;
		}
	};
	

	auto zaxis = vector3d(0, 0, 1);
	for (auto f : msh.faces())
	{
		auto tv = gettriverts(msh, f);


		auto p0 = top__(msh.position(tv.v0));
		auto p1 = top__(msh.position(tv.v1));
		auto p2 = top__(msh.position(tv.v2));

		auto d01 = p1 - p0;
		auto d02 = p2 - p0;
		auto n_ = d01.cross(d02);
		if (n_.dot(zaxis) < 0.0)
			continue;
	

		AABB aabb;
		aabb.expandBy(p0);
		aabb.expandBy(p1);
		aabb.expandBy(p2);

		int i0 = int(aabb.min[0] / rect);
		int j0 = int(aabb.min[1] / rect);
		int i1 = int(aabb.max[0] / rect) + 1;
		int j1 = int(aabb.max[1] / rect) + 1;

		for (int i = i0; i < i1; i++) {
			for (int j = j0; j < j1; j++) {
				udpatebin(int(f.idx()), i, j, p0, p1, p2);
			}
		}
	}

	int vidx = 0;
	for (int j = 0; j < vcnt; j++)
	{
		for (int i = 0; i < ucnt; i++)
		{
			auto& bin = bins[bindid(i, j)];

			double depth = bboxscale;
			if (std::numeric_limits<double>::max() != bin.mind)
				depth = bin.mind;

			double x = i * rect;
			double y = j * rect;
			double z = bboxscale - depth;

			vector3d p0(x, y, z);
			vector3d p1(x + rect, y, z);
			vector3d p2(x + rect, y + rect, z);
			vector3d p3(x, y + rect, z);

			vector3d n_(0, 0, 1);

			auto p4 = p0;// + n_ * 1000;
			auto p5 = p1;// + n_ * 1000;
			auto p6 = p2;// + n_ * 1000;
			auto p7 = p3;// + n_ * 1000;

			p4[2] = bboxscale;
			p5[2] = bboxscale;
			p6[2] = bboxscale;
			p7[2] = bboxscale;

			popintpos(p0);
			popintpos(p1);
			popintpos(p2);
			popintpos(p3);

			popintpos(p4);
			popintpos(p5);
			popintpos(p6);
			popintpos(p7);

			for (int j = 0; j < 5; j++)
			{
				pmp::Vertex v0(cutterfaces[j][0] + vidx * 8);
				pmp::Vertex v1(cutterfaces[j][1] + vidx * 8);
				pmp::Vertex v2(cutterfaces[j][2] + vidx * 8);
				pmp::Vertex v3(cutterfaces[j][3] + vidx * 8);

				outputmsh.add_quad(v0, v1, v2, v3);
			}
			vidx++;


		}
	}

	std::cout << "cutter total face:\t" << outputmsh.n_faces() << std::endl;
	return;
	// the cut logic need more check
	for (int j = 0; j < vcnt; j++)
	{
		for (int i = 0; i < ucnt; i++)
		{
			auto& bin = bins[bindid(i, j)];

			double depth = bboxscale;
			if (std::numeric_limits<double>::max() != bin.mind)
				depth = bin.mind;

			if (depth != bboxscale) {

				auto tv = gettriverts(msh, pmp::Face(bin.triid));


				auto p0 = top__(msh.position(tv.v0));
				auto p1 = top__(msh.position(tv.v1));
				auto p2 = top__(msh.position(tv.v2));

				auto p_ = (p0 + p1 + p2) / 3.0;

				auto d01 = p1 - p0;
				auto d02 = p2 - p0;
				auto n_ = d01.cross(d02);
				n_.normalize();

				// vector3d n_(n[0], n[1], n[2]), x, y, p_(p[0], p[1], p[2]);
				vector3d x, y;
				buildframe(n_, x, y);

				{
					auto p0 = p_ - recthalf * x - recthalf * y;
					auto p1 = p0 + rect * x;
					auto p2 = p1 + rect * y;
					auto p3 = p0 + rect * y;

					auto p4 = p0 + n_ * recthalf;
					auto p5 = p1 + n_ * recthalf;
					auto p6 = p2 + n_ * recthalf;
					auto p7 = p3 + n_ * recthalf;


					// duplicate:
					popintpos(p0);
					popintpos(p1);
					popintpos(p2);
					popintpos(p3);

					popintpos(p4);
					popintpos(p5);
					popintpos(p6);
					popintpos(p7);

					for (int j = 0; j < 5; j++)
					{
						pmp::Vertex v0(cutterfaces[j][0] + vidx * 8);
						pmp::Vertex v1(cutterfaces[j][1] + vidx * 8);
						pmp::Vertex v2(cutterfaces[j][2] + vidx * 8);
						pmp::Vertex v3(cutterfaces[j][3] + vidx * 8);

						outputmsh.add_quad(v0, v1, v2, v3);
					}
					vidx++;
				}
			}
		}
	}

	
}

void mesh2millcutter(const Mesh &sh, const std::vector<Eigen::Vector3i> &intposlist,
	pmp::SurfaceMesh &outputmsh, std::vector<Eigen::Vector3i> &intcuterposlist)
{
	// for bunny 60 
	const double iscale = 1000;
	const double avglen = 60;
	const double avglenhalf = avglen / 2.0;
	auto msh = sh.msh;

	// pmp::SurfaceMesh outputmsh;
	// std::vector<Eigen::Vector3i> intcuterposlist;

	intcuterposlist.resize(msh.n_vertices() * 8);

	// update to the scale position
	for (auto v : msh.vertices())
	{
		msh.position(v)[0] = intposlist[v.idx()][0];
		msh.position(v)[1] = intposlist[v.idx()][1];
		msh.position(v)[2] = intposlist[v.idx()][2];
	}

	vertex_normals__(msh);
	auto vnormal = msh.vertex_property<pmp::Normal>("v:normal");




	auto popintpos = [&](const vector3d &p_) {
		pmp::Point p(p_[0], p_[1], p_[2]);
		auto v  = outputmsh.add_vertex(p);
		intcuterposlist[v.idx()] = { int(p_[0]), int(p_[1]), int(p_[2]) };
	};

	for (auto v : msh.vertices())
	{
		auto p = msh.position(v);
		auto n = vnormal[v];
		vector3d n_(n[0], n[1], n[2]), x, y, p_(p[0], p[1], p[2]);
		buildframe(n_, x, y);

		auto p0 = p_ - avglenhalf * x - avglenhalf * y;
		auto p1 = p0 + avglen * x;
		auto p2 = p1 + avglen * y;
		auto p3 = p0 + avglen * y;

		auto p4 = p0 + n_ * 1000;
		auto p5 = p1 + n_ * 1000;
		auto p6 = p2 + n_ * 1000;
		auto p7 = p3 + n_ * 1000;
		
		popintpos(p0);
		popintpos(p1);
		popintpos(p2);
		popintpos(p3);

		popintpos(p4);
		popintpos(p5);
		popintpos(p6);
		popintpos(p7);

		for (int j = 0; j < 5; j++)
		{
			pmp::Vertex v0(cutterfaces[j][0] + v.idx() * 8);
			pmp::Vertex v1(cutterfaces[j][1] + v.idx() * 8);
			pmp::Vertex v2(cutterfaces[j][2] + v.idx() * 8);
			pmp::Vertex v3(cutterfaces[j][3] + v.idx() * 8);

			outputmsh.add_quad(v0, v1, v2, v3);
		}
	}
}

// https://stackoverflow.com/questions/8518743/get-directory-from-file-path-c
std::string getfolder(const std::string& filename)
{
	std::string directory;
	const size_t last_slash_idx = filename.rfind('\\');
	if (std::string::npos != last_slash_idx)
	{
		directory = filename.substr(0, last_slash_idx);
	}

	return directory;
}

int mesh2bsp(Mesh& sh_, bool cutter, bool iscuttercheckpoint, const char* outputfn)
{
	std::string directory;
	if (nullptr != outputfn) {
		std::string outputfnstr = outputfn;
		directory = getfolder(outputfnstr);
	}

	pmp::BoundingBox bb;
	for (auto v : sh_.msh.vertices())
		bb += sh_.msh.position(v);

	auto dir = bb.max() - bb.min();
	const int axis = max_axis(dir);

	int iscale = 1000;
	double scale = iscale;


	std::vector<Eigen::Vector3i> intposlist;
	intposlist.resize(sh_.msh.n_vertices());
	for (int i = 0; i < intposlist.size(); i++)
	{
		auto p0 = sh_.msh.position(pmp::Vertex(i));
		p0 = p0 - bb.min();

		for (int j = 0; j < 3; j++)
		{
			intposlist[i][j] = int(p0[j] / dir[axis] * scale);
		}
	}

	Mesh sh = sh_;
	if (cutter) {
		std::vector<Eigen::Vector3i> intcutterposlist;
		sh.msh = pmp::SurfaceMesh();
		// mesh2millcutter(sh_, intposlist, sh.msh, intcutterposlist);
		mesh2millbin2d(sh_, intposlist, sh.msh, intcutterposlist);
		intposlist = intcutterposlist;
	}

	if (0)
	{
		Mesh sh1 = sh;
		for (auto v : sh1.msh.vertices())
		{
			sh1.msh.position(v)[0] = intposlist[v.idx()][0];
			sh1.msh.position(v)[1] = intposlist[v.idx()][1];
			sh1.msh.position(v)[2] = intposlist[v.idx()][2];
		}

		sh1.write("im.obj");
	}


	std::vector<BspPlane> mshplanes(sh.msh.n_faces());
	{
		for (auto f : sh.msh.faces())
		{
			auto tv = gettriverts(sh, f);

			int* p0 = intposlist[tv.v0.idx()].data();
			int* p1 = intposlist[tv.v1.idx()].data();
			int* p2 = intposlist[tv.v2.idx()].data();

			plane_from_points(p0, p1, p2, mshplanes[f.idx()]);
		}
	}



	{
		static int totalfcnt = 0;
		Bsp bsp;
		Mesh cubesh;
		build_cube(bsp, cubesh, iscale);
		if (0)
			cubesh.write("imcube.obj");

		auto dumpsolidfaces = [&](std::string& fn)
		{
			std::vector<pmp::Face> solidfaces;
			for (auto& node : bsp.nodes)
			{
				if (node.isright) {
					auto seedv = pmp::Vertex(node.right);
					auto h = cubesh.msh.halfedge(seedv);
					auto f = cubesh.msh.face(h);

					auto convexfaces = collectfacefromseed(cubesh.msh, f);
					solidfaces.insert(solidfaces.end(), convexfaces.begin(), convexfaces.end());
				}
			}


			dump(fn, cubesh.msh, solidfaces, true);
		};

		for (auto f : sh.msh.faces())
		{
			totalfcnt++;
			const int fid = f.idx();
			if (0 == fid % 1000) {
				// printf("%d\n", fid);

				if (cutter&& iscuttercheckpoint)
				{

					std::string fn;
					if (!directory.empty())
						fn += directory;
					fn += "/solid";
					fn += std::to_string(fid);
					fn += ".obj";
					dumpsolidfaces(fn);
				}
			}
			try
			{

				auto tv = gettriverts(sh, f);

				int* p0 = intposlist[tv.v0.idx()].data();
				int* p1 = intposlist[tv.v1.idx()].data();
				int* p2 = intposlist[tv.v2.idx()].data();

				const int plid = int(bsp.planes.size());
				auto pl = mshplanes[f.idx()];
				if (pl.a == 0 && pl.b == 0 && pl.c == 0 && pl.d == 0)
					continue;

				bsp.planes.push_back(pl);

				bsp.pc.init(p0, p1, p2);

				BspSpliter spliter = { cubesh, bsp, pl, plid, p0, p1, p2 };
				spliter.split(0, 0);


				if (0)
					dumpmesh(int(f.idx()), cubesh);
			}
			catch (const std::exception&)
			{
				printf("---\n");
			}

		}


#if 0 // debug
		static int cutrltcnt = 10000;
		// dumpmesh(cutrltcnt++, cubesh);'


		static int cutrltcnt100 = 100000;

		if (cutter)
		{
			std::string fn = "solid";
			fn += std::to_string(cutrltcnt100);
			fn += ".obj";
			dumpsolidfaces(fn);
		}

		std::string fn;
		fn += std::to_string(cutrltcnt100++);
		fn += ".obj";

		dumpsolidfaces(fn);
#endif

		if (nullptr != outputfn) {
			std::string outputfnstr = outputfn;
			dumpsolidfaces(outputfnstr);
		}
	}


	return 0;
}

int convexmesh_cut(Bsp& bsp, Mesh& sh, BspPlane& pl)
{
	int newplid = int(bsp.planes.size());
	bsp.planes.push_back(pl);



	auto rtn = convexmesh_cut(bsp, sh, newplid);
	return rtn.first;
}




int bsp_cut(Bsp& bsp, Mesh& sh, BspPlane& pl)
{

	return 0;
}


/*
         
        /|\Z
         |
         |  /
         | /
		 |/
---------/--------->Y
		/|
	   / |
	 |/_ |
	 X   |
	     |

*/
int childoctantcenteroffset[8][3] = {
	{1, -1,  -1},
	{1, 1,   -1},
	{-1, 1,  -1},
	{-1, -1, -1},

	{1, -1,  1},
	{1, 1,   1},
	{-1, 1,  1},
	{-1, -1, 1},
};


void Octree::build(uint32_t parentid, std::vector<uint32_t>& trilist, TriFun fun)
{
	if (trilist.size() < Octant_MAX_TRI) {
		octantlist[parentid].childtris = std::move(trilist);
		return;
	}
		

	std::vector< std::vector<uint32_t> > octanttrislist(8);

	auto parentnode = octantlist[parentid];
	const float childhalfsz = parentnode.extent / 2.0f;

	float boxhalfsize[3] = { childhalfsz , childhalfsz , childhalfsz };

	auto getchildcenter = [&](const int j, float childcenter[3]) {
		childcenter[0] = parentnode.x + childhalfsz * childoctantcenteroffset[j][0];
		childcenter[1] = parentnode.y + childhalfsz * childoctantcenteroffset[j][1];
		childcenter[2] = parentnode.z + childhalfsz * childoctantcenteroffset[j][2];
	};
	
	//TODO seperate axis
	for (auto triid : trilist)
	{
		float triverts[3][3];
		fun(triid, triverts[0], triverts[1], triverts[2]);

		for (int j = 0; j < 8; j++)
		{
			float childcenter[3];
			getchildcenter(j, childcenter);
			
			if (1 == triBoxOverlap(childcenter, boxhalfsize, triverts)) {
				octanttrislist[j].push_back(triid);
			}
		}
	}

	const int offset = int(octantlist.size());
	octantlist.resize(octantlist.size() + 8);

	octantlist[parentid].childid = offset;
	for (int j = 0; j < 8; j++)
	{
		getchildcenter(j, &(octantlist[offset + j].x));
		octantlist[offset + j].extent = childhalfsz;

		const int childid = offset + j;
		build(childid, octanttrislist[j], fun);
	}
}





//---

std::vector<pmp::Face> collectfacefromseed(pmp::SurfaceMesh& msh, pmp::Face f)
{
	std::vector<pmp::Face> faces;
	std::vector<bool> faceflag(msh.n_faces());
	faceflag[f.idx()] = true;
	faces.push_back(f);

	for (size_t i = 0; i < faces.size(); i++)
	{
		auto f = faces[i];
		auto h = msh.halfedge(f);
		auto crrh = h;
		while (true)
		{
			auto oph = msh.opposite_halfedge(crrh);
			auto opf = msh.face(oph);

			if (!faceflag[opf.idx()])
			{
				faceflag[opf.idx()] = true;
				faces.push_back(opf);
			}

			crrh = msh.next_halfedge(crrh);
			if (h == crrh)
				break;
		}
	}

	return faces;
}

std::vector<pmp::Face> collectfacefromseed(Mesh& sh, pmp::Face f)
{
	return collectfacefromseed(sh.msh, f);
}

std::vector<pmp::Face> collectfacefromseedv(Mesh &sh, pmp::Vertex v3)
{
	auto r = sh.msh.halfedges(v3);

	std::vector<pmp::Face> faces;
	for (auto h : r) {
		auto f = sh.msh.face(h);
		faces.push_back(f);
	}

	return collectfacefromseed(sh, faces.front());
}


void dump(std::string fn, pmp::SurfaceMesh& msh, std::vector<pmp::Face>& faces, bool isdump)
{
	pmp::SurfaceMesh oputmesh;
	std::vector<pmp::Vertex> nodemap(msh.n_vertices());
	std::vector<pmp::Vertex> facevertices;
	for (auto f : faces)
	{
		auto h = msh.halfedge(f);

		auto crrh = h;
		for (;;)
		{
			auto v = msh.to_vertex(crrh);
			if (!nodemap[v.idx()].is_valid()) {
				auto p = msh.position(v);
				nodemap[v.idx()] = oputmesh.add_vertex(p);
			}

			facevertices.push_back(nodemap[v.idx()]);

			crrh = msh.next_halfedge(crrh);
			if (crrh == h)
				break;
		}

		oputmesh.add_face(facevertices);
		facevertices.clear();
	}

	if (isdump)
		pmp::write_obj(oputmesh, fn, {});
}


void dumpmesh(Mesh& sh, const std::string& fn__, const int i, bool isdump)
{
	int zonecnt = 0;
	std::vector<bool> facevis(sh.msh.n_faces());
	for (auto f : sh.msh.faces())
	{
		if (facevis[f.idx()])
			continue;

		std::vector<pmp::Face> fq;
		fq.push_back(f);
		facevis[f.idx()] = true;

		for (int i = 0; i < fq.size(); i++)
		{
			auto backf = fq[i];

			auto h = sh.msh.halfedge(backf);
			auto crrh = h;
			for (;;)
			{
				auto opcrrh = sh.msh.opposite_halfedge(crrh);
				auto opf = sh.msh.face(opcrrh);
				if (!facevis[opf.idx()])
				{
					facevis[opf.idx()] = true;
					fq.push_back(opf);
				}

				crrh = sh.msh.next_halfedge(crrh);
				if (h == crrh)
					break;
			}
		}

		std::string fn = fn__;// "cut";
		fn += std::to_string(i);
		fn += "_";
		fn += std::to_string(zonecnt++);
		fn += ".obj";


		dump(fn, sh.msh, fq, isdump);
	}
}

void dumpmesh(const int i, Mesh& sh)
{
	std::string fn = "cut";
	fn += std::to_string(i);
	fn += ".obj";
	sh.write(fn);

	if (0)
	{
		fn = "cut";
		dumpmesh(sh, fn,  i, false);
	}
}


void dumpplane(const BspPlane& pl, const Eigen::Vector3d& refpt, const double r)
{
	pmp::SurfaceMesh msh;
	auto rq = plane2quad(pl, refpt, r);

	std::vector<pmp::Vertex> vertices;
	vertices.push_back(addvertex0(msh, rq.p0.data()));
	vertices.push_back(addvertex0(msh, rq.p1.data()));
	vertices.push_back(addvertex0(msh, rq.p2.data()));
	vertices.push_back(addvertex0(msh, rq.p3.data()));

	msh.add_face(vertices);

	pmp::write_obj(msh, "plane0.obj", {});
}