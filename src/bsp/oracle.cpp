#include "bsp.h"

#include <slimcpplib/long_int.h>
#include "plywriter.h"
#include "pmp/io/write_obj.h"


using namespace slim::literals;
using uint128_t = slim::uint128_t;
using int128_t = slim::int128_t;
using uint256_t = slim::uint256_t;
using int256_t = slim::int256_t;


void printedgeimpl(pmp::SurfaceMesh &msh, pmp::Halfedge h)
{
	static int edgecnt = 0;
	std::string fn = "e_";
	fn += std::to_string(edgecnt++);
	fn += ".ply";

	std::vector<double> crds3;
	std::vector<int> edges2;

	auto h1 = msh.opposite_halfedge(h);
	auto v0 = msh.to_vertex(h1);
	auto v1 = msh.to_vertex(h);

	auto p0 = msh.position(v0);
	auto p1 = msh.position(v1);

	crds3.push_back(p0[0]);
	crds3.push_back(p0[1]);
	crds3.push_back(p0[2]);

	crds3.push_back(p1[0]);
	crds3.push_back(p1[1]);
	crds3.push_back(p1[2]);

	edges2.push_back(0);
	edges2.push_back(1);

	writeplot(fn.c_str(), crds3, edges2);
};

int classify_vertex(const BspPlane& s, int* p)
{
	auto hd = tohd4(p);

	return classify_vertex(s, hd);
}

#define LOWBIT 0
int classify_vertex(const BspPlane& s, const Homogeneous4DPoint& p)
{
#if LOWBIT
	int64_t d0 = s.a * p.x1 + s.b * p.x2 + s.c * p.x3 + s.d * p.x4;

#if 0
	double dd0 = double(d0);
	double d1 = (dd0 / p.x4 / std::sqrt(s.a * s.a + s.b * s.b + s.c * s.c));
#endif

	return sgn(d0) * sgn(p.x4);
#else // LOWBIT

	int128_t sval = s.a;
	int128_t pval = p.x1;

	int128_t d0 = sval * pval;

	sval = s.b;
	pval = p.x2;
	d0 += sval * pval;

	sval = s.c;
	pval = p.x3;
	d0 += sval * pval;

	sval = s.d;
	pval = p.x4;
	d0 += sval * pval;

	return sgn(d0) * sgn(p.x4);

#endif
}


int virtualtriangleplanes(int* p0, int* p1, int* p2,
	BspPlane& pl0, BspPlane& pl1, BspPlane& pl2)
{
	BspPlane plane;
	plane_from_points(p0, p1, p2, plane);
	vector3d n(plane.a, plane.b, plane.c);
	n.normalize();

	vector3d v0(p0[0], p0[1], p0[2]);
	vector3d v1(p1[0], p1[1], p1[2]);
	vector3d v2(p2[0], p2[1], p2[2]);

	vector3d v01 = v1 - v0;
	vector3d v12 = v2 - v1;
	vector3d v20 = v0 - v2;

	double argl = (v01.norm() + v12.norm() + v20.norm()) / 3.0;

	vector3d center = (v0 + v1 + v2) / 3.0;

	vector3d tip = center + argl * n;
	Vector3i itip(int(tip.x()), int(tip.y()), int(tip.z()));

	if (0)
	{

		auto truncate = [](int val) {
			if (val > 1000)return 999;
			if (val < 0) return 0;
			return val;
		};
		itip[0] = truncate(itip[0]);
		itip[1] = truncate(itip[1]);
		itip[2] = truncate(itip[2]);
	}

	assert(1 == classify_vertex(plane, itip.data()));

	plane_from_points(p0, p1, itip.data(), pl0);
	plane_from_points(p1, p2, itip.data(), pl1);
	plane_from_points(p2, p0, itip.data(), pl2);


	if (0)
	{
		pmp::SurfaceMesh msh0;
		auto i0 = addvertex0(msh0, v0.data());
		auto i1 = addvertex0(msh0, v1.data());
		auto i2 = addvertex0(msh0, v2.data());
		auto dtp = tip;
		dtp[0] = itip[0];
		dtp[1] = itip[1];
		dtp[2] = itip[2];
		auto i3 = addvertex0(msh0, dtp.data());
		msh0.add_triangle(i1, i0, i2);
		msh0.add_triangle(i0, i1, i3);
		msh0.add_triangle(i1, i2, i3);
		msh0.add_triangle(i2, i0, i3);

		pmp::write_obj(msh0, "t0.obj", {});
	}

	return 0;
}

int planelinesegmentintersection(int* p0, int* p1, const BspPlane& pl, Homogeneous4DPoint &hd4)
{
	Vector3i p0_(p0[0], p0[1], p0[2]);
	Vector3i p1_(p1[0], p1[1], p1[2]);

	Vector3i d01 = p0_ - p1_;
	Vector3i n(pl.a, pl.b, pl.c);

	int64_t l1 = n.dot(p0_) + pl.d;
	int64_t l2 = n.dot(d01);

	Vector3i p = l2 * p0_ + l1 * d01;

	hd4.x1 = p[0];
	hd4.x2 = p[1];
	hd4.x3 = p[2];
	hd4.x4 = l2;
	return 0;
}


enum CutType
{
	CT_VERT = 0,
	CT_EDGE = 1
};

struct ScopedCleanVertexEdgeFalg
{
	Bsp& bsp;
	Mesh& sh;

	std::vector<pmp::Vertex> taggedvertex;
	std::vector<pmp::Edge> taggededge;
	std::vector<pmp::Face> taggedfaces;
	
	~ScopedCleanVertexEdgeFalg()
	{
		for (auto v : taggedvertex)
		{
			bsp.vertexflag[pmp::Vertex(v)] = 0;
		}

		for (auto e : taggededge)
		{
			bsp.edgeflag[pmp::Edge(e)] = false;
		}

		for (auto f : taggedfaces)
			bsp.faceflag[f] = false;
	}

	void mark(pmp::Vertex v) {
		bsp.vertexflag[v] = 1;
		taggedvertex.push_back(v);
	}

	void mark(pmp::Edge e) {
		bsp.edgeflag[e] = true;
		taggededge.push_back(e);
	}

	void mark(pmp::Halfedge h) {
		auto e = sh.msh.edge(h);
		mark(e);
	}

	void mark(pmp::Face f) {
		bsp.faceflag[f] = true;
		taggedfaces.push_back(f);
	}

	bool ismarked(pmp::Face f) {
		return bsp.faceflag[f];
	}

	bool ismarked(pmp::Vertex v) {
		return bsp.vertexflag[v]==1;
	}

	bool ismarked(pmp::Edge e) {
		return bsp.edgeflag[e];
	}

	bool ismarked(pmp::Halfedge h) {
		auto e = sh.msh.edge(h);
		return ismarked(e);
	}
};
std::vector< CutEntity> convexmeshcut_traversal(Bsp& bsp, Mesh& sh, BspPlane& pl, pmp::Vertex v_, pmp::Halfedge nxth_, int type)
{
	ScopedCleanVertexEdgeFalg evflag = { bsp, sh };
	std::vector< CutEntity> entities;

	auto printedge = [&](pmp::Halfedge h) {
		if (h.is_valid())
			printedgeimpl(sh.msh, h);

		std::vector<double> crds3;
		for (int i = 0; i < entities.size(); i++)
		{
			if (0 == entities[i].type) {
				auto v = pmp::Vertex(entities[i].entity);
				auto p = sh.msh.position(v);
				crds3.push_back(p[0]);
				crds3.push_back(p[1]);
				crds3.push_back(p[2]);
			}
			else {
				auto e = sh.msh.edge(pmp::Halfedge(entities[i].entity));
				auto plid0 = bsp.edgeplanepair[e].p0;
				auto plid1 = bsp.edgeplanepair[e].p1;

				auto pl0 = bsp.planes[plid0];
				auto pl1 = bsp.planes[plid1];
				auto pl2 = pl;

				auto rwpl0 = bspplane2plane(pl0);
				auto rwpl1 = bspplane2plane(pl1);
				auto rwpl2 = bspplane2plane(pl2);

				auto ipt0 = IntersectPlanesLU(rwpl0, rwpl1, rwpl2);
				crds3.push_back(ipt0[0]);
				crds3.push_back(ipt0[1]);
				crds3.push_back(ipt0[2]);
			}
		}

		std::vector<int> edges2;
		for (int i = 0; i < int(entities.size()) - 1; i++)
		{
			edges2.push_back(i);
			edges2.push_back(i + 1);
		}

		if (!edges2.empty())
			writeplot("cut.ply", crds3, edges2);
	};

	int prevtype_ = 0;
	if (0 != type) {
		auto preh = sh.msh.prev_halfedge(nxth_);
		auto prev = sh.msh.to_vertex(preh);

		prevtype_ = classify_vertex(pl, bsp.hdpts[prev]);
	}

	pmp::Face preface;
	
	auto getnxth = [&](pmp::Halfedge nxth,  pmp::Vertex v) {

		bool isfromvert = false;
		int prevtype = nxth.is_valid() ? 1 : 0;
		if (0 == prevtype) {
			isfromvert = true;
			assert(!nxth.is_valid());

			auto range = sh.msh.halfedges(v);
			int pretype = -1;
			bool initpretype = false;


			auto starth = sh.msh.halfedge(v);
			auto h = starth;
			//for (auto h : range)
			while (true)
			{
				auto v1 = sh.msh.to_vertex(h);

				const int type1 = classify_vertex(pl, bsp.hdpts[v1]);

				if (!initpretype) {
					initpretype = true;
					pretype = type1;
				}
				else if (
					preface!=sh.msh.face(h)&&
					!evflag.ismarked(sh.msh.edge(h))&&
					!evflag.ismarked(sh.msh.face(h))) {
					if (type1 != pretype && 0 != type1 && 0 != pretype)
					{
						nxth = h;
						break;
					}
				}

				// auto preh = sh.msh.prev_halfedge(h);
				// auto opppreh = sh.msh.opposite_halfedge(preh);
				// h = opppreh;
				auto oph = sh.msh.opposite_halfedge(h);
				h = sh.msh.next_halfedge(oph);
				pretype = type1;// need update skip case
				if (starth == h)
					break;
			}

			if (!nxth.is_valid())
			{
				const auto startv = sh.msh.to_vertex(starth);
				const auto startvtype = classify_vertex(pl, bsp.hdpts[startv]);
				if (startvtype != pretype && 0 != startvtype && 0 != pretype)
				{
					if (!evflag.ismarked(sh.msh.face(h)))
					{
						// assert(preface != sh.msh.face(h));
						nxth = h;// circle back
					}

				}
				
				if (!nxth.is_valid())
				{
					// assert(false);// return pmp::Halfedge();// convext cut this plane coninside exiting face
					for (auto h : range)
					{
						auto v1 = sh.msh.to_vertex(h);
						const int type1 = classify_vertex(pl, bsp.hdpts[v1]);
						if (0 == type1) {
							if (v1.idx() != entities.back().entity&&!evflag.ismarked(h))// check not bite back
							{
								evflag.mark(h);
								return std::make_pair(pmp::Halfedge(), v1);
							}
						}
					}
				}

				if (!nxth.is_valid()) {
					// assert(false);
					return std::make_pair(pmp::Halfedge(), pmp::Vertex());
				}
			}
		}
		else {
			assert(!v.is_valid());
		}

		evflag.mark(nxth);
		auto nxtv = sh.msh.to_vertex(nxth);
		prevtype = classify_vertex(pl, bsp.hdpts[nxtv]);
		assert(nxth.is_valid());
		auto crrh = sh.msh.next_halfedge(nxth);
		while (true)
		{
			if (crrh == nxth) {
				assert(false);
				break;
			}
			nxtv = sh.msh.to_vertex(crrh);
			auto nxttype = classify_vertex(pl, bsp.hdpts[nxtv]);
			assert(prevtype != 0);
			if (nxttype != prevtype) {

				evflag.mark(sh.msh.face(crrh));

				if (0 == nxttype) {// point cut
					return std::make_pair(pmp::Halfedge(), nxtv);
				}
				else {
					pmp::Halfedge opcrrh = sh.msh.opposite_halfedge(crrh);
					return std::make_pair(opcrrh, pmp::Vertex());
				}
			}

			crrh = sh.msh.next_halfedge(crrh);

		}

		return std::make_pair(pmp::Halfedge(), pmp::Vertex());
	};

	auto startv = v_;
	auto starth = nxth_;

	auto crrh = starth;
	auto crrv = startv;

	int isprintedge = 0;
	for (;;)
	{
		if (isprintedge)
			printedge(crrh);

		if (crrh.is_valid()) {
			assert(!crrv.is_valid());

			CutEntity ent = { int(crrh.idx()), CT_EDGE };
			entities.push_back(ent);
			// evflag.mark(sh.msh.edge(crrh));

		}
		else if (crrv.is_valid()) {
			assert(!crrh.is_valid());

			CutEntity ent = { int(crrv.idx()), CT_VERT };
			entities.push_back(ent);
			// evflag.mark(crrv);
		}

		auto nxt = getnxth(crrh, crrv);

		if (crrh.is_valid())
			preface = sh.msh.face(crrh);
		else
			preface = pmp::Face();

		crrh = nxt.first;
		crrv = nxt.second;
		if (!crrh.is_valid() && !crrv.is_valid())
			return {};

		if ((starth.is_valid()&&crrh == starth) || (startv.is_valid()&&crrv == startv))
			break;
	}

	entities.push_back(entities.front());
	return entities;
}