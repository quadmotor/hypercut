#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <vector>


inline void writeplyTQ(const char *fn, std::vector<double> &crds3,
	std::vector<int> &tris3, std::vector<int>& quads4, std::vector<int> &edges2)
{
	FILE* fp = fopen(fn, "w");

	const size_t nn = size_t(crds3.size()/3);
	const size_t nt = size_t(tris3.size()/3);
	const size_t nq = size_t(quads4.size() / 4);
	const size_t nedge = size_t(edges2.size() / 2);;

	const size_t nf = nt + nq;

	fprintf(fp, "ply\n");
	fprintf(fp, "format ascii 1.0\n");
	fprintf(fp, "comment author : Greg Turk\n");
	fprintf(fp, "comment object : another cube\n");
	fprintf(fp, "element vertex %d\n", int(nn));
	fprintf(fp, "property float x\n");
	fprintf(fp, "property float y\n");
	fprintf(fp, "property float z\n");
	fprintf(fp, "element face %d\n", int(nf));
	fprintf(fp, "property list uchar int vertex_index\n");
	fprintf(fp, "element edge %d\n", int(nedge));
	fprintf(fp, "property int vertex1\n");
	fprintf(fp, "property int vertex2\n");
	fprintf(fp, "end_header\n");

	for (size_t i = 0; i < nn; i++)
	{
		fprintf(fp, "%lf %lf %lf\n", crds3[i * 3 + 0], crds3[i * 3 + 1], crds3[i * 3 + 2]);
	}

	for (size_t i = 0; i < nt; i++)
	{
		fprintf(fp, "3 %d %d %d\n", tris3[i * 3 + 0], tris3[i * 3 + 1], tris3[i * 3 + 2]);
	}

	for (size_t i = 0; i < nq; i++)
	{
		fprintf(fp, "4 %d %d %d %d\n", 
			quads4[i * 4 + 0], 
			quads4[i * 4 + 1], 
			quads4[i * 4 + 2],
			quads4[i * 4 + 3]);
	}

	for (size_t i = 0; i < nedge; i++)
	{
		fprintf(fp, "%d %d\n", edges2[i * 2 + 0], edges2[i * 2 + 1]);
	}


	fclose(fp);


}


inline void writeply(const char* fn, std::vector<double>& crds3,
	std::vector<int>& tris3, std::vector<int>& edges2)
{
	std::vector<int> quads4;
	writeplyTQ(fn, crds3, tris3, quads4, edges2);
}

inline void writeplot(const char* fn, std::vector<double>& crds3, std::vector<int>& edges2)
{
	std::vector<int> tris3;
	writeply(fn, crds3, tris3, edges2);
}
/*
ply
format ascii 1.0
comment author : Greg Turk
comment object : another cube
element vertex 8
property float x
property float y
property float z
property red uchar { start of vertex color }
property green uchar
property blue uchar
element face 7
property list uchar int vertex_index { number of vertices for each face }
element edge 5 { five edges in object }
property int vertex1 { index to first vertex of edge }
property int vertex2 { index to second vertex }
property uchar red { start of edge color }
property uchar green
property uchar blue end_header
0 0 0 255 0 0 { start of vertex list }
0 0 1 255 0 0
0 1 1 255 0 0
0 1 0 255 0 0
1 0 0 0 0 255
1 0 1 0 0 255
1 1 1 0 0 255
1 1 0 0 0 255
3 0 1 2 { start of face list, begin with a triangle }
3 0 2 3 { another triangle }
4 7 6 5 4 { now some quadrilaterals }
4 0 4 5 1
4 1 5 6 2
4 2 6 7 3
4 3 7 4 0
0 1 255 255 255 { start of edge list, begin with white edge }
1 2 255 255 255
2 3 255 255 255
3 0 255 255 255
2 0 0 0 0 { end with a single black line }
*/