#pragma once


#include <pmp/surface_mesh.h>


struct Mesh
{
	pmp::SurfaceMesh msh;

    void load(std::string fn);
	void write(std::string fn);
};




