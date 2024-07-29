#include "mesh.h"
#include <pmp/io/read_obj.h>
#include <pmp/io/write_obj.h>


void Mesh::load(std::string fn)
{
    pmp::read_obj(this->msh, fn);
}


void Mesh::write(std::string fn)
{
    pmp::write_obj(this->msh, fn, {});
}


