#ifndef main_H
#define main_H

#include "primitives.h"
#include "mesh.h"
#include "string_utils.h"

#include <cstdlib>

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;
using std::istringstream;


indexed_mesh mesh;


vector< vector<size_t> > tri_neighbours;
vector<vertex_3> tri_normals;

size_t curr_tri_index = 0;




#endif
