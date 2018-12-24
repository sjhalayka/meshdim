#ifndef MESH_H
#define MESH_H

#include "primitives.h"

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <iomanip>
using std::setiosflags;

#include <ios>
using std::ios_base;
using std::ios;

#include <set>
using std::set;

#include <vector>
using std::vector;

#include <limits>
using std::numeric_limits;

#include <cstring> // for memcpy()
#include <cctype>


class indexed_mesh
{
public:
	vector<vertex_3> vertices;
	vector<indexed_triangle> triangles;
	vector< vector<size_t> > vertex_to_vertex_indices;
	vector< vector<size_t> > vertex_to_triangle_indices;
	vector<vertex_3> vertex_normals;
	vector<vertex_3> triangle_normals;

	void clear(void);
	bool operator==(const indexed_mesh &right);
	bool operator!=(const indexed_mesh &right);

	void fix_cracks(void);

	bool load_from_binary_stereo_lithography_file(const char *const file_name, const bool generate_normals = true, const size_t buffer_width = 65536);
	bool save_to_binary_stereo_lithography_file(const char *const file_name, const size_t buffer_width = 65536);
	bool save_to_povray_mesh2_file(const char *const file_name, const bool write_vertex_normals = false);

	float get_x_extent(void);
	float get_y_extent(void);
	float get_z_extent(void);
	float get_max_extent(void);
	float set_max_extent(const float target_max_extent);
	float get_triangle_area(const size_t tri_index);
	float get_vertex_neighbourhood_area(const size_t vertex_index);
	float get_area(void);
	float set_area(const float target_area);
	float get_triangle_volume(const size_t tri_index);
	float get_volume(void);
	float set_volume(const float target_volume);

	// See: Geometric Signal Processing on Polygonal Meshes by G. Taubin
	void laplace_smooth(const float scale);
	void laplace_smooth_angle(const float scale);
	void laplace_smooth_fujiwara(const float scale);
	void laplace_smooth_cn(const float scale);

	void taubin_smooth(const float lambda, const float mu, const size_t steps);
	void taubin_smooth_angle(const float lambda, const float mu, const size_t steps);
	void taubin_smooth_fujiwara(const float lambda, const float mu, const size_t steps);
	void taubin_smooth_cn(const float lambda, const float mu, const size_t steps);

	void get_degenerate_triangles(vector<size_t> &degenerates);
	size_t get_degenerate_triangle_count(void);

	void get_tri_neighbours(const size_t tri_index, vector<size_t> &neighbours)
	{
		set<size_t> temp_neighbours;

		vector<size_t> tri_indices;

		get_triangles_shared_by_vertex_pair(triangles[tri_index].vertex_indices[0], triangles[tri_index].vertex_indices[1], tri_indices);

		for(size_t i = 0; i < tri_indices.size(); i++)
			if(tri_indices[i] != tri_index)
				temp_neighbours.insert(tri_indices[i]);

		get_triangles_shared_by_vertex_pair(triangles[tri_index].vertex_indices[0], triangles[tri_index].vertex_indices[2], tri_indices);

		for(size_t i = 0; i < tri_indices.size(); i++)
			if(tri_indices[i] != tri_index)
				temp_neighbours.insert(tri_indices[i]);

		get_triangles_shared_by_vertex_pair(triangles[tri_index].vertex_indices[1], triangles[tri_index].vertex_indices[2], tri_indices);

		for(size_t i = 0; i < tri_indices.size(); i++)
			if(tri_indices[i] != tri_index)
				temp_neighbours.insert(tri_indices[i]);
		
		neighbours.clear();
		
		for(set<size_t>::const_iterator ci = temp_neighbours.begin(); ci != temp_neighbours.end(); ci++)
			neighbours.push_back(*ci);
	}

	vertex_3 get_tri_normal(size_t i)
	{
		vertex_3 v0 = vertices[triangles[i].vertex_indices[1]] - vertices[triangles[i].vertex_indices[0]];
		vertex_3 v1 = vertices[triangles[i].vertex_indices[2]] - vertices[triangles[i].vertex_indices[0]];
		vertex_3 normal = v0.cross(v1);
		normal.normalize();

		return normal;
	}

	void colour_tris(vector<vertex_3> &tri_colours)
	{
		tri_colours.resize(triangles.size());

		for(size_t i = 0; i < triangles.size(); i++)
		{
			vector<size_t> neighbours;
			get_tri_neighbours(i, neighbours);

			// instead of getting edge-based neighbours,
			// get vertex-based neighbours

			bool flag = false;

			vertex_3 tri_normal = get_tri_normal(i);

			for(size_t j = 0; j < neighbours.size(); j++)
			{
				vertex_3 neighbour_normal = get_tri_normal(neighbours[j]);

				float dot = tri_normal.dot(neighbour_normal);

				cout << i << ' ' << neighbours[j] << ' ' << dot << endl;

				if(dot < 0.0f)
				{
					flag = true;
//					break;
				}
			}

			if(true == flag)
			{
				tri_colours[i].x = 1.0f;
				tri_colours[i].y = 0.0f;
				tri_colours[i].z = 0.0f;
			}
			else
			{
				tri_colours[i].x = 1.0f;
				tri_colours[i].y = 1.0f;
				tri_colours[i].z = 1.0f;
			}
		}
	}







	float get_vertex_angle_total(size_t vertex_index)
	{
		float angle = 0;

		if(vertex_index >= vertices.size())
			return angle;

		// Get sum of angle at vertex[i] (also try doing the ring around vertex[i])
		for(size_t i = 0; i < vertex_to_triangle_indices[vertex_index].size(); i++)
		{
			size_t tri_index = vertex_to_triangle_indices[vertex_index][i];

			size_t v0, v1;

			bool found_first_vertex = false;

			for(size_t j = 0; j < 3; j++)
			{
				if(vertex_index != triangles[tri_index].vertex_indices[j])
				{
					if(false == found_first_vertex)
					{
						v0 = triangles[tri_index].vertex_indices[j];
						found_first_vertex = true;
					}
					else
					{
						v1 = triangles[tri_index].vertex_indices[j];
						break;
					}
				}
			}

			vertex_3 a = vertices[v0] - vertices[vertex_index];
			vertex_3 b = vertices[v1] - vertices[vertex_index];
			a.normalize();
			b.normalize();

			float dot = a.dot(b);

			if(-1 > dot)
				dot = -1;
			else if(1 < dot)
				dot = 1;

			angle += acosf(dot);
		}

		return angle;
	}

protected:
	void generate_vertex_normals(void);
	void generate_triangle_normals(void);
	void generate_vertex_and_triangle_normals(void);
	void regenerate_vertex_and_triangle_normals_if_exists(void);
	template<typename T> void eliminate_vector_duplicates(vector<T> &v);
	bool merge_vertex_pair(const size_t keeper, const size_t goner);
	void get_triangles_shared_by_vertex_pair(const size_t v0, const size_t v1, vector<size_t> &triangle_indices);
	float get_angle_opposite_of_vertex_pair(const size_t v0, const size_t v1, const size_t tri_index);
};


#endif
