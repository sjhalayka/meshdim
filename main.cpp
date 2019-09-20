#include "main.h"

int main(int argc, char **argv)
{
	// Check for commandline arguments
	if(2 != argc)
	{
		cout << "Example usage: meshdim filename.stl" << endl;
		return 1;
	}

	// Create a mesh
	indexed_mesh mesh;

	// Load from binary STL file
	if(false == mesh.load_from_binary_stereo_lithography_file(argv[1]))
	{
		cout << "Error: Could not properly read file " << argv[1] << endl;
		return 2;
	}

	// Store triangle neighbour indices for later use
	vector< vector<size_t> > tri_neighbours;

	// Store triangle normals for later use
	vector<vertex_3> tri_normals;

	tri_neighbours.resize(mesh.triangles.size());
	tri_normals.resize(mesh.triangles.size());

	// For each triangle in the mesh
	for(size_t i = 0; i < mesh.triangles.size(); i++)
	{
		mesh.get_tri_neighbours(i, tri_neighbours[i]);

		// Make sure that the mesh is closed
		if (3 != tri_neighbours[i].size())
		{
			cout << "Error: Mesh is not closed." << endl;
			return 3;
		}

		tri_normals[i] = mesh.get_tri_normal(i);
	}

	float final_measure = 0;

	// For normalizing the measure
	// Thanks to JoeJ on gamedev.net for the idea
	const float largest_area = mesh.get_largest_triangle_area();

	// For each triangle in the mesh
	for (size_t i = 0; i < mesh.triangles.size(); i++)
	{
		// Get current triangle's face normal
		vertex_3 n_i = tri_normals[i];

		// Get neighbouring triangles' face normals
		// Assume that there are three neighbouring triangles
		// This means that the mesh must be closed
		// (e.g. no holes or cracks).
		vertex_3 o_1 = tri_normals[tri_neighbours[i][0]];
		vertex_3 o_2 = tri_normals[tri_neighbours[i][1]];
		vertex_3 o_3 = tri_normals[tri_neighbours[i][2]];

		// Average the dot products
		float d_i = (n_i.dot(o_1) + n_i.dot(o_2) + n_i.dot(o_3)) / 3.0f;

		// Normalize the average dot product
		float measure = (1.0f - d_i) / 2.0f;

		// Get current triangle area
		const float triangle_area = mesh.get_triangle_area(i);

		// Normalize the measure by area
		final_measure += measure * (triangle_area / largest_area);
	}

	// Average the measure
	float x = final_measure / mesh.triangles.size();

	// Print the dimension
	cout << "Dim: " << 2.0 + x << endl;

	return 0;
}





