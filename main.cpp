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

	float sum = 0;

	// For normalizing the measure
	// Thanks to JoeJ on gamedev.net for the idea
	const float largest_area = mesh.get_largest_triangle_area();

	double largest_len = 0;

	// For each triangle in the mesh
	for (size_t i = 0; i < mesh.triangles.size(); i++)
	{
		// Get current triangle's face normal
		vertex_3 n_i = tri_normals[i];

		// Get neighbouring triangles' face normals
		vertex_3 o_1 = tri_normals[tri_neighbours[i][0]];
		vertex_3 o_2 = tri_normals[tri_neighbours[i][1]];
		vertex_3 o_3 = tri_normals[tri_neighbours[i][2]];

		// Average the dot products
		float d_i = (n_i.dot(o_1) + n_i.dot(o_2) + n_i.dot(o_3)) / 3.0f;

		// Normalize the average dot product
		float m_i = (1.0f - d_i) / 2.0f;

		// Get current triangle area
		const float triangle_area = mesh.get_triangle_area(i);

		// Normalize the measure by area
		sum += m_i * (triangle_area / largest_area);

		// Find longest edge
		vertex_3 a = mesh.vertices[mesh.triangles[i].vertex_indices[1]] - mesh.vertices[mesh.triangles[i].vertex_indices[0]];
		vertex_3 b = mesh.vertices[mesh.triangles[i].vertex_indices[2]] - mesh.vertices[mesh.triangles[i].vertex_indices[0]];
		vertex_3 c = mesh.vertices[mesh.triangles[i].vertex_indices[2]] - mesh.vertices[mesh.triangles[i].vertex_indices[1]];

		float a_len = a.length();
		float b_len = b.length();
		float c_len = c.length();

		if (a_len > largest_len)
			largest_len = a_len;

		if (b_len > largest_len)
			largest_len = b_len;

		if (c_len > largest_len)
			largest_len = c_len;
	}

	cout << "Longest edge: " << largest_len << endl;

	// Average the measure
	float lambda = sum / mesh.triangles.size();

	// Print the dimension
	cout << "Dim: " << 2.0 + lambda << endl;

	return 0;
}





