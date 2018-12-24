#include "main.h"
#include "mesh.h"

int main(int argc, char **argv)
{
	if(argc != 2)
	{
		cout << "Example usage: meshdim filename.stl" << endl;
		return 1;
	}

	indexed_mesh mesh;

	if(false == mesh.load_from_binary_stereo_lithography_file(argv[1]))
	{
		cout << "Error: Could not properly read file " << argv[1] << endl;
		return 2;
	}

	vector< vector<size_t> > tri_neighbours;
	vector<vertex_3> tri_normals;

	tri_neighbours.resize(mesh.triangles.size());
	tri_normals.resize(mesh.triangles.size());

	for(size_t i = 0; i < mesh.triangles.size(); i++)
	{
		mesh.get_tri_neighbours(i, tri_neighbours[i]);
		tri_normals[i] = mesh.get_tri_normal(i);
	}

	float final_measure = 0;

	for (size_t i = 0; i < mesh.triangles.size(); i++)
	{
		vertex_3 n0 = tri_normals[i];
		vertex_3 n1 = tri_normals[tri_neighbours[i][0]];
		vertex_3 n2 = tri_normals[tri_neighbours[i][1]];
		vertex_3 n3 = tri_normals[tri_neighbours[i][2]];

		float dot1 = n0.dot(n1);
		float dot2 = n0.dot(n2);
		float dot3 = n0.dot(n3);

		float d = (dot1 + dot2 + dot3) / 3.0f;
		float measure = (1.0f - d) / 2.0f;

		final_measure += measure;
	}

	cout << "Dim: " << 2.0 + final_measure/mesh.triangles.size() << endl;

	return 0;
}





