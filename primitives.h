// Source code by Shawn Halayka
// Source code is in the public domain

#ifndef PRIMITIVES_H
#define PRIMITIVES_H


#include <cmath>
#include <cstddef> // g++ chokes on size_t without this

#include <vector>
using std::vector;

#include <iostream>
using namespace std;


class vertex_3
{
public:
	inline vertex_3(void) : x(0.0f), y(0.0f), z(0.0f) { /*default constructor*/ }
	inline vertex_3(const float src_x, const float src_y, const float src_z) : x(src_x), y(src_y), z(src_z) { /* custom constructor */ }
	inline vertex_3(const vertex_3 &rhs) : x(rhs.x), y(rhs.y), z(rhs.z) { /* custom constructor */ }

	inline bool operator==(const vertex_3 &right) const
	{
		if(x == right.x && y == right.y && z == right.z)
			return true;

		return false;
	}

	inline bool operator<(const vertex_3 &right) const
	{
		if(x < right.x)
			return true;
		else if(x > right.x)
			return false;

		if(y < right.y)
			return true;
		else if(y > right.y)
			return false;

		if(z < right.z)
			return true;
		else if(z > right.z)
			return false;

		return false;
	}

	inline bool operator>(const vertex_3 &right) const
	{
		if(x > right.x)
			return true;
		else if(x < right.x)
			return false;

		if(y > right.y)
			return true;
		else if(y < right.y)
			return false;

		if(z > right.z)
			return true;
		else if(z < right.z)
			return false;

		return false;
	}


	inline vertex_3& operator+=(const vertex_3 &right)
	{
		x += right.x;
		y += right.y;
		z += right.z;

		return *this;
	}

	inline vertex_3& operator*=(const float &right)
	{
		x *= right;
		y *= right;
		z *= right;

		return *this;
	}


	inline vertex_3& operator=(const vertex_3 &right)
	{
		x = right.x;
		y = right.y;
		z = right.z;

		return *this;
	}

	inline vertex_3 operator-(const vertex_3 &right) const
	{
		vertex_3 temp;

		temp.x = x - right.x;
		temp.y = y - right.y;
		temp.z = z - right.z;

		return temp;
	}

	inline vertex_3 operator+(const vertex_3 &right) const
	{
		vertex_3 temp;

		temp.x = x + right.x;
		temp.y = y + right.y;
		temp.z = z + right.z;

		return temp;
	}

	inline vertex_3 operator*(const float &right) const
	{
		vertex_3 temp;

		temp.x = x * right;
		temp.y = y * right;
		temp.z = z * right;

		return temp;
	}

	inline vertex_3 operator/(const float &right) const
	{
		vertex_3 temp;

		temp.x = x / right;
		temp.y = y / right;
		temp.z = z / right;

		return temp;
	}


	inline vertex_3 cross(const vertex_3 &right) const
	{
		vertex_3 temp;

		temp.x = y*right.z - z*right.y;
		temp.y = z*right.x - x*right.z;
		temp.z = x*right.y - y*right.x;

		return temp;
	}

	inline float dot(const vertex_3 &right) const
	{
		return x*right.x + y*right.y + z*right.z;
	}

	inline float self_dot(void) const
	{
		return x*x + y*y + z*z;
	}

	inline float length(void) const
	{
		return sqrt(self_dot());
	}

	inline float distance(const vertex_3 &right) const
	{
		return sqrt((right.x - x)*(right.x - x) + (right.y - y)*(right.y - y) + (right.z - z)*(right.z - z));
	}

	inline float distance_sq(const vertex_3 &right) const
	{
		return (right.x - x)*(right.x - x) + (right.y - y)*(right.y - y) + (right.z - z)*(right.z - z);
	}

	inline void normalize(void)
	{
		float len = length();

		if(0.0f != len)
		{
			x /= len;
			y /= len;
			z /= len;
		}
	}

	inline void zero(void)
	{
		x = y = z = 0;
	}

	inline void rotate_x(const float &radians)
	{
		float t_y = y;

		y = t_y*cos(radians) + z*sin(radians);
		z = t_y*-sin(radians) + z*cos(radians);
	}

	inline void rotate_y(const float &radians)
	{
		float t_x = x;

		x = t_x*cos(radians) + z*-sin(radians);
		z = t_x*sin(radians) + z*cos(radians);
	}

	float x, y, z;
};

class indexed_vertex_3 : public vertex_3
{
public:
	inline indexed_vertex_3(void) { x = y = z = 0; index = 0; }
	inline indexed_vertex_3(const float src_x, const float src_y, const float src_z, const size_t src_index) { x = src_x; y = src_y; z = src_z; index = src_index; }
	inline indexed_vertex_3(const float src_x, const float src_y, const float src_z) { x = src_x; y = src_y; z = src_z; index = 0; }

	inline bool operator<(const vertex_3 &right) const
	{
		if(right.x > x)
			return true;
		else if(right.x < x)
			return false;

		if(right.y > y)
			return true;
		else if(right.y < y)
			return false;

		if(right.z > z)
			return true;
		else if(right.z < z)
			return false;

		return false;
	}

	inline bool operator>(const vertex_3 &right) const
	{
		if(right.x < x)
			return true;
		else if(right.x > x)
			return false;

		if(right.y < y)
			return true;
		else if(right.y > y)
			return false;

		if(right.z < z)
			return true;
		else if(right.z > z)
			return false;

		return false;
	}

	size_t index;
};

class indexed_triangle
{
public:
	size_t vertex_indices[3];

	inline bool operator==(const indexed_triangle &right) const
	{
		if( right.vertex_indices[0] == vertex_indices[0] && 
			right.vertex_indices[1] == vertex_indices[1] &&
			right.vertex_indices[2] == vertex_indices[2] )
		{
			return true;
		}

		return false;
	}
};


class ordered_size_t_pair
{
public:
	ordered_size_t_pair(const size_t &a, const size_t &b)
	{
		if(a < b)
		{
			indices[0] = a;
			indices[1] = b;
		}
		else
		{
			indices[0] = b;
			indices[1] = a;
		}
	}

	bool operator<(const ordered_size_t_pair &right) const
	{
		if(indices[0] < right.indices[0])
			return true;
		else if(indices[0] > right.indices[0])
			return false;

		if(indices[1] < right.indices[1])
			return true;
		else if(indices[1] > right.indices[1])
			return false;

		return false;
	}

	size_t indices[2];
};

// Should probably use inheritance here, but whatever...
class ordered_indexed_edge
{
public:
	ordered_indexed_edge(const indexed_vertex_3 &a, const indexed_vertex_3 &b)
	{
		if(a.index < b.index)
		{
			indices[0] = a.index;
			indices[1] = b.index;
		}
		else
		{
			indices[0] = b.index;
			indices[1] = a.index;
		}

		centre_point.x = (a.x + b.x)*0.5f;
		centre_point.y = (a.y + b.y)*0.5f;
		centre_point.z = (a.z + b.z)*0.5f;
	}

	bool operator<(const ordered_indexed_edge &right) const
	{
		if(indices[0] < right.indices[0])
			return true;
		else if(indices[0] > right.indices[0])
			return false;

		if(indices[1] < right.indices[1])
			return true;
		else if(indices[1] > right.indices[1])
			return false;

		return false;
	}

	size_t indices[2];
	vertex_3 centre_point;
	size_t id;
};


#endif
