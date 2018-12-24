#include "string_utils.h"

bool is_real_number(const string &src_string)
{
	//ie: 
	//1
	//-23e4
	//1.E2
	//-2.717
	//.31415e1
	//-7.53e-9
	//7.53e+9

	if(src_string == "")
		return false;

	string temp = lower_string(src_string);

	bool found_dot = false;
	bool found_e = false;
	bool found_digit = false;

	for(size_t i = 0; i < temp.size(); i++)
	{
		if(isdigit(temp[i]))
		{
			if(found_digit == false)
				found_digit = true;	
		}
		else if(temp[i] == 'e')
		{
			if(found_e == true || found_digit == false || i == temp.size() - 1)
				return false;
			else
				found_e = true;
		}
		else if(temp[i] == '-' || temp[i] == '+')
		{
			if(!(i == 0 || temp[i-1] == 'e') || i == temp.size() - 1)
				return false;
		}
		else if(temp[i] == '.')
		{
			if(found_dot == true || (i != 0 && temp[i-1] == 'e'))
				return false;
			else
				found_dot = true;
		}
		else
		{
			return false;
		}
	}

	return true;
}

string lower_string(const string &src_string)
{
	string temp = src_string;

	for(string::iterator i = temp.begin(); i != temp.end(); i++)
		*i = tolower(*i);

	return temp;
}
