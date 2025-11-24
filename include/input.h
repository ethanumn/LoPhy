#include <vector> 
#include <string>
#include <fstream>
#include <iostream>

using std::vector;
using std::string;
using std::ifstream;


/* Reads a file where each line has a single value on it */
template<typename T>
int read_lines(const string & file_name, 
                const uint32_t & num_lines,
                vector<T> & lines)
{

	ifstream file(file_name);

    // if file doesn't exist then return -1
	if (!file) {
	    return -1;
	} 
    else 
    {
        // extract each of the gene names from the file
        T temp;
        for (size_t i = 0; i < num_lines; ++i) {
            file >> temp;
            lines.push_back(temp);
        }
        return 0;
    }
}