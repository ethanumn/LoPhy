#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <algorithm>

using namespace std;

/* DataMatrix class used to store/process comma separated data */
template<typename T>
class DataMatrix {
    public:
        void set_columns(const vector<string> & columns) 
        {
            this->columns = columns;
            column_map.clear();
            for (int i = 0; i < static_cast<int>(columns.size()); ++i) column_map[columns[i]] = i;
        }
        void set_index(vector<string> index) {this->index = index;}
        void add_row(vector<T> row) {this->data.push_back(row);}
        vector<string> get_columns(void) {return this->columns;}
        vector<T> get_column(const string & col) const 
        {
            if (column_map.find(col) == column_map.end()) {throw invalid_argument("Column name does not exist.");}
            int col_idx = column_map.at(col);
            vector<T> column_data;
            for (const auto & row : data) column_data.push_back(row[col_idx]);
            return column_data;
        }
        vector<string> get_index(void) {return this->index;}
        vector<vector<T>> get_data(void) {return this->data;}
        void clear(void) {this->columns.clear(); this->index.clear(); this->data.clear();}
        void print(void) {
            auto n = this->index.size();
            for(auto col : columns)
                cout << col << " ";
            cout << "\n";
            for(auto i = 0; i < n; ++i)
            {
                cout << this->index[i] << " ";

                for(auto v : this->data[i])
                    cout << v << " ";
                cout << "\n";
            }
        }

    private:
        vector<string> columns;
        vector<string> index;
        vector<vector<T>> data;
        unordered_map<string, size_t> column_map;
};

/* Extracts region names from input matrix */
inline vector<string> extract_region_names(const vector<string> & string_vec)
{
    vector<string> region_names;
    for (const auto & val : string_vec) 
    {
        // Find position of first underscore
        size_t underscore_pos = val.find('_');
        if (underscore_pos == string::npos || underscore_pos + 1 >= val.size()) {
            region_names.push_back(""); // or skip or push val itself
            continue;
        }

        // Extract substring after underscore
        string name = val.substr(underscore_pos + 1);
        
        // Remove "Region" or "region" prefix if present
        const string prefix1 = "Region";
        const string prefix2 = "region";

        if (name.rfind(prefix1, 0) == 0) 
            name.erase(0, prefix1.length());
        else if (name.rfind(prefix2, 0) == 0) 
            name.erase(0, prefix2.length());

        region_names.push_back(name);
    }
    
    return region_names;
}

/* Extracts chromosome number from input matrix */
inline vector<int> extract_chromosomes(const vector<string> & string_vec)
{
    vector<int> chromosomes;
    for (const auto & val : string_vec) 
    {
        try {
            string chr = val.substr(0, val.find("_"));
            transform(chr.begin(), chr.end(), chr.begin(), ::toupper);
            if(chr == "X")
                chromosomes.push_back(23);
            else if(chr == "Y")
                chromosomes.push_back(24);
            else
                chromosomes.push_back(static_cast<int>(stoul(chr)));
        } catch (const invalid_argument& e) {
            throw invalid_argument("Invalid chromosome value to convert to int: " + val);
        } catch (const out_of_range& e) {
            throw out_of_range("Value out of range for converting to int: " + val);
        }
    }
    
    return chromosomes;
}

/* Convert from a vector<string> to vector<int> */
inline vector<int> convert_to_uint32(const vector<string> & string_vec) 
{
    vector<int> int_vec;
    for (const auto & val : string_vec) {
        try {
            int_vec.push_back(static_cast<int>(stoul(val)));
        } catch (const invalid_argument& e) {
            throw invalid_argument("Invalid value in column for int conversion: " + val);
        } catch (const out_of_range& e) {
            throw out_of_range("Value out of range for converting to int: " + val);
        }
    }
    return int_vec;
}

/* Convert from vector<string> to vector<bool>*/
inline vector<bool> convert_to_bool(const vector<string> & string_vec) 
{
    vector<bool> bool_vec;
    for (const auto & val : string_vec) {
        if (val == "1" || val == "true" || val == "True") {
            bool_vec.push_back(true);
        } else if (val == "0" || val == "false" || val == "False") {
            bool_vec.push_back(false);
        } else {
            throw invalid_argument("Invalid value in column: " + val);
        }
    }
    return bool_vec;
}

/* reads csv file and places context in the DataMatrix */
template<typename T>
void read_csv(string file_name, 
              DataMatrix<T> & data_matrix, 
              int index_col = 0, 
              int header = 0)
{

    // base case
    if(file_name == "")
    {
        return; 
    }

    ifstream infile(file_name);

    // we'll fill these and copy them into the DataMatrix
    vector<string> index;
    vector<string> columns;

    // keep track of which row/column we're at
    auto row = 0;
    auto col = 0;

    while(infile)
    {
        string s; 
        if (!getline(infile, s)) break;
    
        istringstream ss(s);
        vector<T> record;
        col = 0; 

        while (ss)
        {
            string s;
            if (!getline(ss, s, ',')) break;
            if(row == header)
                if(index_col == col) {} // do nothing
                else columns.push_back(s);
            else
            {
                if(col == index_col)
                {
                    if(row == header) {} // do nothing
                    else 
                        index.push_back(s);
                }
                else 
                {
                    T value;
                    std::istringstream iss(s);
                    iss >> value;
                    record.push_back(value);
                }
            }
            col +=1;
        }

        if(row == header) {} // do nothing 
        else data_matrix.add_row(record);
        row += 1;


    }
    if (!infile.eof())
    {
        cerr << "There may be an issue with the input file!\n";
    }

    // set all elements of the data matrix
    data_matrix.set_columns(columns);
    data_matrix.set_index(index);
}

/* write csv where input is a DataMatrix */
template<typename T>
void write_csv(string file_name, 
              DataMatrix<T> & data_matrix, 
              vector<string> index = {}, 
              vector<string> header = {},
              string sep = ",")
{
    write_data(file_name, 
               data_matrix.get_data(), 
               index, 
               header,
               sep);

}

/* write csv where input is a vector of vectors */
template<typename T>
void write_csv(string file_name, 
              vector<vector<T>> & data, 
              vector<string> index = {}, 
              vector<string> header = {},
              string sep = ",")
{
    write_data(file_name, 
               data, 
               index, 
               header,
               sep);

}

template<typename T>
void write_data(string file_name, 
                vector<vector<T>> & data, 
                vector<string> index, 
                vector<string> header,
                string sep)
{
    // Open the file in write mode
    std::ofstream file(file_name);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << file_name << std::endl;
        return;
    }

    // Write the header if specified and it is not empty
    if (!header.empty()) {
        if (!index.empty()) {
            file << "Index" << sep;  // Add "Index" column header if the index is present
        }
        for (size_t i = 0; i < header.size(); ++i) {
            file << header[i];
            if (i < header.size() - 1) file << sep;
        }
        file << "\n";  // End the header row
    }

    // Write the data
    for (int row = 0; row < static_cast<int>(data.size()); ++row) {
        // Write index if specified
        if (!index.empty()) {
            file << index[row] << sep; // Writing the custom index
        }

        // Write the row data
        for (int col = 0; col < static_cast<int>(data[row].size()); ++col) {
            file << data[row][col];
            if (col < data[row].size() - 1) {
                file << sep; // Adding separator between columns
            }
        }

        file << "\n";  // New line after each row
    }

    // Close the file
    file.close();

}

