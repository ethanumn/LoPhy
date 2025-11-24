#include "utils.h"

#include <numeric>   
#include <string>    

using std::reduce;
using std::prev_permutation;
using std::string;

/* Combinations of elements in array */
void comb(vector<vector<int>> & choices, 
          const vector<int> & arr, 
          const int & k)
{
    int n = arr.size();
    string bitmask(k, 1); 
    bitmask.resize(n, 0);
 
    // permute bitmask and use this to select k elements from arr
    do {
        vector<int> c;
        for (int i = 0; i < n; ++i) 
        {
            if (bitmask[i]) 
                c.push_back(arr[i]);
        }
        choices.push_back(c);
    } while (prev_permutation(bitmask.begin(), bitmask.end()));
}
 
/* Get all descendants of a node in a parents vector */
 void get_descendants(unordered_set<int> & descendants, 
                     const vector<int> & parents, 
                     const int & node)
{
    int n = parents.size();
    descendants.clear();
    vector<int> Q = {node};
    while(Q.size() > 0)
    {
        const int d = Q.back();
        Q.pop_back();

        // record descendant
        if(d != node)
            descendants.insert(d);

        for(int i = 0; i < n; ++i)
        {
            if(parents[i] == NO_PARENT)
                continue;
            else if((int)parents[i] == d)
                Q.push_back(i+1);
        }
    }
}

/* Get all ancestors of a node in a parents vector */
void get_ancestors(unordered_set<int> & ancestors, 
                   const vector<int> & parents, 
                   const int & node)
{
    ancestors.clear();
    if(node != ROOT)
    {
        auto a = parents[node - 1];
        while(a != ROOT)
        {
            ancestors.insert(a);
            a = parents[a - 1];
        }
    }
}