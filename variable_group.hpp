#ifndef VARIABLE_GROUP
#define VARIABLE_GROUP
#include <set>
#include <unordered_map>
#include <vector>
#include "container_typedef.hpp"
#include "variable_quad4.hpp"
#include "variable_tri3.hpp"

namespace FEM2D
{

class VariableGroup
{
    /*

    Groups variables that are applied to the same group.

    Variables
    =========
    variable_l2_ptr_vec_in : vector<VariableLine2*>
        vector with pointers to VariableLine2 objects.
    
    */

    public:

    // number of unique points in group
    int num_point = 0;

    // point IDs
    VectorInt point_pfid_to_pgid_vec;  // key: group ID; value: global ID
    MapIntInt point_pgid_to_pfid_map;  // key: global ID; value: group ID

    // variables and domains
    std::vector<VariableTri3*> variable_t3_ptr_vec;  // vector of variables
    std::vector<VariableQuad4*> variable_q4_ptr_vec;  // vector of variables

    // starting column of variables in matrix equation
    int start_col = -1;

    // default constructor
    VariableGroup() {}

    // constructor
    VariableGroup(std::vector<VariableTri3*> variable_t3_ptr_vec_in, std::vector<VariableQuad4*> variable_q4_ptr_vec_in)
    {
        
        // store vector of variables
        variable_t3_ptr_vec = variable_t3_ptr_vec_in;
        variable_q4_ptr_vec = variable_q4_ptr_vec_in;

        // get set of global IDs
        // map global IDs and group IDs

        // initialize set of global IDs
        std::set<int> point_pgid_set;  

        // iterate through each variable and get set of global IDs
        for (auto variable_ptr : variable_t3_ptr_vec)
        {
            for (auto &pgid : variable_ptr->domain_ptr->point_pdid_to_pgid_vec)
            {
                point_pgid_set.insert(pgid);
            }
        }
        for (auto variable_ptr : variable_q4_ptr_vec)
        {
            for (auto &pgid : variable_ptr->domain_ptr->point_pdid_to_pgid_vec)
            {
                point_pgid_set.insert(pgid);
            }
        }

        // initialize group ID
        int pfid = 0;

        // iterate through each global ID and assign a group ID
        for (auto pgid : point_pgid_set)
        {

            // skip if global ID is already recorded
            if (point_pgid_to_pfid_map.count(pgid))
            {
                continue;
            }

            // map global ID to group ID and vice versa
            point_pgid_to_pfid_map[pgid] = pfid;
            point_pfid_to_pgid_vec.push_back(pgid);
            
            // increment group ID
            pfid++;

        }

        // total number of group points
        num_point = pfid;

    }

};

}

#endif
