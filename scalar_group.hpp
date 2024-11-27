#ifndef SCALAR_GROUP
#define SCALAR_GROUP
#include <set>
#include <unordered_map>
#include <vector>
#include "container_typedef.hpp"
#include "domain_quad4.hpp"
#include "domain_tri3.hpp"
#include "scalar_quad4.hpp"
#include "scalar_tri3.hpp"

namespace FEM2D
{

class ScalarGroup
{
    /*

    Groups scalars that are applied to the same group.

    Variables
    =========
    scalar_t3_ptr_vec_in : vector<ScalarTri3*>
        vector with pointers to ScalarTri3 objects.
    scalar_q4_ptr_vec_in : vector<ScalarQuad4*>
        vector with pointers to ScalarQuad4 objects.

    Functions
    =========
    update_value : void
        Recalculates non-constant values.

    */

    public:

    // number of unique points in group
    int num_point_group = 0;

    // point IDs
    VectorInt point_pfid_to_pgid_vec;  // key: group ID; value: global ID
    MapIntInt point_pgid_to_pfid_map;  // key: global ID; value: group ID

    // scalars and domains
    std::vector<ScalarTri3*> scalar_t3_ptr_vec;  // vector of scalars
    std::vector<ScalarQuad4*> scalar_q4_ptr_vec;  // vector of scalars

    // functions
    void update_value();

    // default constructor
    ScalarGroup() {}

    // constructor
    ScalarGroup(std::vector<ScalarTri3*> scalar_t3_ptr_vec_in, std::vector<ScalarQuad4*> scalar_q4_ptr_vec_in)
    {
        
        // store vector of scalars
        scalar_t3_ptr_vec = scalar_t3_ptr_vec_in;
        scalar_q4_ptr_vec = scalar_q4_ptr_vec_in;

        // get set of global IDs
        // map global IDs and group IDs

        // initialize set of global IDs
        std::set<int> point_pgid_set;  

        // iterate through each variable and get set of global IDs
        for (auto scalar_ptr : scalar_t3_ptr_vec)
        {
            for (auto &pgid : scalar_ptr->domain_ptr->point_pdid_to_pgid_vec)
            {
                point_pgid_set.insert(pgid);
            }
        }
        for (auto scalar_ptr : scalar_q4_ptr_vec)
        {
            for (auto &pgid : scalar_ptr->domain_ptr->point_pdid_to_pgid_vec)
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
            pfid++;

        }

        // total number of group points
        num_point_group = pfid;

    }

};

void ScalarGroup::update_value()
{
    /*

    Recalculates non-constant values.

    Arguments
    =========
    (none)
    
    Returns
    =========
    (none)

    */

    // iterate through each scalar
    for (auto scalar_ptr : scalar_t3_ptr_vec)
    {
        scalar_ptr->update_value();
    }
    for (auto scalar_ptr : scalar_q4_ptr_vec)
    {
        scalar_ptr->update_value();
    }

}

}

#endif
