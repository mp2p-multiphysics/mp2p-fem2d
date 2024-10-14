#ifndef SCALAR_FIELDGROUP
#define SCALAR_FIELDGROUP
#include <set>
#include <unordered_map>
#include <vector>
#include "mesh_quad4.hpp"
#include "scalar_quad4.hpp"

class ScalarFieldGroup
{
    /*

    Groups scalars that are applied to the same field.

    Variables
    =========
    scalar_ptr_vec_in : vector<ScalarQuad4*>
        vector with pointers to ScalarQuad4 objects.
    
    */

    public:

    // number of unique points in field
    int num_point_field = 0;

    // point IDs
    VectorInt point_gid_vec;  // key: field ID; value: global ID
    MapIntInt point_gid_to_fid_map;  // key: global ID; value: field ID

    // scalars and meshes
    std::vector<ScalarQuad4*> scalar_ptr_vec;  // vector of scalars
    std::unordered_map<MeshQuad4Struct*, ScalarQuad4*> scalar_ptr_map;  // key: mesh; value: scalar

    // default constructor
    ScalarFieldGroup()
    {

    }

    // constructor
    ScalarFieldGroup(std::vector<ScalarQuad4*> scalar_ptr_vec_in)
    {
        
        // store vector of scalars
        scalar_ptr_vec = scalar_ptr_vec_in;

        // map mesh to scalars
        for (auto scalar_ptr : scalar_ptr_vec)
        {
            scalar_ptr_map[scalar_ptr->mesh_q4_ptr] = scalar_ptr;
        }

        // get set of global IDs
        // map global IDs and field IDs

        // initialize set of global IDs
        std::set<int> point_gid_set;  

        // iterate through each variable and get set of global IDs
        for (auto scalar_ptr : scalar_ptr_vec)
        {
            for (auto &point_gid : scalar_ptr->mesh_q4_ptr->point_gid_vec)
            {
                point_gid_set.insert(point_gid);
            }
        }

        // initialize field ID
        int point_fid = 0;

        // iterate through each global ID and assign a field ID
        for (auto point_gid : point_gid_set)
        {

            // skip if global ID is already recorded
            if (point_gid_to_fid_map.count(point_gid))
            {
                continue;
            }

            // map global ID to field ID and vice versa
            point_gid_to_fid_map[point_gid] = point_fid;
            point_gid_vec.push_back(point_gid);
            point_fid++;

        }

        // total number of field points
        num_point_field = point_fid;

    }

};

#endif