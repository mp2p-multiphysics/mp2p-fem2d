#ifndef SCALAR_FIELDGROUP
#define SCALAR_FIELDGROUP
#include <set>
#include <unordered_map>
#include <vector>
#include "mesh_quad4.hpp"
#include "mesh_tri3.hpp"
#include "scalar_quad4.hpp"
#include "scalar_tri3.hpp"

class ScalarFieldGroup
{
    /*

    Groups scalars that are applied to the same field.

    Variables
    =========
    scalar_ptr_vec_in : vector<ScalarBase*>
        vector with pointers to ScalarBase objects.
    
    */

    public:

    // number of unique points in field
    int num_point_field = 0;

    // point IDs
    VectorInt point_gid_vec;  // key: field ID; value: global ID
    MapIntInt point_gid_to_fid_map;  // key: global ID; value: field ID

    // scalars and meshes
    std::vector<ScalarTri3*> scalar_t3_ptr_vec;  // vector of scalars
    std::vector<ScalarQuad4*> scalar_q4_ptr_vec;  // vector of scalars
    std::unordered_map<MeshTri3*, ScalarTri3*> mesh_to_scalar_t3_ptr_map;  // key: mesh; value: scalar
    std::unordered_map<MeshQuad4*, ScalarQuad4*> mesh_to_scalar_q4_ptr_map;  // key: mesh; value: scalar

    // default constructor
    ScalarFieldGroup()
    {

    }

    // constructor
    ScalarFieldGroup(std::vector<ScalarTri3*> scalar_t3_ptr_vec_in, std::vector<ScalarQuad4*> scalar_q4_ptr_vec_in)
    {
        
        // store vector of scalars
        scalar_t3_ptr_vec = scalar_t3_ptr_vec_in;
        scalar_q4_ptr_vec = scalar_q4_ptr_vec_in;

        // map mesh to scalars
        for (auto scalar_ptr : scalar_t3_ptr_vec)
        {
            mesh_to_scalar_t3_ptr_map[scalar_ptr->mesh_t3_ptr] = scalar_ptr;
        }
        for (auto scalar_ptr : scalar_q4_ptr_vec)
        {
            mesh_to_scalar_q4_ptr_map[scalar_ptr->mesh_q4_ptr] = scalar_ptr;
        }

        // get set of global IDs
        // map global IDs and field IDs

        // initialize set of global IDs
        std::set<int> point_gid_set;  

        // iterate through each variable and get set of global IDs
        for (auto scalar_ptr : scalar_t3_ptr_vec)
        {
            for (auto &point_gid : scalar_ptr->mesh_t3_ptr->point_gid_vec)
            {
                point_gid_set.insert(point_gid);
            }
        }
        for (auto scalar_ptr : scalar_q4_ptr_vec)
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
