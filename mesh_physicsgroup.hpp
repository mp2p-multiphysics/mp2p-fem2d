#ifndef MESH_PHYSICSGROUP
#define MESH_PHYSICSGROUP
#include <vector>
#include "mesh_quad4.hpp"
#include "mesh_tri3.hpp"

class MeshPhysicsGroup
{
    /*

    Groups mesh domains that are used in the same physics.

    Variables
    =========
    mesh_t3_ptr_vec_in : vector<MeshTri3*>
        vector with pointers to MeshTri3 objects.
    mesh_q4_ptr_vec_in : vector<MeshQuad4*>
        vector with pointers to MeshQuad4 objects.

    */

    public:

    // vector with meshes in group
    std::vector<MeshTri3*> mesh_t3_ptr_vec;
    std::vector<MeshQuad4*> mesh_q4_ptr_vec;

    // default constructor
    MeshPhysicsGroup()
    {

    }

    // constructor
    MeshPhysicsGroup(std::vector<MeshTri3*> mesh_t3_ptr_vec_in, std::vector<MeshQuad4*> mesh_q4_ptr_vec_in)
    {
        
        // store variables
        mesh_t3_ptr_vec = mesh_t3_ptr_vec_in;
        mesh_q4_ptr_vec = mesh_q4_ptr_vec_in;

    }

};

#endif
