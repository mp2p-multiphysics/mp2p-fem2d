#ifndef MESH_PHYSICSGROUP
#define MESH_PHYSICSGROUP
#include <vector>
#include "mesh_quad4.hpp"

class MeshPhysicsGroup
{
    /*

    Groups mesh domains that are used in the same physics.

    Variables
    =========
    mesh_ptr_vec_in : vector<MeshQuad4*>
        vector with pointers to MeshQuad4 objects.
    
    */

    public:

    // vector with meshes in group
    std::vector<MeshQuad4Struct*> mesh_ptr_vec;

    // default constructor
    MeshPhysicsGroup()
    {

    }

    // constructor
    MeshPhysicsGroup(std::vector<MeshQuad4Struct*> mesh_ptr_vec_in)
    {
        
        // store variables
        mesh_ptr_vec = mesh_ptr_vec_in;

    }

};

#endif
