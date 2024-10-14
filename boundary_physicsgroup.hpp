#ifndef BOUNDARY_PHYSICSGROUP
#define BOUNDARY_PHYSICSGROUP
#include <vector>
#include "boundary_quad4.hpp"

class BoundaryPhysicsGroup
{
    /*

    Groups boundary conditions (BC) that are used in the same physics.

    Variables
    =========
    boundary_ptr_vec_in : vector<BoundaryQuad4*>
        vector with pointers to BoundaryQuad4 objects.
    
    */

    public:

    // vector with boundaries in group
    std::vector<BoundaryQuad4*> boundary_ptr_vec;

    // default constructor
    BoundaryPhysicsGroup()
    {

    }

    // constructor
    BoundaryPhysicsGroup(std::vector<BoundaryQuad4*> boundary_ptr_vec_in)
    {
        boundary_ptr_vec = boundary_ptr_vec_in;
    }

};

#endif
