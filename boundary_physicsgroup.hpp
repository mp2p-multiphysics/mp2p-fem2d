#ifndef BOUNDARY_PHYSICSGROUP
#define BOUNDARY_PHYSICSGROUP
#include <vector>
#include "boundary_quad4.hpp"
#include "boundary_tri3.hpp"

class BoundaryPhysicsGroup
{
    /*

    Groups boundary conditions (BC) that are used in the same physics.

    Variables
    =========
    boundary_t3_ptr_vec_in : vector<BoundaryTri3*>
        vector with pointers to BoundaryTri3 objects.
    boundary_q4_ptr_vec_in : vector<BoundaryQuad4*>
        vector with pointers to BoundaryQuad4 objects.

    */

    public:

    // vector with boundaries in group
    std::vector<BoundaryTri3*> boundary_t3_ptr_vec;
    std::vector<BoundaryQuad4*> boundary_q4_ptr_vec;

    // default constructor
    BoundaryPhysicsGroup()
    {

    }

    // constructor
    BoundaryPhysicsGroup(std::vector<BoundaryTri3*> boundary_t3_ptr_vec_in, std::vector<BoundaryQuad4*> boundary_q4_ptr_vec_in)
    {
        boundary_t3_ptr_vec = boundary_t3_ptr_vec_in;
        boundary_q4_ptr_vec = boundary_q4_ptr_vec_in;
    }

};

#endif
