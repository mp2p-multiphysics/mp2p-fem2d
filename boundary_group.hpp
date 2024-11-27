#ifndef BOUNDARY_GROUP
#define BOUNDARY_GROUP
#include <vector>
#include "boundary_quad4.hpp"
#include "boundary_tri3.hpp"

namespace FEM2D
{

class BoundaryGroup
{
    /*

    Groups boundary conditions (BC) that are used in the same physics.

    Variables
    =========
    boundary_l2_ptr_vec_in : vector<BoundaryLine2*>
        vector with pointers to BoundaryLine2 objects.
    
    Functions
    =========
    set_boundary_type : void
        Classifies BCs as essential or natural.
    update_parameter : void
        Recalculates non-constant BC parameters.

    */

    public:

    // vector with boundaries in group
    std::vector<BoundaryTri3*> boundary_t3_ptr_vec;
    std::vector<BoundaryQuad4*> boundary_q4_ptr_vec;

    // functions
    void set_boundary_type(VectorInt boundarytype_essential_vec, VectorInt boundarytype_natural_vec);
    void update_parameter();

    // default constructor
    BoundaryGroup() {}

    // constructor
    BoundaryGroup(std::vector<BoundaryTri3*> boundary_t3_ptr_vec_in, std::vector<BoundaryQuad4*> boundary_q4_ptr_vec_in)
    {
        boundary_t3_ptr_vec = boundary_t3_ptr_vec_in;
        boundary_q4_ptr_vec = boundary_q4_ptr_vec_in;
    }

};


void BoundaryGroup::set_boundary_type(VectorInt boundarytype_essential_vec, VectorInt boundarytype_natural_vec)
{
    /*

    Classifies BCs as essential or natural.

    Arguments
    =========
    boundarytype_essential_vec : VectorInt
        vector with BC type IDs that denote essential BCs.
    boundarytype_essential_vec : VectorInt
        vector with BC type IDs that denote natural BCs.

    Returns
    =======
    (none)

    */

    // iterate through each boundary
    for (auto boundary_ptr : boundary_t3_ptr_vec)
    {
        boundary_ptr->set_boundary_type(boundarytype_essential_vec, boundarytype_natural_vec);
    }
    for (auto boundary_ptr : boundary_q4_ptr_vec)
    {
        boundary_ptr->set_boundary_type(boundarytype_essential_vec, boundarytype_natural_vec);
    }

}

void BoundaryGroup::update_parameter()
{
    /*

    Recalculates non-constant BC parameters.

    Arguments
    =========
    (none)
    
    Returns
    =========
    (none)

    */

    // iterate through each boundary
    for (auto boundary_ptr : boundary_t3_ptr_vec)
    {
        boundary_ptr->update_parameter();
    }
    for (auto boundary_ptr : boundary_q4_ptr_vec)
    {
        boundary_ptr->update_parameter();
    }

}

}

#endif
