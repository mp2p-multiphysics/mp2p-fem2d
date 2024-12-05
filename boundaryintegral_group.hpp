#ifndef BOUNDARYINTEGRAL_GROUP
#define BOUNDARYINTEGRAL_GROUP
#include <vector>
#include "boundaryintegral_unit.hpp"

namespace FEM2D
{

class BoundaryIntegralGroup
{
    /*

    Groups test function (N) integrals that are used in the same physics.

    Variables
    =========
    integral_ptr_vec_in : vector<BoundaryIntegralUnit*>
        vector with pointers to BoundaryIntegralUnit objects.

    Functions
    =========
    evaluate_integral_Ni : void
        Calculates the integral of Ni at the boundary.
    evaluate_integral_Ni_Nj : void
        Calculates the integral of Ni * Nj at the boundary.

    */

    public:

    // variables
    std::vector<BoundaryIntegralUnit*> integral_ptr_vec;

    // functions for computing boundary integrals
    void evaluate_integral_Ni();
    void evaluate_integral_Ni_Nj();

    // default constructor
    BoundaryIntegralGroup() {}

    // constructor
    BoundaryIntegralGroup(std::vector<BoundaryIntegralUnit*> integral_ptr_vec_in)
    {
        integral_ptr_vec = integral_ptr_vec_in;
    }

};

void BoundaryIntegralGroup::evaluate_integral_Ni()
{
    /*

    Calculates the integral of Ni at the boundary.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni();
    }

}

void BoundaryIntegralGroup::evaluate_integral_Ni_Nj()
{
    /*

    Calculates the integral of Ni * Nj at the boundary.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_Nj();
    }

}

}

#endif
