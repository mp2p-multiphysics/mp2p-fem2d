#ifndef INTEGRAL_GROUP
#define INTEGRAL_GROUP
#include <vector>
#include "integral_unit.hpp"

namespace FEM2D
{

class IntegralGroup
{
    /*

    Groups test function (N) integrals that are used in the same physics.

    Variables
    =========
    integral_ptr_vec_in : vector<IntegralTri3*>
        vector with pointers to IntegralTri3 objects.
    integral_q4_ptr_vec_in : vector<IntegralQuad4*>
        vector with pointers to IntegralQuad4 objects.

    Functions
    =========
    evaluate_Ni : void
        Calculates test function values and other properties.
        Must be called before domain integrals are evaluated.
    evaluate_integral_Ni : void
        Calculates the integral of Ni.
    evaluate_integral_derivative_Ni_x : void
        Calculates the integral of d(Ni)/dx.
    evaluate_integral_derivative_Ni_y : void
        Calculates the integral of d(Ni)/dy.
    evaluate_integral_Ni_Nj : void
        Calculates the integral of Ni * Nj.
    evaluate_integral_Ni_derivative_Nj_x : void
        Calculates the integral of Ni * d(Nj)/dx.
    evaluate_integral_Ni_derivative_Nj_y : void
        Calculates the integral of Ni * d(Nj)/dy.
    evaluate_integral_div_Ni_dot_div_Nj : void
        Calculates the integral of div(Ni) dot div(Nj).
    evaluate_integral_Ni_Nj_derivative_Nk_x : void
        Calculates the integral of Ni * Nj * d(Nk)/dx.
    evaluate_integral_Ni_Nj_derivative_Nk_y : void
        Calculates the integral of Ni * Nj * d(Nk)/dy.
    evaluate_boundary_Ni : void
        Calculates test functions values and other properties at the boundary.
        Must be called before boundary integrals are evaluated.
    evaluate_integral_boundary_Ni : void
        Calculates the integral of Ni at the boundary.
    evaluate_integral_boundary_Ni_Nj : void
        Calculates the integral of Ni * Nj at the boundary.

    */

    public:

    // variables
    std::vector<IntegralUnit*> integral_ptr_vec;

    // functions for computing domain integrals
    void evaluate_integral_Ni();
    void evaluate_integral_derivative_Ni_x();
    void evaluate_integral_derivative_Ni_y();
    void evaluate_integral_Ni_Nj();
    void evaluate_integral_Ni_derivative_Nj_x();
    void evaluate_integral_Ni_derivative_Nj_y();
    void evaluate_integral_div_Ni_dot_div_Nj();
    void evaluate_integral_Ni_Nj_derivative_Nk_x();
    void evaluate_integral_Ni_Nj_derivative_Nk_y();

    // functions for computing boundary integrals
    void evaluate_integral_boundary_Ni();
    void evaluate_integral_boundary_Ni_Nj();

    // default constructor
    IntegralGroup() {}

    // constructor
    IntegralGroup(std::vector<IntegralUnit*> integral_ptr_vec_in)
    {
        integral_ptr_vec = integral_ptr_vec_in;
    }

};

void IntegralGroup::evaluate_integral_Ni()
{
    /*

    Calculates the integral of Ni.

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

void IntegralGroup::evaluate_integral_derivative_Ni_x()
{
    /*

    Calculates the integral of d(Ni)/dx.

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
        integral_ptr->evaluate_integral_derivative_Ni_x();
    }

}

void IntegralGroup::evaluate_integral_derivative_Ni_y()
{
    /*

    Calculates the integral of d(Ni)/dy.

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
        integral_ptr->evaluate_integral_derivative_Ni_y();
    }

}

void IntegralGroup::evaluate_integral_Ni_Nj()
{
    /*

    Calculates the integral of Ni * Nj.

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

void IntegralGroup::evaluate_integral_Ni_derivative_Nj_x()
{
    /*

    Calculates the integral of Ni * d(Nj)/dx.

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
        integral_ptr->evaluate_integral_Ni_derivative_Nj_x();
    }

}

void IntegralGroup::evaluate_integral_Ni_derivative_Nj_y()
{
    /*

    Calculates the integral of Ni * d(Nj)/dy.

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
        integral_ptr->evaluate_integral_Ni_derivative_Nj_y();
    }

}

void IntegralGroup::evaluate_integral_div_Ni_dot_div_Nj()
{
    /*

    Calculates the integral of div(Ni) dot div(Nj).

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
        integral_ptr->evaluate_integral_div_Ni_dot_div_Nj();
    }

}

void IntegralGroup::evaluate_integral_Ni_Nj_derivative_Nk_x()
{
    /*

    Calculates the integral of Ni * Nj * d(Nk)/dx.

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
        integral_ptr->evaluate_integral_Ni_Nj_derivative_Nk_x();
    }

}

void IntegralGroup::evaluate_integral_Ni_Nj_derivative_Nk_y()
{
    /*

    Calculates the integral of Ni * Nj * d(Nk)/dy.

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
        integral_ptr->evaluate_integral_Ni_Nj_derivative_Nk_y();
    }

}

void IntegralGroup::evaluate_integral_boundary_Ni()
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
        integral_ptr->evaluate_integral_boundary_Ni();
    }

}

void IntegralGroup::evaluate_integral_boundary_Ni_Nj()
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
        integral_ptr->evaluate_integral_boundary_Ni_Nj();
    }

}

}

#endif
