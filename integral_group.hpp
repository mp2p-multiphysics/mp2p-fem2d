#ifndef INTEGRAL_GROUP
#define INTEGRAL_GROUP
#include <vector>
#include "integral_quad4.hpp"
#include "integral_tri3.hpp"

namespace FEM2D
{

class IntegralGroup
{
    /*

    Groups test function (N) integrals that are used in the same physics.

    Variables
    =========
    integral_t3_ptr_vec_in : vector<IntegralTri3*>
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
    std::vector<IntegralTri3*> integral_t3_ptr_vec;
    std::vector<IntegralQuad4*> integral_q4_ptr_vec;

    // functions for computing domain integrals
    void evaluate_Ni();
    void evaluate_integral_Ni();
    void evaluate_integral_derivative_Ni_x();
    void evaluate_integral_Ni_Nj();
    void evaluate_integral_Ni_derivative_Nj_x();
    void evaluate_integral_div_Ni_dot_div_Nj();
    void evaluate_integral_Ni_Nj_derivative_Nk_x();

    // functions for computing boundary integrals
    void evaluate_boundary_Ni();
    void evaluate_integral_boundary_Ni();
    void evaluate_integral_boundary_Ni_Nj();

    // default constructor
    IntegralGroup() {}

    // constructor
    IntegralGroup(std::vector<IntegralTri3*> integral_t3_ptr_vec_in, std::vector<IntegralQuad4*> integral_q4_ptr_vec_in)
    {
        integral_t3_ptr_vec = integral_t3_ptr_vec_in;
        integral_q4_ptr_vec = integral_q4_ptr_vec_in;
    }

};

void IntegralGroup::evaluate_Ni()
{
    /*

    Calculates test function values and other properties.
    Must be called before domain integrals are evaluated.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_Ni();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_Ni();
    }

}

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
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
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
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_integral_derivative_Ni_x();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_integral_derivative_Ni_x();
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
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_Nj();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
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
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_derivative_Nj_x();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_derivative_Nj_x();
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
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_integral_div_Ni_dot_div_Nj();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
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
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_Nj_derivative_Nk_x();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_Nj_derivative_Nk_x();
    }

}

void IntegralGroup::evaluate_boundary_Ni()
{
    /*

    Calculates test function values and other properties at the boundary.
    Must be called before bounary integrals are evaluated.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // evaluate integrals in each domain
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_boundary_Ni();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_boundary_Ni();
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
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_integral_boundary_Ni();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
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
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_integral_boundary_Ni_Nj();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_integral_boundary_Ni_Nj();
    }

}

}

#endif
