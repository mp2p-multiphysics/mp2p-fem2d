#ifndef INTEGRAL_PHYSICSGROUP
#define INTEGRAL_PHYSICSGROUP
#include <vector>
#include "integral_quad4.hpp"
#include "integral_tri3.hpp"

class IntegralPhysicsGroup
{
    /*

    Groups test function integrals (N) that are used in the same physics.

    Variables
    =========
    integral_t3_ptr_vec_in : vector<IntegralTri3*>
        vector with pointers to IntegralTri3 objects.
    integral_q4_ptr_vec_in : vector<IntegralQuad4*>
        vector with pointers to IntegralQuad4 objects.

    Functions
    =========
    evaluate_Ni_derivative : void
        Calculates test functions (N) and their derivatives.
        Must be called before integrals are evaluated.
    evaluate_integral_Ni : void
        Calculates the integral of Ni.
    evaluate_integral_derivative_Ni_x : void
        Calculates the integral of d(Ni)/dy.
    evaluate_integral_derivative_Ni_y : void
        Calculates the integral of d(Ni)/dx.
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
    evaluate_boundary_Ni_derivative : void
        Calculates test functions (N) and their derivatives at the boundaries.
        Must be called before integrals are evaluated.
    evaluate_boundary_normal : void
        Calculates normal vectors at the boundaries.
    evaluate_boundary_integral_Ni : void
        Calculates the integral of Ni along a boundary.
    evaluate_boundary_integral_Ni_Nj
        Calculates the integral of Ni * Nj along a boundary.

    */


    public:

    // variables
    std::vector<IntegralTri3*> integral_t3_ptr_vec;
    std::vector<IntegralQuad4*> integral_q4_ptr_vec;

    // functions for computing integrals
    void evaluate_Ni_derivative();
    void evaluate_integral_Ni();
    void evaluate_integral_derivative_Ni_x();
    void evaluate_integral_derivative_Ni_y();
    void evaluate_integral_Ni_Nj();
    void evaluate_integral_Ni_derivative_Nj_x();
    void evaluate_integral_Ni_derivative_Nj_y();
    void evaluate_integral_div_Ni_dot_div_Nj();
    void evaluate_integral_Ni_Nj_derivative_Nk_x();
    void evaluate_integral_Ni_Nj_derivative_Nk_y();

    // functions for computing integrals for boundaries
    void evaluate_boundary_Ni_derivative();
    void evaluate_boundary_normal();
    void evaluate_boundary_integral_Ni();
    void evaluate_boundary_integral_Ni_Nj();

    // default constructor
    IntegralPhysicsGroup()
    {

    }

    // constructor
    IntegralPhysicsGroup(std::vector<IntegralTri3*> integral_t3_ptr_vec_in, std::vector<IntegralQuad4*> integral_q4_ptr_vec_in)
    {
        integral_t3_ptr_vec = integral_t3_ptr_vec_in;
        integral_q4_ptr_vec = integral_q4_ptr_vec_in;
    }

};

void IntegralPhysicsGroup::evaluate_Ni_derivative()
{
    /*

    Calculates test functions (N) and their derivatives.
    Must be called before integrals are evaluated.

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
        integral_ptr->evaluate_Ni_derivative();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_Ni_derivative();
    }

}

void IntegralPhysicsGroup::evaluate_boundary_Ni_derivative()
{
    /*

    Calculates test functions (N) and their derivatives at the boundaries.
    Must be called before integrals are evaluated.

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
        integral_ptr->evaluate_boundary_Ni_derivative();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_boundary_Ni_derivative();
    }

}


void IntegralPhysicsGroup::evaluate_boundary_normal()
{
    /*

    Calculates normal vectors at the boundaries.

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
        integral_ptr->evaluate_boundary_normal();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_boundary_normal();
    }

}

void IntegralPhysicsGroup::evaluate_integral_Ni()
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

void IntegralPhysicsGroup::evaluate_integral_derivative_Ni_x()
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


void IntegralPhysicsGroup::evaluate_integral_derivative_Ni_y()
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
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_integral_derivative_Ni_y();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_integral_derivative_Ni_y();
    }

}

void IntegralPhysicsGroup::evaluate_integral_Ni_Nj()
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

void IntegralPhysicsGroup::evaluate_integral_Ni_derivative_Nj_x()
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


void IntegralPhysicsGroup::evaluate_integral_Ni_derivative_Nj_y()
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
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_derivative_Nj_y();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_derivative_Nj_y();
    }

}

void IntegralPhysicsGroup::evaluate_integral_div_Ni_dot_div_Nj()
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

void IntegralPhysicsGroup::evaluate_integral_Ni_Nj_derivative_Nk_x()
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

void IntegralPhysicsGroup::evaluate_integral_Ni_Nj_derivative_Nk_y()
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
    for (auto integral_ptr : integral_t3_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_Nj_derivative_Nk_y();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_Nj_derivative_Nk_y();
    }

}

void IntegralPhysicsGroup::evaluate_boundary_integral_Ni()
{
    /*

    Calculates the integral of Ni along a boundary.

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
        integral_ptr->evaluate_boundary_integral_Ni();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_boundary_integral_Ni();
    }

}

void IntegralPhysicsGroup::evaluate_boundary_integral_Ni_Nj()
{
    /*

    Calculates the integral of Ni * Nj along a boundary.

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
        integral_ptr->evaluate_boundary_integral_Ni_Nj();
    }
    for (auto integral_ptr : integral_q4_ptr_vec)
    {
        integral_ptr->evaluate_boundary_integral_Ni_Nj();
    }

}

#endif
