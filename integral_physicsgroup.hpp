#ifndef INTEGRAL_PHYSICSGROUP
#define INTEGRAL_PHYSICSGROUP
#include <vector>
#include "integral_quad4.hpp"

class IntegralPhysicsGroup
{
    /*

    Groups test function integrals (N) that are used in the same physics.

    Variables
    =========
    integral_ptr_vec_in : vector<IntegralQuad4*>
        vector with pointers to BoundaryQuad4 objects.

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

    */


    public:

    // variables
    std::vector<IntegralQuad4*> integral_ptr_vec;

    // functions
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

    // default constructor
    IntegralPhysicsGroup()
    {

    }

    // constructor
    IntegralPhysicsGroup(std::vector<IntegralQuad4*> integral_ptr_vec_in)
    {
        integral_ptr_vec = integral_ptr_vec_in;
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
    for (auto integral_ptr : integral_ptr_vec)
    {
        integral_ptr->evaluate_Ni_derivative();
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
    for (auto integral_ptr : integral_ptr_vec)
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
    for (auto integral_ptr : integral_ptr_vec)
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
    for (auto integral_ptr : integral_ptr_vec)
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
    for (auto integral_ptr : integral_ptr_vec)
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
    for (auto integral_ptr : integral_ptr_vec)
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
    for (auto integral_ptr : integral_ptr_vec)
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
    for (auto integral_ptr : integral_ptr_vec)
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
    for (auto integral_ptr : integral_ptr_vec)
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
    for (auto integral_ptr : integral_ptr_vec)
    {
        integral_ptr->evaluate_integral_Ni_Nj_derivative_Nk_y();
    }

}

#endif
