#ifndef INTEGRAL_QUAD4
#define INTEGRAL_QUAD4
#include <utility>
#include <vector>
#include "Eigen/Eigen"
#include "boundary_quad4.hpp"
#include "mesh_quad4.hpp"
#include "container_typedef.hpp"

class IntegralQuad4
{
    /*

    Test function (N) integrals for quad4 mesh elements.

    Variables
    =========
    mesh_q4_in : MeshQuad4
        struct with mesh data.
    boundary_q4_in : BoundaryQuad4
        object with boundary condition data.

    Functions
    =========
    evaluate_Ni_derivative : void
        Calculates test functions and their derivatives.
        Must be called before integrals are evaluated.
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

    Notes
    ====
    The calculated integrals are stored in nested vectors.
    Values can be accessed from each vector using the following pattern:
        integral_vec[element_did][i][j]...
    wherein element_did is the domain element ID and i, j, ... are indices.
    The calculated boundary integrals are stored in nested maps.
    Values can be accessed from each map using the following pattern:
        integral_map[element_did][boundary_key][i][j]...
    wherein boundary_key depends on the local point IDs along a given boundary.

    */

    public:
    
    // mesh
    MeshQuad4 *mesh_q4_ptr;
    BoundaryQuad4 *boundary_q4_ptr;

    // vectors with test functions and derivatives
    Vector2D jacobian_determinant_vec;
    Vector3D N_vec;
    Vector3D derivative_N_x_vec;
    Vector3D derivative_N_y_vec;

    // maps with test functions for boundaries
    MapIntIntVector1D boundary_jacobian_determinant_map;
    MapIntIntVector2D boundary_N_map;

    // vectors with integrals
    Vector2D integral_Ni_vec;
    Vector2D integral_derivative_Ni_x_vec;
    Vector2D integral_derivative_Ni_y_vec;
    Vector3D integral_Ni_Nj_vec;
    Vector3D integral_Ni_derivative_Nj_x_vec;
    Vector3D integral_Ni_derivative_Nj_y_vec;
    Vector3D integral_div_Ni_dot_div_Nj_vec;
    Vector4D integral_Ni_Nj_derivative_Nk_x_vec;
    Vector4D integral_Ni_Nj_derivative_Nk_y_vec;

    // maps with integrals for boundaries
    MapIntIntVector1D boundary_integral_Ni_map;
    MapIntIntVector2D boundary_integral_Ni_Nj_map;

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
    void evaluate_boundary_integral_Ni();
    void evaluate_boundary_integral_Ni_Nj();

    // default constructor
    IntegralQuad4()
    {

    }

    // constructor
    IntegralQuad4(MeshQuad4 &mesh_q4_in, BoundaryQuad4 &boundary_q4_in)
    {
        mesh_q4_ptr = &mesh_q4_in;
        boundary_q4_ptr = &boundary_q4_in;
    }

};

void IntegralQuad4::evaluate_Ni_derivative()
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

    // integration points
    // dimensionless coordinates if element is scaled to [-1, 1] * [-1, 1]
    const double M_1_SQRT_3 = 1./sqrt(3);
    double a_arr[4] = {+M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3};
    double b_arr[4] = {+M_1_SQRT_3, -M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3};

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_q4_ptr->num_element_domain; element_did++)
    {

        // initialize
        Vector1D jacobian_determinant_part_ml_vec;
        Vector2D N_part_ml_vec;
        Vector2D derivative_N_x_part_ml_vec;
        Vector2D derivative_N_y_part_ml_vec;

        // get global ID of points around element
        int p0_gid = mesh_q4_ptr->element_p0_gid_vec[element_did];
        int p1_gid = mesh_q4_ptr->element_p1_gid_vec[element_did];
        int p2_gid = mesh_q4_ptr->element_p2_gid_vec[element_did];
        int p3_gid = mesh_q4_ptr->element_p3_gid_vec[element_did];

        // get domain ID of points
        int p0_did = mesh_q4_ptr->point_gid_to_did_map[p0_gid];
        int p1_did = mesh_q4_ptr->point_gid_to_did_map[p1_gid];
        int p2_did = mesh_q4_ptr->point_gid_to_did_map[p2_gid];
        int p3_did = mesh_q4_ptr->point_gid_to_did_map[p3_gid];

        // get x values of points
        double x0 = mesh_q4_ptr->point_position_x_vec[p0_did];
        double x1 = mesh_q4_ptr->point_position_x_vec[p1_did];
        double x2 = mesh_q4_ptr->point_position_x_vec[p2_did];
        double x3 = mesh_q4_ptr->point_position_x_vec[p3_did];

        // get y values of points
        double y0 = mesh_q4_ptr->point_position_y_vec[p0_did];
        double y1 = mesh_q4_ptr->point_position_y_vec[p1_did];
        double y2 = mesh_q4_ptr->point_position_y_vec[p2_did];
        double y3 = mesh_q4_ptr->point_position_y_vec[p3_did];

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 4; indx_l++)
        {

            // initialize
            Vector1D N_part_mli_vec;
            Vector1D derivative_N_x_part_mli_vec;
            Vector1D derivative_N_y_part_mli_vec;

            // get a and b values where function is evaluated
            double a = a_arr[indx_l];
            double b = b_arr[indx_l];

            // get derivatives of x and y with respect to a and b
            double derivative_x_a = 0.25*(b*x0 - b*x1 + b*x2 - b*x3 - x0 - x1 + x2 + x3);
            double derivative_x_b = 0.25*(a*x0 - a*x1 + a*x2 - a*x3 - x0 + x1 + x2 - x3);
            double derivative_y_a = 0.25*(b*y0 - b*y1 + b*y2 - b*y3 - y0 - y1 + y2 + y3);
            double derivative_y_b = 0.25*(a*y0 - a*y1 + a*y2 - a*y3 - y0 + y1 + y2 - y3);

            // get jacobian and its inverse and determinant
            Eigen::Matrix2d jacobian_mat;
            jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
            Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
            double jacobian_determinant = jacobian_mat.determinant();

            // iterate for each test function
            for (int indx_i = 0; indx_i < 4; indx_i++)
            {
        
                // get test function N
                double N = 0.;
                switch (indx_i)
                {
                    case 0: N = 0.25*(1. - a)*(1. - b); break;
                    case 1: N = 0.25*(1. - a)*(1. + b); break;
                    case 2: N = 0.25*(1. + a)*(1. + b); break;
                    case 3: N = 0.25*(1. + a)*(1. - b); break;
                }

                // get derivatives of test function N
                double derivative_N_a = 0.;
                double derivative_N_b = 0.;
                switch (indx_i)
                {
                    case 0: derivative_N_a = +0.25*(b - 1.); derivative_N_b = +0.25*(a - 1.); break;
                    case 1: derivative_N_a = -0.25*(b + 1.); derivative_N_b = -0.25*(a - 1.); break;
                    case 2: derivative_N_a = +0.25*(b + 1.); derivative_N_b = +0.25*(a + 1.); break;
                    case 3: derivative_N_a = -0.25*(b - 1.); derivative_N_b = -0.25*(a + 1.); break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector2d derivative_N_ab_vec;
                derivative_N_ab_vec << derivative_N_a, derivative_N_b;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector2d derivative_N_xy_vec = derivative_N_ab_vec * jacobian_inverse_mat;
                double derivative_N_x = derivative_N_xy_vec.coeffRef(0);
                double derivative_N_y = derivative_N_xy_vec.coeffRef(1);

                // store in vectors
                N_part_mli_vec.push_back(N);
                derivative_N_x_part_mli_vec.push_back(derivative_N_x);
                derivative_N_y_part_mli_vec.push_back(derivative_N_y);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            derivative_N_x_part_ml_vec.push_back(derivative_N_x_part_mli_vec);
            derivative_N_y_part_ml_vec.push_back(derivative_N_y_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_vec.push_back(jacobian_determinant_part_ml_vec);
        N_vec.push_back(N_part_ml_vec);
        derivative_N_x_vec.push_back(derivative_N_x_part_ml_vec);
        derivative_N_y_vec.push_back(derivative_N_y_part_ml_vec);
        
    }

}

void IntegralQuad4::evaluate_boundary_Ni_derivative()
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

    // initialize 2D integration points (within element)
    const double M_1_SQRT_3 = 1./sqrt(3);
    VectorDouble sample_vec = {-M_1_SQRT_3, +M_1_SQRT_3};
    VectorDouble positive_vec = {+1., +1.};
    VectorDouble negative_vec = {-1., -1.};

    // initialize 1D integration points (along edge)
    // use capital letters (e.g., X, A) to denote variables along an edge
    VectorDouble A_vec = {-M_1_SQRT_3, +M_1_SQRT_3};

    // iterate for each element with a flux-type boundary condition
    for (int indx_m = 0; indx_m < boundary_q4_ptr->num_element_flux_domain; indx_m++)
    {

        // get element global ID and point local IDs
        int element_gid = boundary_q4_ptr->element_flux_gid_vec[indx_m];
        int point_lid_a = boundary_q4_ptr->element_flux_pa_lid_vec[indx_m];
        int point_lid_b = boundary_q4_ptr->element_flux_pb_lid_vec[indx_m];

        // get boundary key
        // use a symmetric pairing function
        // this indicates which edge the boundary is applied on
        int helper_num = point_lid_a + point_lid_b + 1;
        int boundary_key = (helper_num*helper_num - helper_num % 2)/4 + std::min(point_lid_a, point_lid_b);

        // get coordinates of points defining the boundary

        // get element local ID
        int element_did = mesh_q4_ptr->element_gid_to_did_map[element_gid];

        // get global ID of points around element
        int p0_gid = mesh_q4_ptr->element_p0_gid_vec[element_did];
        int p1_gid = mesh_q4_ptr->element_p1_gid_vec[element_did];
        int p2_gid = mesh_q4_ptr->element_p2_gid_vec[element_did];
        int p3_gid = mesh_q4_ptr->element_p3_gid_vec[element_did];

        // get domain ID of points
        int p0_did = mesh_q4_ptr->point_gid_to_did_map[p0_gid];
        int p1_did = mesh_q4_ptr->point_gid_to_did_map[p1_gid];
        int p2_did = mesh_q4_ptr->point_gid_to_did_map[p2_gid];
        int p3_did = mesh_q4_ptr->point_gid_to_did_map[p3_gid];

        // get x values of points
        double x0 = mesh_q4_ptr->point_position_x_vec[p0_did];
        double x1 = mesh_q4_ptr->point_position_x_vec[p1_did];
        double x2 = mesh_q4_ptr->point_position_x_vec[p2_did];
        double x3 = mesh_q4_ptr->point_position_x_vec[p3_did];
        double x_arr[4] = {x0, x1, x2, x3};

        // get y values of points
        double y0 = mesh_q4_ptr->point_position_y_vec[p0_did];
        double y1 = mesh_q4_ptr->point_position_y_vec[p1_did];
        double y2 = mesh_q4_ptr->point_position_y_vec[p2_did];
        double y3 = mesh_q4_ptr->point_position_y_vec[p3_did];
        double y_arr[4] = {y0, y1, y2, y3};

        // determine integration points used in gaussian integration

        // get 2D integration points
        VectorDouble a_vec;
        VectorDouble b_vec;
        switch (boundary_key)
        {
            case  1: a_vec = negative_vec; b_vec = sample_vec; break;  // left edge (0-1)
            case  5: a_vec = sample_vec; b_vec = positive_vec; break;  // top edge (1-2)
            case 11: a_vec = positive_vec; b_vec = sample_vec; break;  // right edge (2-3)
            case  4: a_vec = sample_vec; b_vec = negative_vec; break;  // bottom edge (3-0)
        }

        // initialize
        Vector1D jacobian_determinant_part_l_vec;
        Vector2D N_part_l_vec;

        // iterate for each integration point
        for (int indx_l = 0; indx_l < 2; indx_l++)
        {

            // initialize
            Vector1D N_part_li_vec;

            // get jacobian
            // evaluate using points along 1D edge (X; A)

            // transform 2D edge coordinates (x, y) into 1D (X)
            // get (x, y) along edge
            double xa = x_arr[point_lid_a];
            double xb = x_arr[point_lid_b];
            double ya = y_arr[point_lid_a];
            double yb = y_arr[point_lid_b];

            // let Xa = 0. Xb will be the distance from a to b
            double Xa = 0.;
            double Xb = sqrt((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya));

            // get derivatives of X with respect to A
            double derivative_X_A = 0.5*(Xb - Xa);

            // get jacobian and its inverse and determinant
            double jacobian_inverse = 1./derivative_X_A;
            double jacobian_determinant = derivative_X_A;

            // get test function values
            // evaluate using points in 2D element (x, y; a, b)

            // get a and b values where function is evaluated
            double a = a_vec[indx_l];
            double b = b_vec[indx_l];

            // iterate for each test function
            for (int indx_i = 0; indx_i < 4; indx_i++)
            {

                // get test function N
                double N = 0.;
                switch (indx_i)
                {
                    case 0: N = 0.25*(1. - a)*(1. - b); break;
                    case 1: N = 0.25*(1. - a)*(1. + b); break;
                    case 2: N = 0.25*(1. + a)*(1. + b); break;
                    case 3: N = 0.25*(1. + a)*(1. - b); break;
                }

                // store in vectors
                N_part_li_vec.push_back(N);

            }

            // store in vectors
            jacobian_determinant_part_l_vec.push_back(jacobian_determinant);
            N_part_l_vec.push_back(N_part_li_vec);

        }

        // store in vectors
        boundary_jacobian_determinant_map[element_did][boundary_key] = jacobian_determinant_part_l_vec;
        boundary_N_map[element_did][boundary_key] = N_part_l_vec;

    }

}

void IntegralQuad4::evaluate_integral_Ni()
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

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_q4_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector1D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * N_vec[element_did][indx_l][indx_i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_Ni_vec.push_back(integral_part_i_vec);

    }

}

void IntegralQuad4::evaluate_integral_derivative_Ni_x()
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
    
    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_q4_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector1D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * derivative_N_x_vec[element_did][indx_l][indx_i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_derivative_Ni_x_vec.push_back(integral_part_i_vec);

    }

}

void IntegralQuad4::evaluate_integral_derivative_Ni_y()
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
    
    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_q4_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector1D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * derivative_N_y_vec[element_did][indx_l][indx_i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_derivative_Ni_y_vec.push_back(integral_part_i_vec);

    }

}

void IntegralQuad4::evaluate_integral_Ni_Nj()
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

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_q4_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
    Vector1D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 4; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * N_vec[element_did][indx_l][indx_i] * N_vec[element_did][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_Nj_vec.push_back(integral_part_i_vec);

    }

}

void IntegralQuad4::evaluate_integral_Ni_derivative_Nj_x()
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

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_q4_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
    Vector1D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 4; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * N_vec[element_did][indx_l][indx_i] * derivative_N_x_vec[element_did][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Nj_x_vec.push_back(integral_part_i_vec);

    }

}

void IntegralQuad4::evaluate_integral_Ni_derivative_Nj_y()
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

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_q4_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
    Vector1D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 4; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * N_vec[element_did][indx_l][indx_i] * derivative_N_y_vec[element_did][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Nj_y_vec.push_back(integral_part_i_vec);

    }

}

void IntegralQuad4::evaluate_integral_div_Ni_dot_div_Nj()
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

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_q4_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
    Vector1D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 4; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * (derivative_N_x_vec[element_did][indx_l][indx_i] * derivative_N_x_vec[element_did][indx_l][indx_j] + derivative_N_y_vec[element_did][indx_l][indx_i] * derivative_N_y_vec[element_did][indx_l][indx_j]);
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_div_Ni_dot_div_Nj_vec.push_back(integral_part_i_vec);

    }

}

void IntegralQuad4::evaluate_integral_Ni_Nj_derivative_Nk_x()
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

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_q4_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector3D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
    Vector2D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 4; indx_j++){
    Vector1D integral_part_ijk_vec;
    for (int indx_k = 0; indx_k < 4; indx_k++){

        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * N_vec[element_did][indx_l][indx_i] * N_vec[element_did][indx_l][indx_j] * derivative_N_x_vec[element_did][indx_l][indx_k];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_Nj_derivative_Nk_x_vec.push_back(integral_part_i_vec);

    }

}

void IntegralQuad4::evaluate_integral_Ni_Nj_derivative_Nk_y()
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

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_q4_ptr->num_element_domain; element_did++){  
    
    // iterate for each test function combination
    Vector3D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
    Vector2D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 4; indx_j++){
    Vector1D integral_part_ijk_vec;
    for (int indx_k = 0; indx_k < 4; indx_k++){

        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[element_did][indx_l] * N_vec[element_did][indx_l][indx_i] * N_vec[element_did][indx_l][indx_j] * derivative_N_y_vec[element_did][indx_l][indx_k];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_Nj_derivative_Nk_x_vec.push_back(integral_part_i_vec);

    }

}

void IntegralQuad4::evaluate_boundary_integral_Ni()
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

    // iterate for each element with a flux-type boundary condition
    for (int indx_m = 0; indx_m < boundary_q4_ptr->num_element_flux_domain; indx_m++)
    {

        // get element global ID and boundary key

        // get element global ID and point local IDs
        int element_gid = boundary_q4_ptr->element_flux_gid_vec[indx_m];
        int point_lid_a = boundary_q4_ptr->element_flux_pa_lid_vec[indx_m];
        int point_lid_b = boundary_q4_ptr->element_flux_pb_lid_vec[indx_m];

        // get boundary key
        // use a symmetric pairing function
        // this indicates which edge the boundary is applied on
        int helper_num = point_lid_a + point_lid_b + 1;
        int boundary_key = (helper_num*helper_num - helper_num % 2)/4 + std::min(point_lid_a, point_lid_b);

        // get element local ID
        int element_did = mesh_q4_ptr->element_gid_to_did_map[element_gid];
    
        // calculate integrals

        // iterate for each test function combination
        Vector1D integral_part_i_vec;
        for (int indx_i = 0; indx_i < 4; indx_i++){

            // iterate for each integration point
            double integral_value = 0;
            for (int indx_l = 0; indx_l < 2; indx_l++) 
            {
                integral_value += boundary_jacobian_determinant_map[element_did][boundary_key][indx_l] * boundary_N_map[element_did][boundary_key][indx_l][indx_i];
            }
            integral_part_i_vec.push_back(integral_value);
        
        }
        boundary_integral_Ni_map[element_did][boundary_key] = integral_part_i_vec;

    }

}

void IntegralQuad4::evaluate_boundary_integral_Ni_Nj()
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

    // iterate for each element with a flux-type boundary condition
    for (int indx_m = 0; indx_m < boundary_q4_ptr->num_element_flux_domain; indx_m++)
    {

        // get element global ID and boundary key

        // get element global ID and point local IDs
        int element_gid = boundary_q4_ptr->element_flux_gid_vec[indx_m];
        int point_lid_a = boundary_q4_ptr->element_flux_pa_lid_vec[indx_m];
        int point_lid_b = boundary_q4_ptr->element_flux_pb_lid_vec[indx_m];

        // get boundary key
        // use a symmetric pairing function
        // this indicates which edge the boundary is applied on
        int helper_num = point_lid_a + point_lid_b + 1;
        int boundary_key = (helper_num*helper_num - helper_num % 2)/4 + std::min(point_lid_a, point_lid_b);

        // get element local ID
        int element_did = mesh_q4_ptr->element_gid_to_did_map[element_gid];
    
        // calculate integrals

        // iterate for each test function combination
        Vector2D integral_part_i_vec;
        for (int indx_i = 0; indx_i < 4; indx_i++){
        Vector1D integral_part_ij_vec;
        for (int indx_j = 0; indx_j < 4; indx_j++){

            // iterate for each integration point
            double integral_value = 0;
            for (int indx_l = 0; indx_l < 2; indx_l++) 
            {
                integral_value += boundary_jacobian_determinant_map[element_did][boundary_key][indx_l] * boundary_N_map[element_did][boundary_key][indx_l][indx_i] * boundary_N_map[element_did][boundary_key][indx_l][indx_j];;
            }
            integral_part_ij_vec.push_back(integral_value);
        
        }
        integral_part_i_vec.push_back(integral_part_ij_vec);
        }
        boundary_integral_Ni_Nj_map[element_did][boundary_key] = integral_part_i_vec;

    }

}

#endif
