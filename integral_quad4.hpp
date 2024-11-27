#ifndef INTEGRAL_QUAD4
#define INTEGRAL_QUAD4
#include <vector>
#include "Eigen/Eigen"
#include "boundary_quad4.hpp"
#include "container_typedef.hpp"
#include "domain_quad4.hpp"

namespace FEM2D
{

class IntegralQuad4
{
    /*

    Test function (N) integrals for quad4 elements.

    Variables
    =========
    domain_in : DomainQuad4
        Domain where element integrals are calculated.
    boundary_in : BoundaryQuad4
        Boundaries where element integrals are calculated.

    Functions
    =========
    evaluate_Ni : void
        Calculates test function values and other properties.
        Must be called before domain integrals are evaluated.
    evaluate_integral_Ni : void
        Calculates the integral of Ni.
    evaluate_integral_derivative_Ni_x : void
        Calculates the integral of d(Ni)/dx.
    evaluate_integral_Ni_Nj : void
        Calculates the integral of Ni * Nj.
    evaluate_integral_Ni_derivative_Nj_x : void
        Calculates the integral of Ni * d(Nj)/dx.
    evaluate_integral_div_Ni_dot_div_Nj : void
        Calculates the integral of div(Ni) dot div(Nj).
    evaluate_integral_Ni_Nj_derivative_Nk_x : void
        Calculates the integral of Ni * Nj * d(Nk)/dx.
    evaluate_boundary_Ni : void
        Calculates test functions values and other properties at the boundary.
        Must be called before boundary integrals are evaluated.
    evaluate_integral_boundary_Ni : void
        Calculates the integral of Ni at the boundary.
    evaluate_integral_boundary_Ni_Nj : void
        Calculates the integral of Ni * Nj at the boundary.

    Notes
    ====
    The calculated integrals are stored in nested maps and vectors.
    Values of the domain integrals can be accessed from each vector using the following pattern:
        integral_vec[edid][i][j]...
        wherein edid is the domain element ID and i, j, ... are indices.
    Values of the boundary integrals can be accessed from each vector using the following pattern:
        integral_vec[edid][boundary_key][i][j]...
        wherein boundary_key is an int denoting the location of the boundary.
        For quad4 elements, the boundary_key is just the local ID of the point (0 or 1).

    */

    public:
    
    // domain and boundary
    DomainQuad4 *domain_ptr;
    BoundaryQuad4 *boundary_ptr;

    // vectors with domain test functions
    // index as follows: [edid][integration_point][i]
    Vector3D Ni_vec;
    Vector3D derivative_Ni_x_vec;
    Vector3D derivative_Ni_y_vec;
    Vector2D jacobian_determinant_vec;

    // vectors with domain integrals
    // index as follows: [edid][i][j][k]
    Vector2D integral_Ni_vec;
    Vector2D integral_derivative_Ni_x_vec;
    Vector2D integral_derivative_Ni_y_vec;
    Vector3D integral_Ni_Nj_vec;
    Vector3D integral_Ni_derivative_Nj_x_vec;
    Vector3D integral_Ni_derivative_Nj_y_vec;
    Vector3D integral_div_Ni_dot_div_Nj_vec;
    Vector4D integral_Ni_Nj_derivative_Nk_x_vec;
    Vector4D integral_Ni_Nj_derivative_Nk_y_vec;

    // vectors with boundary test functions
    MapVector4D boundary_Ni_vec;  // [edid][boundary_key][integration_point][i]
    MapVector2D boundary_normal_x_vec;  // [edid][boundary_key]
    MapVector2D boundary_normal_y_vec;  // [edid][boundary_key]
    MapVector2D boundary_jacobian_determinant_vec;  // [edid][boundary_key]

    // vectors with boundary integrals
    // index as follows: [edid][boundary_key][i][j]
    MapVector3D integral_boundary_Ni_vec;
    MapVector4D integral_boundary_Ni_Nj_vec;

    // functions for computing domain integrals
    void evaluate_Ni();
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
    void evaluate_boundary_Ni();
    void evaluate_integral_boundary_Ni();
    void evaluate_integral_boundary_Ni_Nj();

    // default constructor
    IntegralQuad4() {}

    // constructor
    IntegralQuad4(DomainQuad4 &domain_in, BoundaryQuad4 &boundary_in)
    {
        
        // store domain and boundaries
        domain_ptr = &domain_in;
        boundary_ptr = &boundary_in;

        // evaluate test functions
        evaluate_Ni();
        evaluate_boundary_Ni();

    }

};

void IntegralQuad4::evaluate_Ni()
{
    /*

    Calculates test function values and other properties.
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // initialize
        Vector1D jacobian_determinant_part_ml_vec;
        Vector2D N_part_ml_vec;
        Vector2D derivative_N_x_part_ml_vec;
        Vector2D derivative_N_y_part_ml_vec;

        // get global ID of points around element
        int p0_pgid = domain_ptr->element_p0_pgid_vec[edid];
        int p1_pgid = domain_ptr->element_p1_pgid_vec[edid];
        int p2_pgid = domain_ptr->element_p2_pgid_vec[edid];
        int p3_pgid = domain_ptr->element_p3_pgid_vec[edid];

        // get domain ID of points
        int p0_pdid = domain_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_ptr->point_pgid_to_pdid_map[p3_pgid];

        // get x values of points
        double x0 = domain_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_ptr->point_position_x_vec[p3_pdid];

        // get y values of points
        double y0 = domain_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_ptr->point_position_y_vec[p3_pdid];

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
        Ni_vec.push_back(N_part_ml_vec);
        derivative_Ni_x_vec.push_back(derivative_N_x_part_ml_vec);
        derivative_Ni_y_vec.push_back(derivative_N_y_part_ml_vec);
 
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector1D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i];
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector1D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * derivative_Ni_x_vec[edid][indx_l][indx_i];
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector1D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * derivative_Ni_y_vec[edid][indx_l][indx_i];
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
    Vector1D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 4; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * Ni_vec[edid][indx_l][indx_j];
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
    Vector1D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 4; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * derivative_Ni_x_vec[edid][indx_l][indx_j];
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
    Vector1D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 4; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * derivative_Ni_y_vec[edid][indx_l][indx_j];
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    Vector2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < 4; indx_i++){  
    Vector1D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < 4; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < 4; indx_l++) 
        {
            integral_value += jacobian_determinant_vec[edid][indx_l] * (derivative_Ni_x_vec[edid][indx_l][indx_i] * derivative_Ni_x_vec[edid][indx_l][indx_j] + derivative_Ni_y_vec[edid][indx_l][indx_i] * derivative_Ni_y_vec[edid][indx_l][indx_j]);
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
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
            integral_value += jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * Ni_vec[edid][indx_l][indx_j] * derivative_Ni_x_vec[edid][indx_l][indx_k];
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
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
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
            integral_value += jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * Ni_vec[edid][indx_l][indx_j] * derivative_Ni_y_vec[edid][indx_l][indx_k];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_Nj_derivative_Nk_y_vec.push_back(integral_part_i_vec);

    }

}

void IntegralQuad4::evaluate_boundary_Ni()
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

    // initialize 2D integration points (within element)
    const double M_1_SQRT_3 = 1./sqrt(3);
    VectorDouble sample_vec = {-M_1_SQRT_3, +M_1_SQRT_3};
    VectorDouble positive_vec = {+1., +1.};
    VectorDouble negative_vec = {-1., -1.};

    // initialize 1D integration points (along edge)
    // use capital letters (e.g., X, A) to denote variables along an edge
    VectorDouble A_vec = {-M_1_SQRT_3, +M_1_SQRT_3};

    // iterate through each boundary
    for (int bid = 0; bid < boundary_ptr->num_boundary; bid++)
    {

        // get the element ID
        int egid = boundary_ptr->boundary_egid_vec[bid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get points around element
        int p0_pgid = domain_ptr->element_p0_pgid_vec[edid];
        int p1_pgid = domain_ptr->element_p1_pgid_vec[edid];
        int p2_pgid = domain_ptr->element_p2_pgid_vec[edid];
        int p3_pgid = domain_ptr->element_p3_pgid_vec[edid];
        int p0_pdid = domain_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_ptr->point_pgid_to_pdid_map[p3_pgid];

        // get boundary key
        // use a symmetric pairing function
        // this indicates which edge the boundary is applied on
        int pa_plid = boundary_ptr->boundary_pa_plid_vec[bid];
        int pb_plid = boundary_ptr->boundary_pb_plid_vec[bid];
        int helper_num = pa_plid + pb_plid + 1;
        int boundary_key = (helper_num*helper_num - helper_num % 2)/4 + std::min(pa_plid, pb_plid);

        // get x values of points
        double x0 = domain_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_ptr->point_position_x_vec[p3_pdid];
        double x_arr[4] = {x0, x1, x2, x3};

        // get y values of points
        double y0 = domain_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_ptr->point_position_y_vec[p3_pdid];
        double y_arr[4] = {y0, y1, y2, y3};

        // get jacobian
        // evaluate using points along 1D edge (X; A)

        // transform 2D edge coordinates (x, y) into 1D (X)
        // get (x, y) along edge
        double xa = x_arr[pa_plid];
        double xb = x_arr[pb_plid];
        double ya = y_arr[pa_plid];
        double yb = y_arr[pb_plid];

        // let Xa = 0. Xb will be the distance from a to b
        double Xa = 0.;
        double Xb = sqrt((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya));

        // get derivatives of X with respect to A
        double derivative_X_A = 0.5*(Xb - Xa);

        // get jacobian and its inverse and determinant
        double jacobian_inverse = 1./derivative_X_A;
        double jacobian_determinant = derivative_X_A;

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
            N_part_l_vec.push_back(N_part_li_vec);

        }

        // store in vectors
        boundary_jacobian_determinant_vec[edid][boundary_key] = jacobian_determinant;
        boundary_Ni_vec[edid][boundary_key] = N_part_l_vec;
        
        // determine outward unit normal vector

        // get vector along boundary pointing counterclocwise
        double x_start = 0.; double y_start = 0.;
        double x_end = 0.; double y_end = 0.;    
        switch (boundary_key)
        {
            case  1:  // left edge (0-1)
                x_start = x0; y_start = y0;
                x_end = x1; y_end = y1;
            break;  
            case  5:  // top edge (1-2)
                x_start = x1; y_start = y1;
                x_end = x2; y_end = y2;
            break;  
            case 11:  // right edge (2-3)
                x_start = x2; y_start = y2;
                x_end = x3; y_end = y3;
            break;
            case  4:  // bottom edge (3-0)
                x_start = x3; y_start = y3;
                x_end = x0; y_end = y0;
            break; 
        }

        // calculate normal to boundary vector
        // should be pointing out of element
        double boundarynormal_x = y_end - y_start;
        double boundarynormal_y = x_start - x_end;

        // calculate unit normal vector components
        double boundarynormal_mag = sqrt(boundarynormal_x*boundarynormal_x + boundarynormal_y*boundarynormal_y);
        boundary_normal_x_vec[edid][boundary_key] = boundarynormal_x/boundarynormal_mag;
        boundary_normal_y_vec[edid][boundary_key] = boundarynormal_y/boundarynormal_mag;

    }

}

void IntegralQuad4::evaluate_integral_boundary_Ni()
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

    // iterate through each boundary
    for (int bid = 0; bid < boundary_ptr->num_boundary; bid++)
    {

        // get the element ID
        int egid = boundary_ptr->boundary_egid_vec[bid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get boundary key
        // use a symmetric pairing function
        // this indicates which edge the boundary is applied on
        int pa_plid = boundary_ptr->boundary_pa_plid_vec[bid];
        int pb_plid = boundary_ptr->boundary_pb_plid_vec[bid];
        int helper_num = pa_plid + pb_plid + 1;
        int boundary_key = (helper_num*helper_num - helper_num % 2)/4 + std::min(pa_plid, pb_plid);

        // iterate for each test function combination
        Vector1D integral_part_i_vec;
        for (int indx_i = 0; indx_i < 4; indx_i++){

            // iterate for each integration point
            double integral_value = 0;
            for (int indx_l = 0; indx_l < 2; indx_l++) 
            {
                integral_value += boundary_jacobian_determinant_vec[edid][boundary_key] * boundary_Ni_vec[edid][boundary_key][indx_l][indx_i];
            }
            integral_part_i_vec.push_back(integral_value);
        
        }
        integral_boundary_Ni_vec[edid][boundary_key] = integral_part_i_vec;

    }

}

void IntegralQuad4::evaluate_integral_boundary_Ni_Nj()
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

    // iterate through each boundary
    for (int bid = 0; bid < boundary_ptr->num_boundary; bid++)
    {

        // get the element ID
        int egid = boundary_ptr->boundary_egid_vec[bid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get boundary key
        // use a symmetric pairing function
        // this indicates which edge the boundary is applied on
        int pa_plid = boundary_ptr->boundary_pa_plid_vec[bid];
        int pb_plid = boundary_ptr->boundary_pb_plid_vec[bid];
        int helper_num = pa_plid + pb_plid + 1;
        int boundary_key = (helper_num*helper_num - helper_num % 2)/4 + std::min(pa_plid, pb_plid);

        // iterate for each test function combination
        Vector2D integral_part_i_vec;
        for (int indx_i = 0; indx_i < 4; indx_i++){
        Vector1D integral_part_ij_vec;
        for (int indx_j = 0; indx_j < 4; indx_j++){

            // iterate for each integration point
            double integral_value = 0;
            for (int indx_l = 0; indx_l < 2; indx_l++) 
            {
                integral_value += boundary_jacobian_determinant_vec[edid][boundary_key] * boundary_Ni_vec[edid][boundary_key][indx_l][indx_i] * boundary_Ni_vec[edid][boundary_key][indx_l][indx_j];
            }
            integral_part_ij_vec.push_back(integral_value);
        
        }
        integral_part_i_vec.push_back(integral_part_ij_vec);
        }
        integral_boundary_Ni_Nj_vec[edid][boundary_key] = integral_part_i_vec;

    }

}

}

#endif
