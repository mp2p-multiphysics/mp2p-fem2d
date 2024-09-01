#ifndef MODULE_NAVIERSTOKES_STEADY_QUAD4QUAD8
#define MODULE_NAVIERSTOKES_STEADY_QUAD4QUAD8
#include <vector>
#include "Eigen/Eigen"
#include "integral_quad4quad8.hpp"
#include "grid_quad4.hpp"
#include "grid_quad8.hpp"
#include "scalar_quad8.hpp"

struct NavierStokesSteadyQuad4Quad8IntegralStruct
{
    
    // integrals with one shape function
    std::vector<std::vector<double>> integral_Mi_quad8_vec;
    
    // integrals with two shape functions
    std::vector<std::vector<std::vector<double>>> integral_Ni_quad4_derivative_Mj_quad8_x_vec;
    std::vector<std::vector<std::vector<double>>> integral_Ni_quad4_derivative_Mj_quad8_y_vec;
    std::vector<std::vector<std::vector<double>>> integral_Mi_quad8_derivative_Nj_quad4_x_vec;
    std::vector<std::vector<std::vector<double>>> integral_Mi_quad8_derivative_Nj_quad4_y_vec;
    std::vector<std::vector<std::vector<double>>> integral_div_Mi_quad8_dot_div_Mj_quad8_vec;
    std::vector<std::vector<std::vector<double>>> integral_Mi_quad8_Mj_quad8_vec;

    // integrals with three shape functions
    std::vector<std::vector<std::vector<std::vector<double>>>> integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_vec;
    std::vector<std::vector<std::vector<std::vector<double>>>> integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_vec;

};

class NavierStokesSteadyQuad4Quad8Class
{

    public:

    // variables
    GridQuad8Struct gq8s;
    GridQuad4Struct gq4s;
    BoundaryQuad8Struct u_bq8s;
    BoundaryQuad8Struct v_bq8s;
    BoundaryQuad4Struct p_bq4s;

    // functions
    void integral_evaluate();
    void matrix_fill(
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec,
        ScalarQuad8Class &density_sq8c, ScalarQuad8Class &viscosity_sq8c, double gravity_x, double gravity_y,
        int start_id_u, int start_id_v, int start_id_p
    );

    // struct with integrals
    NavierStokesSteadyQuad4Quad8IntegralStruct nssq4q8is;

    // constructor
    NavierStokesSteadyQuad4Quad8Class(GridQuad8Struct &gq8s_in, GridQuad4Struct &gq4s_in, BoundaryQuad8Struct &u_bq8s_in, BoundaryQuad8Struct &v_bq8s_in, BoundaryQuad4Struct &p_bq4s_in)
    {
        gq8s = gq8s_in;
        gq4s = gq4s_in;
        u_bq8s = u_bq8s_in;
        v_bq8s = v_bq8s_in;
        p_bq4s = p_bq4s_in;
    }

};

void NavierStokesSteadyQuad4Quad8Class::integral_evaluate()
{
    
    // iterate for each element
    for (int m = 0; m < gq4s.num_element; m++)
    {

        // get id of points around quad8 element
        int n0_quad8 = gq8s.element_p0_id_vec[m];
        int n1_quad8 = gq8s.element_p1_id_vec[m];
        int n2_quad8 = gq8s.element_p2_id_vec[m];
        int n3_quad8 = gq8s.element_p3_id_vec[m];
        int n4_quad8 = gq8s.element_p4_id_vec[m];
        int n5_quad8 = gq8s.element_p5_id_vec[m];
        int n6_quad8 = gq8s.element_p6_id_vec[m];
        int n7_quad8 = gq8s.element_p7_id_vec[m];

        // get id of points around quad4 element
        int n0_quad4 = gq4s.element_p0_id_vec[m];
        int n1_quad4 = gq4s.element_p1_id_vec[m];
        int n2_quad4 = gq4s.element_p2_id_vec[m];
        int n3_quad4 = gq4s.element_p3_id_vec[m];

        // get x coordinates of points around quad8 element
        double x0_quad8 = gq8s.point_pos_x_vec[n0_quad8];
        double x1_quad8 = gq8s.point_pos_x_vec[n1_quad8];
        double x2_quad8 = gq8s.point_pos_x_vec[n2_quad8];
        double x3_quad8 = gq8s.point_pos_x_vec[n3_quad8];
        double x4_quad8 = gq8s.point_pos_x_vec[n4_quad8];
        double x5_quad8 = gq8s.point_pos_x_vec[n5_quad8];
        double x6_quad8 = gq8s.point_pos_x_vec[n6_quad8];
        double x7_quad8 = gq8s.point_pos_x_vec[n7_quad8];
        double x_quad8_arr[8] = {x0_quad8, x1_quad8, x2_quad8, x3_quad8, x4_quad8, x5_quad8, x6_quad8, x7_quad8};

        // get y coordinates of points around quad8 element
        double y0_quad8 = gq8s.point_pos_y_vec[n0_quad8];
        double y1_quad8 = gq8s.point_pos_y_vec[n1_quad8];
        double y2_quad8 = gq8s.point_pos_y_vec[n2_quad8];
        double y3_quad8 = gq8s.point_pos_y_vec[n3_quad8];
        double y4_quad8 = gq8s.point_pos_y_vec[n4_quad8];
        double y5_quad8 = gq8s.point_pos_y_vec[n5_quad8];
        double y6_quad8 = gq8s.point_pos_y_vec[n6_quad8];
        double y7_quad8 = gq8s.point_pos_y_vec[n7_quad8];
        double y_quad8_arr[8] = {y0_quad8, y1_quad8, y2_quad8, y3_quad8, y4_quad8, y5_quad8, y6_quad8, y7_quad8};

        // temporarily stores values
        double part_val = 0;

        // evaluate integral_Mi_quad8
        std::vector<double> integral_Mi_quad8_part_i_vec;
        for (int i = 0; i < 8; i++)
        {
            part_val = integral_Mi_quad8(i, x_quad8_arr, y_quad8_arr);
            integral_Mi_quad8_part_i_vec.push_back(part_val);
        }

        // evaluate integral_Ni_quad4_derivative_Mj_quad8_x
        std::vector<std::vector<double>> integral_Ni_quad4_derivative_Mj_quad8_x_part_i_vec;
        for (int i = 0; i < 4; i++)
        {
            std::vector<double> integral_Ni_quad4_derivative_Mj_quad8_x_part_ij_vec;
            for (int j = 0; j < 8; j++)
            {
                part_val = integral_Ni_quad4_derivative_Mj_quad8_x(i, j, x_quad8_arr, y_quad8_arr);
                integral_Ni_quad4_derivative_Mj_quad8_x_part_ij_vec.push_back(part_val);
            }
            integral_Ni_quad4_derivative_Mj_quad8_x_part_i_vec.push_back(integral_Ni_quad4_derivative_Mj_quad8_x_part_ij_vec);
        }

        // evaluate integral_Ni_quad4_derivative_Mj_quad8_y
        std::vector<std::vector<double>> integral_Ni_quad4_derivative_Mj_quad8_y_part_i_vec;
        for (int i = 0; i < 4; i++)
        {
            std::vector<double> integral_Ni_quad4_derivative_Mj_quad8_y_part_ij_vec;
            for (int j = 0; j < 8; j++)
            {
                part_val = integral_Ni_quad4_derivative_Mj_quad8_y(i, j, x_quad8_arr, y_quad8_arr);
                integral_Ni_quad4_derivative_Mj_quad8_y_part_ij_vec.push_back(part_val);
            }
            integral_Ni_quad4_derivative_Mj_quad8_y_part_i_vec.push_back(integral_Ni_quad4_derivative_Mj_quad8_y_part_ij_vec);
        }

        // evaluate integral_Mi_quad8_derivative_Nj_quad4_x
        std::vector<std::vector<double>> integral_Mi_quad8_derivative_Nj_quad4_x_part_i_vec;
        for (int i = 0; i < 8; i++)
        {
            std::vector<double> integral_Mi_quad8_derivative_Nj_quad4_x_part_ij_vec;
            for (int j = 0; j < 4; j++)
            {
                part_val = integral_Mi_quad8_derivative_Nj_quad4_x(i, j, x_quad8_arr, y_quad8_arr);
                integral_Mi_quad8_derivative_Nj_quad4_x_part_ij_vec.push_back(part_val);
            }
            integral_Mi_quad8_derivative_Nj_quad4_x_part_i_vec.push_back(integral_Mi_quad8_derivative_Nj_quad4_x_part_ij_vec);
        }

        // evaluate integral_Mi_quad8_derivative_Nj_quad4_y
        std::vector<std::vector<double>> integral_Mi_quad8_derivative_Nj_quad4_y_part_i_vec;
        for (int i = 0; i < 8; i++)
        {
            std::vector<double> integral_Mi_quad8_derivative_Nj_quad4_y_part_ij_vec;
            for (int j = 0; j < 4; j++)
            {
                part_val = integral_Mi_quad8_derivative_Nj_quad4_y(i, j, x_quad8_arr, y_quad8_arr);
                integral_Mi_quad8_derivative_Nj_quad4_y_part_ij_vec.push_back(part_val);
            }
            integral_Mi_quad8_derivative_Nj_quad4_y_part_i_vec.push_back(integral_Mi_quad8_derivative_Nj_quad4_y_part_ij_vec);
        }

        // evaluate integral_div_Mi_quad8_dot_div_Mj_quad8
        std::vector<std::vector<double>> integral_div_Mi_quad8_dot_div_Mj_quad8_part_i_vec;
        for (int i = 0; i < 8; i++)
        {
            std::vector<double> integral_div_Mi_quad8_dot_div_Mj_quad8_part_ij_vec;
            for (int j = 0; j < 8; j++)
            {
                part_val = integral_div_Mi_quad8_dot_div_Mj_quad8(i, j, x_quad8_arr, y_quad8_arr);
                integral_div_Mi_quad8_dot_div_Mj_quad8_part_ij_vec.push_back(part_val);
            }
            integral_div_Mi_quad8_dot_div_Mj_quad8_part_i_vec.push_back(integral_div_Mi_quad8_dot_div_Mj_quad8_part_ij_vec);
        }

        // evaluate integral_Mi_Mj_quad8
        std::vector<std::vector<double>> integral_Mi_quad8_Mj_quad8_part_i_vec;
        for (int i = 0; i < 8; i++)
        {
            std::vector<double> integral_Mi_quad8_Mj_quad8_part_ij_vec;
            for (int j = 0; j < 8; j++)
            {
                part_val = integral_Mi_quad8_Mj_quad8(i, j, x_quad8_arr, y_quad8_arr);
                integral_Mi_quad8_Mj_quad8_part_ij_vec.push_back(part_val);
            }
            integral_Mi_quad8_Mj_quad8_part_i_vec.push_back(integral_Mi_quad8_Mj_quad8_part_ij_vec);
        }

        // evaluate integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x
        std::vector<std::vector<std::vector<double>>> integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_part_i_vec;
        for (int i = 0; i < 8; i++)
        {
            std::vector<std::vector<double>> integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_part_ij_vec;
            for (int j = 0; j < 8; j++)
            {
                std::vector<double> integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_part_ijk_vec;
                for (int k = 0; k < 8; k++)
                {
                    part_val = integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x(i, j, k, x_quad8_arr, y_quad8_arr);
                    integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_part_ijk_vec.push_back(part_val);                    
                }
                integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_part_ij_vec.push_back(integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_part_ijk_vec); 
            }
            integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_part_i_vec.push_back(integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_part_ij_vec);
        }

        // evaluate integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y
        std::vector<std::vector<std::vector<double>>> integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_part_i_vec;
        for (int i = 0; i < 8; i++)
        {
            std::vector<std::vector<double>> integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_part_ij_vec;
            for (int j = 0; j < 8; j++)
            {
                std::vector<double> integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_part_ijk_vec;
                for (int k = 0; k < 8; k++)
                {
                    part_val = integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y(i, j, k, x_quad8_arr, y_quad8_arr);
                    integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_part_ijk_vec.push_back(part_val);                    
                }
                integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_part_ij_vec.push_back(integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_part_ijk_vec); 
            }
            integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_part_i_vec.push_back(integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_part_ij_vec);
        }

        // store integrals to struct
        nssq4q8is.integral_Mi_quad8_vec.push_back(integral_Mi_quad8_part_i_vec);
        nssq4q8is.integral_Ni_quad4_derivative_Mj_quad8_x_vec.push_back(integral_Ni_quad4_derivative_Mj_quad8_x_part_i_vec);
        nssq4q8is.integral_Ni_quad4_derivative_Mj_quad8_y_vec.push_back(integral_Ni_quad4_derivative_Mj_quad8_y_part_i_vec);
        nssq4q8is.integral_Mi_quad8_derivative_Nj_quad4_x_vec.push_back(integral_Mi_quad8_derivative_Nj_quad4_x_part_i_vec);
        nssq4q8is.integral_Mi_quad8_derivative_Nj_quad4_y_vec.push_back(integral_Mi_quad8_derivative_Nj_quad4_y_part_i_vec);
        nssq4q8is.integral_div_Mi_quad8_dot_div_Mj_quad8_vec.push_back(integral_div_Mi_quad8_dot_div_Mj_quad8_part_i_vec);
        nssq4q8is.integral_Mi_quad8_Mj_quad8_vec.push_back(integral_Mi_quad8_Mj_quad8_part_i_vec);
        nssq4q8is.integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_vec.push_back(integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_part_i_vec);
        nssq4q8is.integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_vec.push_back(integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_part_i_vec);

    }

}

void NavierStokesSteadyQuad4Quad8Class::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec,
    ScalarQuad8Class &density_sq8c, ScalarQuad8Class &viscosity_sq8c, double gravity_x, double gravity_y,
    int start_id_u, int start_id_v, int start_id_p
)
{

    // fill up matrix with NSE
    for (int m = 0; m < gq8s.num_element; m++)
    {

        // get id of points around quad8 element
        int n0_quad8 = gq8s.element_p0_id_vec[m];
        int n1_quad8 = gq8s.element_p1_id_vec[m];
        int n2_quad8 = gq8s.element_p2_id_vec[m];
        int n3_quad8 = gq8s.element_p3_id_vec[m];
        int n4_quad8 = gq8s.element_p4_id_vec[m];
        int n5_quad8 = gq8s.element_p5_id_vec[m];
        int n6_quad8 = gq8s.element_p6_id_vec[m];
        int n7_quad8 = gq8s.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0_quad8, n1_quad8, n2_quad8, n3_quad8, n4_quad8, n5_quad8, n6_quad8, n7_quad8};

        // get id of points around quad4 element
        int n0_quad4 = gq4s.element_p0_id_vec[m];
        int n1_quad4 = gq4s.element_p1_id_vec[m];
        int n2_quad4 = gq4s.element_p2_id_vec[m];
        int n3_quad4 = gq4s.element_p3_id_vec[m];
        int n_quad4_arr[4] = {n0_quad4, n1_quad4, n2_quad4, n3_quad4};

        // get u velocity values at points around element
        double u0_quad8 = x_last_iteration_vec[start_id_u + n0_quad8];
        double u1_quad8 = x_last_iteration_vec[start_id_u + n1_quad8];
        double u2_quad8 = x_last_iteration_vec[start_id_u + n2_quad8];
        double u3_quad8 = x_last_iteration_vec[start_id_u + n3_quad8];
        double u4_quad8 = x_last_iteration_vec[start_id_u + n4_quad8];
        double u5_quad8 = x_last_iteration_vec[start_id_u + n5_quad8];
        double u6_quad8 = x_last_iteration_vec[start_id_u + n6_quad8];
        double u7_quad8 = x_last_iteration_vec[start_id_u + n7_quad8];
        double u_quad8_arr[8] = {u0_quad8, u1_quad8, u2_quad8, u3_quad8, u4_quad8, u5_quad8, u6_quad8, u7_quad8};

        // get v velocity values at points around element
        double v0_quad8 = x_last_iteration_vec[start_id_v + n0_quad8];
        double v1_quad8 = x_last_iteration_vec[start_id_v + n1_quad8];
        double v2_quad8 = x_last_iteration_vec[start_id_v + n2_quad8];
        double v3_quad8 = x_last_iteration_vec[start_id_v + n3_quad8];
        double v4_quad8 = x_last_iteration_vec[start_id_v + n4_quad8];
        double v5_quad8 = x_last_iteration_vec[start_id_v + n5_quad8];
        double v6_quad8 = x_last_iteration_vec[start_id_v + n6_quad8];
        double v7_quad8 = x_last_iteration_vec[start_id_v + n7_quad8];
        double v_quad8_arr[8] = {v0_quad8, v1_quad8, v2_quad8, v3_quad8, v4_quad8, v5_quad8, v6_quad8, v7_quad8};

        // get density at points around element
        double density0_quad8 = density_sq8c.scalar_vec[n0_quad8];
        double density1_quad8 = density_sq8c.scalar_vec[n1_quad8];
        double density2_quad8 = density_sq8c.scalar_vec[n2_quad8];
        double density3_quad8 = density_sq8c.scalar_vec[n3_quad8];
        double density4_quad8 = density_sq8c.scalar_vec[n4_quad8];
        double density5_quad8 = density_sq8c.scalar_vec[n5_quad8];
        double density6_quad8 = density_sq8c.scalar_vec[n6_quad8];
        double density7_quad8 = density_sq8c.scalar_vec[n7_quad8];
        double density_quad8_arr[8] = {density0_quad8, density1_quad8, density2_quad8, density3_quad8, density4_quad8, density5_quad8, density6_quad8, density7_quad8};

        // get viscosity at points around element
        double viscosity0_quad8 = viscosity_sq8c.scalar_vec[n0_quad8];
        double viscosity1_quad8 = viscosity_sq8c.scalar_vec[n1_quad8];
        double viscosity2_quad8 = viscosity_sq8c.scalar_vec[n2_quad8];
        double viscosity3_quad8 = viscosity_sq8c.scalar_vec[n3_quad8];
        double viscosity4_quad8 = viscosity_sq8c.scalar_vec[n4_quad8];
        double viscosity5_quad8 = viscosity_sq8c.scalar_vec[n5_quad8];
        double viscosity6_quad8 = viscosity_sq8c.scalar_vec[n6_quad8];
        double viscosity7_quad8 = viscosity_sq8c.scalar_vec[n7_quad8];
        double viscosity_quad8_arr[8] = {viscosity0_quad8, viscosity1_quad8, viscosity2_quad8, viscosity3_quad8, viscosity4_quad8, viscosity5_quad8, viscosity6_quad8, viscosity7_quad8};

        // NSE x-component

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable [u0 ... uN v0 ... vN p0 ... pN]

        // coefficients of u
        for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {

            // calculate velocity integrals
            double integral_Mi_quad8_u_derivative_Mj_quad8_x = 0.;
            double integral_Mi_quad8_v_derivative_Mj_quad8_y = 0.;
            for (int k = 0; k < 8; k++)
            {
                integral_Mi_quad8_u_derivative_Mj_quad8_x += u_quad8_arr[k] * nssq4q8is.integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_vec[m][i][j][k];
                integral_Mi_quad8_v_derivative_Mj_quad8_y += v_quad8_arr[k] * nssq4q8is.integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_vec[m][i][j][k];
            }

            // fill up coefficient matrix
            a_mat.coeffRef(start_id_u + n_quad8_arr[i], start_id_u + n_quad8_arr[j]) += (
                  viscosity_quad8_arr[i]*nssq4q8is.integral_div_Mi_quad8_dot_div_Mj_quad8_vec[m][i][j]
                + density_quad8_arr[i]*integral_Mi_quad8_u_derivative_Mj_quad8_x
                + density_quad8_arr[i]*integral_Mi_quad8_v_derivative_Mj_quad8_y
            );

        }}

        // coefficients of p
        for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 4; j++) {
            a_mat.coeffRef(start_id_u + n_quad8_arr[i], start_id_p + n_quad4_arr[j]) += nssq4q8is.integral_Mi_quad8_derivative_Nj_quad4_x_vec[m][i][j];
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 8; i++) {
            b_vec.coeffRef(start_id_u + n_quad8_arr[i]) += density_quad8_arr[i]*gravity_x*nssq4q8is.integral_Mi_quad8_vec[m][i];
        }

        // NSE y-component

        // coefficients of v
        for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {

            // calculate velocity integrals
            double integral_Mi_quad8_u_derivative_Mj_quad8_x = 0.;
            double integral_Mi_quad8_v_derivative_Mj_quad8_y = 0.;
            for (int k = 0; k < 8; k++)
            {
                integral_Mi_quad8_u_derivative_Mj_quad8_x += u_quad8_arr[k] * nssq4q8is.integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x_vec[m][i][j][k];
                integral_Mi_quad8_v_derivative_Mj_quad8_y += v_quad8_arr[k] * nssq4q8is.integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y_vec[m][i][j][k];
            }

            // fill up coefficient matrix
            a_mat.coeffRef(start_id_v + n_quad8_arr[i], start_id_v + n_quad8_arr[j]) += (
                  viscosity_quad8_arr[i]*nssq4q8is.integral_div_Mi_quad8_dot_div_Mj_quad8_vec[m][i][j]
                + density_quad8_arr[i]*integral_Mi_quad8_u_derivative_Mj_quad8_x
                + density_quad8_arr[i]*integral_Mi_quad8_v_derivative_Mj_quad8_y
            );

        }}

        // coefficients of p
        for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 4; j++) {
            a_mat.coeffRef(start_id_v + n_quad8_arr[i], start_id_p + n_quad4_arr[j]) += nssq4q8is.integral_Mi_quad8_derivative_Nj_quad4_y_vec[m][i][j];
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 8; i++) {
            b_vec.coeffRef(start_id_v + n_quad8_arr[i]) += density_quad8_arr[i]*gravity_y*nssq4q8is.integral_Mi_quad8_vec[m][i];
        }

    }

    // fill up matrix with continuity equation
    for (int m = 0; m < gq4s.num_element; m++)
    {

        // get id of points around quad8 element
        int n0_quad8 = gq8s.element_p0_id_vec[m];
        int n1_quad8 = gq8s.element_p1_id_vec[m];
        int n2_quad8 = gq8s.element_p2_id_vec[m];
        int n3_quad8 = gq8s.element_p3_id_vec[m];
        int n4_quad8 = gq8s.element_p4_id_vec[m];
        int n5_quad8 = gq8s.element_p5_id_vec[m];
        int n6_quad8 = gq8s.element_p6_id_vec[m];
        int n7_quad8 = gq8s.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0_quad8, n1_quad8, n2_quad8, n3_quad8, n4_quad8, n5_quad8, n6_quad8, n7_quad8};

        // get id of points around quad4 element
        int n0_quad4 = gq4s.element_p0_id_vec[m];
        int n1_quad4 = gq4s.element_p1_id_vec[m];
        int n2_quad4 = gq4s.element_p2_id_vec[m];
        int n3_quad4 = gq4s.element_p3_id_vec[m];
        int n_quad4_arr[4] = {n0_quad4, n1_quad4, n2_quad4, n3_quad4};

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable [u0 ... uN v0 ... vN p0 ... pN]

        // coefficients of u
        for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            a_mat.coeffRef(start_id_p + n_quad4_arr[i], start_id_u + n_quad8_arr[j]) += nssq4q8is.integral_Ni_quad4_derivative_Mj_quad8_x_vec[m][i][j];
        }}

        // coefficients of v
        for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            a_mat.coeffRef(start_id_p + n_quad4_arr[i], start_id_v + n_quad8_arr[j]) += nssq4q8is.integral_Ni_quad4_derivative_Mj_quad8_y_vec[m][i][j];
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 4; i++) {
            b_vec.coeffRef(start_id_p + n_quad4_arr[i]) += 0.;
        }

    }

    // iterate for each flux boundary element (x-component of NSE)
    for (int k = 0; k < u_bq8s.num_element_flux; k++)
    {

        // get id of element
        int m = u_bq8s.element_flux_id_vec[k];

        // get id of points around quad8 element
        int n0_quad8 = gq8s.element_p0_id_vec[m];
        int n1_quad8 = gq8s.element_p1_id_vec[m];
        int n2_quad8 = gq8s.element_p2_id_vec[m];
        int n3_quad8 = gq8s.element_p3_id_vec[m];
        int n4_quad8 = gq8s.element_p4_id_vec[m];
        int n5_quad8 = gq8s.element_p5_id_vec[m];
        int n6_quad8 = gq8s.element_p6_id_vec[m];
        int n7_quad8 = gq8s.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0_quad8, n1_quad8, n2_quad8, n3_quad8, n4_quad8, n5_quad8, n6_quad8, n7_quad8};

        // get x coordinates of points around quad8 element
        double x0_quad8 = gq8s.point_pos_x_vec[n0_quad8];
        double x1_quad8 = gq8s.point_pos_x_vec[n1_quad8];
        double x2_quad8 = gq8s.point_pos_x_vec[n2_quad8];
        double x3_quad8 = gq8s.point_pos_x_vec[n3_quad8];
        double x4_quad8 = gq8s.point_pos_x_vec[n4_quad8];
        double x5_quad8 = gq8s.point_pos_x_vec[n5_quad8];
        double x6_quad8 = gq8s.point_pos_x_vec[n6_quad8];
        double x7_quad8 = gq8s.point_pos_x_vec[n7_quad8];
        double x_quad8_arr[8] = {x0_quad8, x1_quad8, x2_quad8, x3_quad8, x4_quad8, x5_quad8, x6_quad8, x7_quad8};

        // get y coordinates of points around quad8 element
        double y0_quad8 = gq8s.point_pos_y_vec[n0_quad8];
        double y1_quad8 = gq8s.point_pos_y_vec[n1_quad8];
        double y2_quad8 = gq8s.point_pos_y_vec[n2_quad8];
        double y3_quad8 = gq8s.point_pos_y_vec[n3_quad8];
        double y4_quad8 = gq8s.point_pos_y_vec[n4_quad8];
        double y5_quad8 = gq8s.point_pos_y_vec[n5_quad8];
        double y6_quad8 = gq8s.point_pos_y_vec[n6_quad8];
        double y7_quad8 = gq8s.point_pos_y_vec[n7_quad8];
        double y_quad8_arr[8] = {y0_quad8, y1_quad8, y2_quad8, y3_quad8, y4_quad8, y5_quad8, y6_quad8, y7_quad8};

        // get points where the boundary is applied
        int a = u_bq8s.element_flux_pa_loc_vec[k];  // 0 to 7
        int b = u_bq8s.element_flux_pb_loc_vec[k];  // 0 to 7
        int c = u_bq8s.element_flux_pb_loc_vec[k];  // 0 to 7

        // identify boundary type
        int config_id = u_bq8s.element_flux_config_id_vec[k];
        BoundaryConfigQuad8Struct bcq8s = u_bq8s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcq8s.boundary_type_str == "neumann")
        {

            // calculate distance between a-b-c
            // APPROXIMATION ONLY. FOR REVISION.
            double xa = x_quad8_arr[a];
            double xc = x_quad8_arr[c];
            double ya = y_quad8_arr[a];
            double yc = y_quad8_arr[c];
            double dist_pa_pb_pc = sqrt((xa - xc)*(xa - xc) + (ya - yc)*(ya - yc));

            // add to b_vec
            // NOT YET VERIFIED. PLACEHOLDER ONLY.
            b_vec.coeffRef(start_id_u + n_quad8_arr[a]) += 0.5 * bcq8s.boundary_parameter_vec[0] * dist_pa_pb_pc;
            b_vec.coeffRef(start_id_u + n_quad8_arr[b]) += 0.5 * bcq8s.boundary_parameter_vec[0] * dist_pa_pb_pc;
            b_vec.coeffRef(start_id_u + n_quad8_arr[c]) += 0.5 * bcq8s.boundary_parameter_vec[0] * dist_pa_pb_pc;

        }

    }

    // iterate for each flux boundary element (y-component of NSE)
    for (int k = 0; k < v_bq8s.num_element_flux; k++)
    {

        // get id of element
        int m = v_bq8s.element_flux_id_vec[k];

        // get id of points around quad8 element
        int n0_quad8 = gq8s.element_p0_id_vec[m];
        int n1_quad8 = gq8s.element_p1_id_vec[m];
        int n2_quad8 = gq8s.element_p2_id_vec[m];
        int n3_quad8 = gq8s.element_p3_id_vec[m];
        int n4_quad8 = gq8s.element_p4_id_vec[m];
        int n5_quad8 = gq8s.element_p5_id_vec[m];
        int n6_quad8 = gq8s.element_p6_id_vec[m];
        int n7_quad8 = gq8s.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0_quad8, n1_quad8, n2_quad8, n3_quad8, n4_quad8, n5_quad8, n6_quad8, n7_quad8};

        // get x coordinates of points around quad8 element
        double x0_quad8 = gq8s.point_pos_x_vec[n0_quad8];
        double x1_quad8 = gq8s.point_pos_x_vec[n1_quad8];
        double x2_quad8 = gq8s.point_pos_x_vec[n2_quad8];
        double x3_quad8 = gq8s.point_pos_x_vec[n3_quad8];
        double x4_quad8 = gq8s.point_pos_x_vec[n4_quad8];
        double x5_quad8 = gq8s.point_pos_x_vec[n5_quad8];
        double x6_quad8 = gq8s.point_pos_x_vec[n6_quad8];
        double x7_quad8 = gq8s.point_pos_x_vec[n7_quad8];
        double x_quad8_arr[8] = {x0_quad8, x1_quad8, x2_quad8, x3_quad8, x4_quad8, x5_quad8, x6_quad8, x7_quad8};

        // get y coordinates of points around quad8 element
        double y0_quad8 = gq8s.point_pos_y_vec[n0_quad8];
        double y1_quad8 = gq8s.point_pos_y_vec[n1_quad8];
        double y2_quad8 = gq8s.point_pos_y_vec[n2_quad8];
        double y3_quad8 = gq8s.point_pos_y_vec[n3_quad8];
        double y4_quad8 = gq8s.point_pos_y_vec[n4_quad8];
        double y5_quad8 = gq8s.point_pos_y_vec[n5_quad8];
        double y6_quad8 = gq8s.point_pos_y_vec[n6_quad8];
        double y7_quad8 = gq8s.point_pos_y_vec[n7_quad8];
        double y_quad8_arr[8] = {y0_quad8, y1_quad8, y2_quad8, y3_quad8, y4_quad8, y5_quad8, y6_quad8, y7_quad8};

        // get points where the boundary is applied
        int a = v_bq8s.element_flux_pa_loc_vec[k];  // 0 to 7
        int b = v_bq8s.element_flux_pb_loc_vec[k];  // 0 to 7
        int c = v_bq8s.element_flux_pb_loc_vec[k];  // 0 to 7

        // identify boundary type
        int config_id = v_bq8s.element_flux_config_id_vec[k];
        BoundaryConfigQuad8Struct bcq8s = v_bq8s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcq8s.boundary_type_str == "neumann")
        {

            // calculate distance between a-b-c
            // APPROXIMATION ONLY. FOR REVISION.
            double xa = x_quad8_arr[a];
            double xc = x_quad8_arr[c];
            double ya = y_quad8_arr[a];
            double yc = y_quad8_arr[c];
            double dist_pa_pb_pc = sqrt((xa - xc)*(xa - xc) + (ya - yc)*(ya - yc));

            // add to b_vec
            // NOT YET VERIFIED. PLACEHOLDER ONLY.
            b_vec.coeffRef(start_id_v + n_quad8_arr[a]) += 0.5 * bcq8s.boundary_parameter_vec[0] * dist_pa_pb_pc;
            b_vec.coeffRef(start_id_v + n_quad8_arr[b]) += 0.5 * bcq8s.boundary_parameter_vec[0] * dist_pa_pb_pc;
            b_vec.coeffRef(start_id_v + n_quad8_arr[c]) += 0.5 * bcq8s.boundary_parameter_vec[0] * dist_pa_pb_pc;

        }

    }

    // iterate for each flux boundary element (continuity)
    for (int k = 0; k < p_bq4s.num_element_flux; k++)
    {

        // get id of element
        int m = p_bq4s.element_flux_id_vec[k];

        // get id of points around quad4 element
        int n0_quad4 = gq4s.element_p0_id_vec[m];
        int n1_quad4 = gq4s.element_p1_id_vec[m];
        int n2_quad4 = gq4s.element_p2_id_vec[m];
        int n3_quad4 = gq4s.element_p3_id_vec[m];
        int n_quad4_arr[4] = {n0_quad4, n1_quad4, n2_quad4, n3_quad4};

        // get x coordinates of points around element
        double x0_quad4 = gq4s.point_pos_x_vec[n0_quad4];
        double x1_quad4 = gq4s.point_pos_x_vec[n1_quad4];
        double x2_quad4 = gq4s.point_pos_x_vec[n2_quad4];
        double x3_quad4 = gq4s.point_pos_x_vec[n3_quad4];
        double x_quad4_arr[4] = {x0_quad4, x1_quad4, x2_quad4, x3_quad4};

        // get y coordinates of points around element
        double y0_quad4 = gq4s.point_pos_y_vec[n0_quad4];
        double y1_quad4 = gq4s.point_pos_y_vec[n1_quad4];
        double y2_quad4 = gq4s.point_pos_y_vec[n2_quad4];
        double y3_quad4 = gq4s.point_pos_y_vec[n3_quad4];
        double y_quad4_arr[4] = {y0_quad4, y1_quad4, y2_quad4, y3_quad4};

        // get points where the boundary is applied
        int a = p_bq4s.element_flux_pa_loc_vec[k];  // 0 to 3
        int b = p_bq4s.element_flux_pb_loc_vec[k];  // 0 to 3

        // identify boundary type
        int config_id = p_bq4s.element_flux_config_id_vec[k];
        BoundaryConfigQuad4Struct bcq4s = p_bq4s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcq4s.boundary_type_str == "neumann")
        {

            // calculate distance from point a to b
            double xa = x_quad4_arr[a];
            double xb = x_quad4_arr[b];
            double ya = y_quad4_arr[a];
            double yb = y_quad4_arr[b];
            double dist_pa_pb = sqrt((xa - xb)*(xa - xb) + (ya - yb)*(ya - yb));

            // add to b_vec
            b_vec.coeffRef(start_id_p + n_quad4_arr[a]) += -0.5 * bcq4s.boundary_parameter_vec[0] * dist_pa_pb;
            b_vec.coeffRef(start_id_p + n_quad4_arr[b]) += -0.5 * bcq4s.boundary_parameter_vec[0] * dist_pa_pb;

        }

    }

    // clear rows with value boundary elements (x-component of NSE)
    for (int k = 0; k < u_bq8s.num_element_value; k++)
    {

        // get id of element
        int m = u_bq8s.element_value_id_vec[k];

        // get id of points around element
        int n0_quad8 = gq8s.element_p0_id_vec[m];
        int n1_quad8 = gq8s.element_p1_id_vec[m];
        int n2_quad8 = gq8s.element_p2_id_vec[m];
        int n3_quad8 = gq8s.element_p3_id_vec[m];
        int n4_quad8 = gq8s.element_p4_id_vec[m];
        int n5_quad8 = gq8s.element_p5_id_vec[m];
        int n6_quad8 = gq8s.element_p6_id_vec[m];
        int n7_quad8 = gq8s.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0_quad8, n1_quad8, n2_quad8, n3_quad8, n4_quad8, n5_quad8, n6_quad8, n7_quad8};

        // get points where the boundary is applied
        int a = u_bq8s.element_value_pa_loc_vec[k];  // 0 to 7
        int b = u_bq8s.element_value_pb_loc_vec[k];  // 0 to 7
        int c = u_bq8s.element_value_pc_loc_vec[k];  // 0 to 7

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            a_mat.row(start_id_u + n_quad8_arr[a]) *= 0.;
            b_vec.coeffRef(start_id_u + n_quad8_arr[a]) = 0.;
        }
        if (b != -1)
        {
            a_mat.row(start_id_u + n_quad8_arr[b]) *= 0.;
            b_vec.coeffRef(start_id_u + n_quad8_arr[b]) = 0.;
        }
        if (c != -1)
        {
            a_mat.row(start_id_u + n_quad8_arr[c]) *= 0.;
            b_vec.coeffRef(start_id_u + n_quad8_arr[c]) = 0.;
        }

    }

    // clear rows with value boundary elements (y-component of NSE)
    for (int k = 0; k < v_bq8s.num_element_value; k++)
    {

        // get id of element
        int m = v_bq8s.element_value_id_vec[k];

        // get id of points around element
        int n0_quad8 = gq8s.element_p0_id_vec[m];
        int n1_quad8 = gq8s.element_p1_id_vec[m];
        int n2_quad8 = gq8s.element_p2_id_vec[m];
        int n3_quad8 = gq8s.element_p3_id_vec[m];
        int n4_quad8 = gq8s.element_p4_id_vec[m];
        int n5_quad8 = gq8s.element_p5_id_vec[m];
        int n6_quad8 = gq8s.element_p6_id_vec[m];
        int n7_quad8 = gq8s.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0_quad8, n1_quad8, n2_quad8, n3_quad8, n4_quad8, n5_quad8, n6_quad8, n7_quad8};

        // get points where the boundary is applied
        int a = v_bq8s.element_value_pa_loc_vec[k];  // 0 to 7
        int b = v_bq8s.element_value_pb_loc_vec[k];  // 0 to 7
        int c = v_bq8s.element_value_pc_loc_vec[k];  // 0 to 7

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            a_mat.row(start_id_v + n_quad8_arr[a]) *= 0.;
            b_vec.coeffRef(start_id_v + n_quad8_arr[a]) = 0.;
        }
        if (b != -1)
        {
            a_mat.row(start_id_v + n_quad8_arr[b]) *= 0.;
            b_vec.coeffRef(start_id_v + n_quad8_arr[b]) = 0.;
        }
        if (c != -1)
        {
            a_mat.row(start_id_v + n_quad8_arr[c]) *= 0.;
            b_vec.coeffRef(start_id_v + n_quad8_arr[c]) = 0.;
        }

    }

    // clear rows with value boundary elements (continuity)
    for (int k = 0; k < p_bq4s.num_element_value; k++)
    {

        // get id of element
        int m = p_bq4s.element_value_id_vec[k];

        // get id of points around element
        int n0_quad4 = gq8s.element_p0_id_vec[m];
        int n1_quad4 = gq8s.element_p1_id_vec[m];
        int n2_quad4 = gq8s.element_p2_id_vec[m];
        int n3_quad4 = gq8s.element_p3_id_vec[m];
        int n_quad4_arr[4] = {n0_quad4, n1_quad4, n2_quad4, n3_quad4};

        // get points where the boundary is applied
        int a = p_bq4s.element_value_pa_loc_vec[k];  // 0 to 3
        int b = p_bq4s.element_value_pb_loc_vec[k];  // 0 to 3

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            a_mat.row(start_id_p + n_quad4_arr[a]) *= 0.;
            b_vec.coeffRef(start_id_p + n_quad4_arr[a]) = 0.;
        }
        if (b != -1)
        {
            a_mat.row(start_id_p + n_quad4_arr[b]) *= 0.;
            b_vec.coeffRef(start_id_p + n_quad4_arr[b]) = 0.;
        }

    }

    // iterate for each value boundary element (x-component of NSE)
    for (int k = 0; k < u_bq8s.num_element_value; k++)
    {

        // get id of element
        int m = u_bq8s.element_value_id_vec[k];

        // get id of points around quad8 element
        int n0_quad8 = gq8s.element_p0_id_vec[m];
        int n1_quad8 = gq8s.element_p1_id_vec[m];
        int n2_quad8 = gq8s.element_p2_id_vec[m];
        int n3_quad8 = gq8s.element_p3_id_vec[m];
        int n4_quad8 = gq8s.element_p4_id_vec[m];
        int n5_quad8 = gq8s.element_p5_id_vec[m];
        int n6_quad8 = gq8s.element_p6_id_vec[m];
        int n7_quad8 = gq8s.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0_quad8, n1_quad8, n2_quad8, n3_quad8, n4_quad8, n5_quad8, n6_quad8, n7_quad8};

        // get points where the boundary is applied
        int a = u_bq8s.element_value_pa_loc_vec[k];  // 0 to 7
        int b = u_bq8s.element_value_pb_loc_vec[k];  // 0 to 7
        int c = u_bq8s.element_value_pc_loc_vec[k];  // 0 to 7

        // identify boundary type
        int config_id = u_bq8s.element_value_config_id_vec[k];
        BoundaryConfigQuad8Struct bcq8s = u_bq8s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcq8s.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            if (a != -1)
            {
                a_mat.coeffRef(start_id_u + n_quad8_arr[a], start_id_u + n_quad8_arr[a]) += 1.;
                b_vec.coeffRef(start_id_u + n_quad8_arr[a]) += bcq8s.boundary_parameter_vec[0];
            }
            if (b != -1)
            {
                a_mat.coeffRef(start_id_u + n_quad8_arr[b], start_id_u + n_quad8_arr[b]) += 1.;
                b_vec.coeffRef(start_id_u + n_quad8_arr[b]) += bcq8s.boundary_parameter_vec[0];
            }
            if (c != -1)
            {
                a_mat.coeffRef(start_id_u + n_quad8_arr[c], start_id_u + n_quad8_arr[c]) += 1.;
                b_vec.coeffRef(start_id_u + n_quad8_arr[c]) += bcq8s.boundary_parameter_vec[0];
            }

        }

    }

    // iterate for each value boundary element (y-component of NSE)
    for (int k = 0; k < v_bq8s.num_element_value; k++)
    {

        // get id of element
        int m = v_bq8s.element_value_id_vec[k];

        // get id of points around quad8 element
        int n0_quad8 = gq8s.element_p0_id_vec[m];
        int n1_quad8 = gq8s.element_p1_id_vec[m];
        int n2_quad8 = gq8s.element_p2_id_vec[m];
        int n3_quad8 = gq8s.element_p3_id_vec[m];
        int n4_quad8 = gq8s.element_p4_id_vec[m];
        int n5_quad8 = gq8s.element_p5_id_vec[m];
        int n6_quad8 = gq8s.element_p6_id_vec[m];
        int n7_quad8 = gq8s.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0_quad8, n1_quad8, n2_quad8, n3_quad8, n4_quad8, n5_quad8, n6_quad8, n7_quad8};

        // get points where the boundary is applied
        int a = v_bq8s.element_value_pa_loc_vec[k];  // 0 to 7
        int b = v_bq8s.element_value_pb_loc_vec[k];  // 0 to 7
        int c = v_bq8s.element_value_pc_loc_vec[k];  // 0 to 7

        // identify boundary type
        int config_id = v_bq8s.element_value_config_id_vec[k];
        BoundaryConfigQuad8Struct bcq8s = v_bq8s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcq8s.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            if (a != -1)
            {
                a_mat.coeffRef(start_id_v + n_quad8_arr[a], start_id_v + n_quad8_arr[a]) += 1.;
                b_vec.coeffRef(start_id_v + n_quad8_arr[a]) += bcq8s.boundary_parameter_vec[0];
            }
            if (b != -1)
            {
                a_mat.coeffRef(start_id_v + n_quad8_arr[b], start_id_v + n_quad8_arr[b]) += 1.;
                b_vec.coeffRef(start_id_v + n_quad8_arr[b]) += bcq8s.boundary_parameter_vec[0];
            }
            if (c != -1)
            {
                a_mat.coeffRef(start_id_v + n_quad8_arr[c], start_id_v + n_quad8_arr[c]) += 1.;
                b_vec.coeffRef(start_id_v + n_quad8_arr[c]) += bcq8s.boundary_parameter_vec[0];
            }

        }

    }

    // iterate for each value boundary element (continuity)
    for (int k = 0; k < p_bq4s.num_element_value; k++)
    {

        // get id of element
        int m = p_bq4s.element_value_id_vec[k];

        // get id of points around quad8 element
        int n0_quad4 = gq4s.element_p0_id_vec[m];
        int n1_quad4 = gq4s.element_p1_id_vec[m];
        int n2_quad4 = gq4s.element_p2_id_vec[m];
        int n3_quad4 = gq4s.element_p3_id_vec[m];
        int n_quad4_arr[4] = {n0_quad4, n1_quad4, n2_quad4, n3_quad4};

        // get points where the boundary is applied
        int a = p_bq4s.element_value_pa_loc_vec[k];  // 0 to 3
        int b = p_bq4s.element_value_pb_loc_vec[k];  // 0 to 3

        // identify boundary type
        int config_id = p_bq4s.element_value_config_id_vec[k];
        BoundaryConfigQuad4Struct bcq4s = p_bq4s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcq4s.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            if (a != -1)
            {
                a_mat.coeffRef(start_id_p + n_quad4_arr[a], start_id_p + n_quad4_arr[a]) += 1.;
                b_vec.coeffRef(start_id_p + n_quad4_arr[a]) += bcq4s.boundary_parameter_vec[0];
            }
            if (b != -1)
            {
                a_mat.coeffRef(start_id_p + n_quad4_arr[b], start_id_p + n_quad4_arr[b]) += 1.;
                b_vec.coeffRef(start_id_p + n_quad4_arr[b]) += bcq4s.boundary_parameter_vec[0];
            }

        }

    }

}

#endif
