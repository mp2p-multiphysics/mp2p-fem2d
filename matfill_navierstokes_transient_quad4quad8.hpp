#ifndef MATFILL_NAVIERSTOKES_STEADY_QUAD4QUAD8
#define MATFILL_NAVIERSTOKES_STEADY_QUAD4QUAD8
#include <vector>
#include "Eigen/Eigen"
#include "helper_integrate_quad4quad8.hpp"
#include "helper_struct.hpp"

void matfill_navierstokes_transient_quad4quad8
(
    GridQuad8Data u_gq8d, GridQuad8Data v_gq8d, GridQuad4Data p_gq4d,
    BoundaryQuad8Data u_bq8d, BoundaryQuad8Data v_bq8d, BoundaryQuad4Data p_bq4d,
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &x_it_vec, Eigen::VectorXd &x_ts_vec, Eigen::VectorXd &b_vec,
    double dt, double density, double viscosity, double gravity_x, double gravity_y,
    int start_id_u, int start_id_v, int start_id_p
)
{

    // initialize a_mat for b vector
    int num_point_all = u_gq8d.num_point + v_gq8d.num_point + p_gq4d.num_point;
    Eigen::SparseMatrix<double> a_mat_helper(num_point_all, num_point_all);

    // fill up matrix with x-component of NSE
    for (int m = 0; m < u_gq8d.num_element; m++)
    {

        // get id of points around quad8 element
        int n0_quad8 = u_gq8d.element_p0_id_vec[m];
        int n1_quad8 = u_gq8d.element_p1_id_vec[m];
        int n2_quad8 = u_gq8d.element_p2_id_vec[m];
        int n3_quad8 = u_gq8d.element_p3_id_vec[m];
        int n4_quad8 = u_gq8d.element_p4_id_vec[m];
        int n5_quad8 = u_gq8d.element_p5_id_vec[m];
        int n6_quad8 = u_gq8d.element_p6_id_vec[m];
        int n7_quad8 = u_gq8d.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0_quad8, n1_quad8, n2_quad8, n3_quad8, n4_quad8, n5_quad8, n6_quad8, n7_quad8};

        // get id of points around quad4 element
        int n0_quad4 = p_gq4d.element_p0_id_vec[m];
        int n1_quad4 = p_gq4d.element_p1_id_vec[m];
        int n2_quad4 = p_gq4d.element_p2_id_vec[m];
        int n3_quad4 = p_gq4d.element_p3_id_vec[m];
        int n_quad4_arr[4] = {n0_quad4, n1_quad4, n2_quad4, n3_quad4};

        // get x coordinates of points around element
        double x0 = u_gq8d.point_pos_x_vec[n0_quad8];
        double x1 = u_gq8d.point_pos_x_vec[n1_quad8];
        double x2 = u_gq8d.point_pos_x_vec[n2_quad8];
        double x3 = u_gq8d.point_pos_x_vec[n3_quad8];
        double x4 = u_gq8d.point_pos_x_vec[n4_quad8];
        double x5 = u_gq8d.point_pos_x_vec[n5_quad8];
        double x6 = u_gq8d.point_pos_x_vec[n6_quad8];
        double x7 = u_gq8d.point_pos_x_vec[n7_quad8];
        double x_quad8_arr[8] = {x0, x1, x2, x3, x4, x5, x6, x7};

        // get y coordinates of points around element
        double y0 = u_gq8d.point_pos_y_vec[n0_quad8];
        double y1 = u_gq8d.point_pos_y_vec[n1_quad8];
        double y2 = u_gq8d.point_pos_y_vec[n2_quad8];
        double y3 = u_gq8d.point_pos_y_vec[n3_quad8];
        double y4 = u_gq8d.point_pos_y_vec[n4_quad8];
        double y5 = u_gq8d.point_pos_y_vec[n5_quad8];
        double y6 = u_gq8d.point_pos_y_vec[n6_quad8];
        double y7 = u_gq8d.point_pos_y_vec[n7_quad8];
        double y_quad8_arr[8] = {y0, y1, y2, y3, y4, y5, y6, y7};

        // get u velocity values at points around element (last iteration)
        double u0_it = x_it_vec[start_id_u + n0_quad8];
        double u1_it = x_it_vec[start_id_u + n1_quad8];
        double u2_it = x_it_vec[start_id_u + n2_quad8];
        double u3_it = x_it_vec[start_id_u + n3_quad8];
        double u4_it = x_it_vec[start_id_u + n4_quad8];
        double u5_it = x_it_vec[start_id_u + n5_quad8];
        double u6_it = x_it_vec[start_id_u + n6_quad8];
        double u7_it = x_it_vec[start_id_u + n7_quad8];
        double u_it_quad8_arr[8] = {u0_it, u1_it, u2_it, u3_it, u4_it, u5_it, u6_it, u7_it};

        // get v velocity values at points around element (last iteration)
        double v0_it = x_it_vec[start_id_v + n0_quad8];
        double v1_it = x_it_vec[start_id_v + n1_quad8];
        double v2_it = x_it_vec[start_id_v + n2_quad8];
        double v3_it = x_it_vec[start_id_v + n3_quad8];
        double v4_it = x_it_vec[start_id_v + n4_quad8];
        double v5_it = x_it_vec[start_id_v + n5_quad8];
        double v6_it = x_it_vec[start_id_v + n6_quad8];
        double v7_it = x_it_vec[start_id_v + n7_quad8];
        double v_it_quad8_arr[8] = {v0_it, v1_it, v2_it, v3_it, v4_it, v5_it, v6_it, v7_it};

        // get u velocity values at points around element (last timestep)
        double u0_ts = x_ts_vec[start_id_u + n0_quad8];
        double u1_ts = x_ts_vec[start_id_u + n1_quad8];
        double u2_ts = x_ts_vec[start_id_u + n2_quad8];
        double u3_ts = x_ts_vec[start_id_u + n3_quad8];
        double u4_ts = x_ts_vec[start_id_u + n4_quad8];
        double u5_ts = x_ts_vec[start_id_u + n5_quad8];
        double u6_ts = x_ts_vec[start_id_u + n6_quad8];
        double u7_ts = x_ts_vec[start_id_u + n7_quad8];
        double u_ts_quad8_arr[8] = {u0_ts, u1_ts, u2_ts, u3_ts, u4_ts, u5_ts, u6_ts, u7_ts};

        // get v velocity values at points around element (last timestep)
        double v0_ts = x_ts_vec[start_id_v + n0_quad8];
        double v1_ts = x_ts_vec[start_id_v + n1_quad8];
        double v2_ts = x_ts_vec[start_id_v + n2_quad8];
        double v3_ts = x_ts_vec[start_id_v + n3_quad8];
        double v4_ts = x_ts_vec[start_id_v + n4_quad8];
        double v5_ts = x_ts_vec[start_id_v + n5_quad8];
        double v6_ts = x_ts_vec[start_id_v + n6_quad8];
        double v7_ts = x_ts_vec[start_id_v + n7_quad8];
        double v_ts_quad8_arr[8] = {v0_ts, v1_ts, v2_ts, v3_ts, v4_ts, v5_ts, v6_ts, v7_ts};

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable [u0 ... uN v0 ... vN p0 ... pN]

        // coefficients of u
        for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            
            // a_mat
            a_mat.coeffRef(start_id_u + n_quad8_arr[i], start_id_u + n_quad8_arr[j]) += (
                  integrate_constant_Mi_Mj_quad8(density/dt, i, j, x_quad8_arr, y_quad8_arr)
                + integrate_constant_Mi_quad8_u_derivative_Mj_quad8_x(0.5*density, i, j, u_it_quad8_arr, x_quad8_arr, y_quad8_arr)
                + integrate_constant_Mi_quad8_u_derivative_Mj_quad8_y(0.5*density, i, j, v_it_quad8_arr, x_quad8_arr, y_quad8_arr)
                + integrate_constant_div_Mi_quad8_dot_div_Mj_quad8(0.5*viscosity, i, j, x_quad8_arr, y_quad8_arr)
            );

            // a_mat_helper
            a_mat_helper.coeffRef(start_id_u + n_quad8_arr[i], start_id_u + n_quad8_arr[j]) += (
                  integrate_constant_Mi_Mj_quad8(density/dt, i, j, x_quad8_arr, y_quad8_arr)
                - integrate_constant_Mi_quad8_u_derivative_Mj_quad8_x(0.5*density, i, j, u_ts_quad8_arr, x_quad8_arr, y_quad8_arr)
                - integrate_constant_Mi_quad8_u_derivative_Mj_quad8_y(0.5*density, i, j, v_ts_quad8_arr, x_quad8_arr, y_quad8_arr)
                - integrate_constant_div_Mi_quad8_dot_div_Mj_quad8(0.5*viscosity, i, j, x_quad8_arr, y_quad8_arr)
            );
            
        }}

        // coefficients of p
        for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 4; j++) {
            a_mat.coeffRef(start_id_u + n_quad8_arr[i], start_id_p + n_quad4_arr[j]) += integrate_constant_Mi_quad8_derivative_Nj_quad4_x(0.5, i, j, x_quad8_arr, y_quad8_arr);
            a_mat_helper.coeffRef(start_id_u + n_quad8_arr[i], start_id_p + n_quad4_arr[j]) += -integrate_constant_Mi_quad8_derivative_Nj_quad4_x(0.5, i, j, x_quad8_arr, y_quad8_arr);
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 8; i++) {
            b_vec.coeffRef(start_id_u + n_quad8_arr[i]) += integrate_constant_Mi_quad8(density*gravity_x, i, x_quad8_arr, y_quad8_arr);
        }

    }

    // fill up matrix with y-component of NSE
    for (int m = 0; m < v_gq8d.num_element; m++)
    {

        // get id of points around quad8 element
        int n0_quad8 = u_gq8d.element_p0_id_vec[m];
        int n1_quad8 = u_gq8d.element_p1_id_vec[m];
        int n2_quad8 = u_gq8d.element_p2_id_vec[m];
        int n3_quad8 = u_gq8d.element_p3_id_vec[m];
        int n4_quad8 = u_gq8d.element_p4_id_vec[m];
        int n5_quad8 = u_gq8d.element_p5_id_vec[m];
        int n6_quad8 = u_gq8d.element_p6_id_vec[m];
        int n7_quad8 = u_gq8d.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0_quad8, n1_quad8, n2_quad8, n3_quad8, n4_quad8, n5_quad8, n6_quad8, n7_quad8};

        // get id of points around quad4 element
        int n0_quad4 = p_gq4d.element_p0_id_vec[m];
        int n1_quad4 = p_gq4d.element_p1_id_vec[m];
        int n2_quad4 = p_gq4d.element_p2_id_vec[m];
        int n3_quad4 = p_gq4d.element_p3_id_vec[m];
        int n_quad4_arr[4] = {n0_quad4, n1_quad4, n2_quad4, n3_quad4};

        // get x coordinates of points around element
        double x0 = u_gq8d.point_pos_x_vec[n0_quad8];
        double x1 = u_gq8d.point_pos_x_vec[n1_quad8];
        double x2 = u_gq8d.point_pos_x_vec[n2_quad8];
        double x3 = u_gq8d.point_pos_x_vec[n3_quad8];
        double x4 = u_gq8d.point_pos_x_vec[n4_quad8];
        double x5 = u_gq8d.point_pos_x_vec[n5_quad8];
        double x6 = u_gq8d.point_pos_x_vec[n6_quad8];
        double x7 = u_gq8d.point_pos_x_vec[n7_quad8];
        double x_quad8_arr[8] = {x0, x1, x2, x3, x4, x5, x6, x7};

        // get y coordinates of points around element
        double y0 = u_gq8d.point_pos_y_vec[n0_quad8];
        double y1 = u_gq8d.point_pos_y_vec[n1_quad8];
        double y2 = u_gq8d.point_pos_y_vec[n2_quad8];
        double y3 = u_gq8d.point_pos_y_vec[n3_quad8];
        double y4 = u_gq8d.point_pos_y_vec[n4_quad8];
        double y5 = u_gq8d.point_pos_y_vec[n5_quad8];
        double y6 = u_gq8d.point_pos_y_vec[n6_quad8];
        double y7 = u_gq8d.point_pos_y_vec[n7_quad8];
        double y_quad8_arr[8] = {y0, y1, y2, y3, y4, y5, y6, y7};

        // get u velocity values at points around element (last iteration)
        double u0_it = x_it_vec[start_id_u + n0_quad8];
        double u1_it = x_it_vec[start_id_u + n1_quad8];
        double u2_it = x_it_vec[start_id_u + n2_quad8];
        double u3_it = x_it_vec[start_id_u + n3_quad8];
        double u4_it = x_it_vec[start_id_u + n4_quad8];
        double u5_it = x_it_vec[start_id_u + n5_quad8];
        double u6_it = x_it_vec[start_id_u + n6_quad8];
        double u7_it = x_it_vec[start_id_u + n7_quad8];
        double u_it_quad8_arr[8] = {u0_it, u1_it, u2_it, u3_it, u4_it, u5_it, u6_it, u7_it};

        // get v velocity values at points around element (last iteration)
        double v0_it = x_it_vec[start_id_v + n0_quad8];
        double v1_it = x_it_vec[start_id_v + n1_quad8];
        double v2_it = x_it_vec[start_id_v + n2_quad8];
        double v3_it = x_it_vec[start_id_v + n3_quad8];
        double v4_it = x_it_vec[start_id_v + n4_quad8];
        double v5_it = x_it_vec[start_id_v + n5_quad8];
        double v6_it = x_it_vec[start_id_v + n6_quad8];
        double v7_it = x_it_vec[start_id_v + n7_quad8];
        double v_it_quad8_arr[8] = {v0_it, v1_it, v2_it, v3_it, v4_it, v5_it, v6_it, v7_it};

        // get u velocity values at points around element (last timestep)
        double u0_ts = x_ts_vec[start_id_u + n0_quad8];
        double u1_ts = x_ts_vec[start_id_u + n1_quad8];
        double u2_ts = x_ts_vec[start_id_u + n2_quad8];
        double u3_ts = x_ts_vec[start_id_u + n3_quad8];
        double u4_ts = x_ts_vec[start_id_u + n4_quad8];
        double u5_ts = x_ts_vec[start_id_u + n5_quad8];
        double u6_ts = x_ts_vec[start_id_u + n6_quad8];
        double u7_ts = x_ts_vec[start_id_u + n7_quad8];
        double u_ts_quad8_arr[8] = {u0_ts, u1_ts, u2_ts, u3_ts, u4_ts, u5_ts, u6_ts, u7_ts};

        // get v velocity values at points around element (last timestep)
        double v0_ts = x_ts_vec[start_id_v + n0_quad8];
        double v1_ts = x_ts_vec[start_id_v + n1_quad8];
        double v2_ts = x_ts_vec[start_id_v + n2_quad8];
        double v3_ts = x_ts_vec[start_id_v + n3_quad8];
        double v4_ts = x_ts_vec[start_id_v + n4_quad8];
        double v5_ts = x_ts_vec[start_id_v + n5_quad8];
        double v6_ts = x_ts_vec[start_id_v + n6_quad8];
        double v7_ts = x_ts_vec[start_id_v + n7_quad8];
        double v_ts_quad8_arr[8] = {v0_ts, v1_ts, v2_ts, v3_ts, v4_ts, v5_ts, v6_ts, v7_ts};

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable [u0 ... uN v0 ... vN p0 ... pN]

        // coefficients of v
        for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            
            // a_mat
            a_mat.coeffRef(start_id_v + n_quad8_arr[i], start_id_v + n_quad8_arr[j]) += (
                  integrate_constant_Mi_Mj_quad8(density/dt, i, j, x_quad8_arr, y_quad8_arr)
                + integrate_constant_Mi_quad8_u_derivative_Mj_quad8_x(0.5*density, i, j, u_it_quad8_arr, x_quad8_arr, y_quad8_arr)
                + integrate_constant_Mi_quad8_u_derivative_Mj_quad8_y(0.5*density, i, j, v_it_quad8_arr, x_quad8_arr, y_quad8_arr)
                + integrate_constant_div_Mi_quad8_dot_div_Mj_quad8(0.5*viscosity, i, j, x_quad8_arr, y_quad8_arr)
            );

            // a_mat_helper
            a_mat_helper.coeffRef(start_id_v + n_quad8_arr[i], start_id_v + n_quad8_arr[j]) += (
                  integrate_constant_Mi_Mj_quad8(density/dt, i, j, x_quad8_arr, y_quad8_arr)
                - integrate_constant_Mi_quad8_u_derivative_Mj_quad8_x(0.5*density, i, j, u_ts_quad8_arr, x_quad8_arr, y_quad8_arr)
                - integrate_constant_Mi_quad8_u_derivative_Mj_quad8_y(0.5*density, i, j, v_ts_quad8_arr, x_quad8_arr, y_quad8_arr)
                - integrate_constant_div_Mi_quad8_dot_div_Mj_quad8(0.5*viscosity, i, j, x_quad8_arr, y_quad8_arr)
            );
            
        }}

        // coefficients of p
        for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 4; j++) {
            a_mat.coeffRef(start_id_v + n_quad8_arr[i], start_id_p + n_quad4_arr[j]) += integrate_constant_Mi_quad8_derivative_Nj_quad4_y(0.5, i, j, x_quad8_arr, y_quad8_arr);
            a_mat_helper.coeffRef(start_id_v + n_quad8_arr[i], start_id_p + n_quad4_arr[j]) += -integrate_constant_Mi_quad8_derivative_Nj_quad4_y(0.5, i, j, x_quad8_arr, y_quad8_arr);
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 8; i++) {
            b_vec.coeffRef(start_id_v + n_quad8_arr[i]) += integrate_constant_Mi_quad8(density*gravity_y, i, x_quad8_arr, y_quad8_arr);
        }

    }

    // fill up matrix with continuity equation
    for (int m = 0; m < p_gq4d.num_element; m++)
    {

        // get id of points around quad8 element
        int n0_quad8 = u_gq8d.element_p0_id_vec[m];
        int n1_quad8 = u_gq8d.element_p1_id_vec[m];
        int n2_quad8 = u_gq8d.element_p2_id_vec[m];
        int n3_quad8 = u_gq8d.element_p3_id_vec[m];
        int n4_quad8 = u_gq8d.element_p4_id_vec[m];
        int n5_quad8 = u_gq8d.element_p5_id_vec[m];
        int n6_quad8 = u_gq8d.element_p6_id_vec[m];
        int n7_quad8 = u_gq8d.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0_quad8, n1_quad8, n2_quad8, n3_quad8, n4_quad8, n5_quad8, n6_quad8, n7_quad8};

        // get id of points around quad4 element
        int n0_quad4 = p_gq4d.element_p0_id_vec[m];
        int n1_quad4 = p_gq4d.element_p1_id_vec[m];
        int n2_quad4 = p_gq4d.element_p2_id_vec[m];
        int n3_quad4 = p_gq4d.element_p3_id_vec[m];
        int n_quad4_arr[4] = {n0_quad4, n1_quad4, n2_quad4, n3_quad4};

        // get x coordinates of points around element
        double x0 = u_gq8d.point_pos_x_vec[n0_quad8];
        double x1 = u_gq8d.point_pos_x_vec[n1_quad8];
        double x2 = u_gq8d.point_pos_x_vec[n2_quad8];
        double x3 = u_gq8d.point_pos_x_vec[n3_quad8];
        double x4 = u_gq8d.point_pos_x_vec[n4_quad8];
        double x5 = u_gq8d.point_pos_x_vec[n5_quad8];
        double x6 = u_gq8d.point_pos_x_vec[n6_quad8];
        double x7 = u_gq8d.point_pos_x_vec[n7_quad8];
        double x_quad8_arr[8] = {x0, x1, x2, x3, x4, x5, x6, x7};

        // get y coordinates of points around element
        double y0 = u_gq8d.point_pos_y_vec[n0_quad8];
        double y1 = u_gq8d.point_pos_y_vec[n1_quad8];
        double y2 = u_gq8d.point_pos_y_vec[n2_quad8];
        double y3 = u_gq8d.point_pos_y_vec[n3_quad8];
        double y4 = u_gq8d.point_pos_y_vec[n4_quad8];
        double y5 = u_gq8d.point_pos_y_vec[n5_quad8];
        double y6 = u_gq8d.point_pos_y_vec[n6_quad8];
        double y7 = u_gq8d.point_pos_y_vec[n7_quad8];
        double y_quad8_arr[8] = {y0, y1, y2, y3, y4, y5, y6, y7};

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable [u0 ... uN v0 ... vN p0 ... pN]

        // coefficients of u
        for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            a_mat.coeffRef(start_id_p + n_quad4_arr[i], start_id_u + n_quad8_arr[j]) += integrate_constant_Ni_quad4_derivative_Mj_quad8_x(1., i, j, x_quad8_arr, y_quad8_arr);
        }}

        // coefficients of v
        for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            a_mat.coeffRef(start_id_p + n_quad4_arr[i], start_id_v + n_quad8_arr[j]) += integrate_constant_Ni_quad4_derivative_Mj_quad8_y(1., i, j, x_quad8_arr, y_quad8_arr);
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 4; i++) {
            b_vec.coeffRef(start_id_p + n_quad4_arr[i]) += 0.;
        }

    }

    // iterate for each flux boundary element (x-component of NSE)
    for (int k = 0; u_bq8d.num_element_flux; k++)
    {

        // get id of element
        int m = u_bq8d.element_value_id_vec[k];

        // get id of points around quad8 element
        int n0 = u_gq8d.element_p0_id_vec[m];
        int n1 = u_gq8d.element_p1_id_vec[m];
        int n2 = u_gq8d.element_p2_id_vec[m];
        int n3 = u_gq8d.element_p3_id_vec[m];
        int n4 = u_gq8d.element_p4_id_vec[m];
        int n5 = u_gq8d.element_p5_id_vec[m];
        int n6 = u_gq8d.element_p6_id_vec[m];
        int n7 = u_gq8d.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0, n1, n2, n3, n4, n5, n6, n7};

        // get x coordinates of points around element
        double x0 = u_gq8d.point_pos_x_vec[n0];
        double x1 = u_gq8d.point_pos_x_vec[n1];
        double x2 = u_gq8d.point_pos_x_vec[n2];
        double x3 = u_gq8d.point_pos_x_vec[n3];
        double x4 = u_gq8d.point_pos_x_vec[n4];
        double x5 = u_gq8d.point_pos_x_vec[n5];
        double x6 = u_gq8d.point_pos_x_vec[n6];
        double x7 = u_gq8d.point_pos_x_vec[n7];
        double x_quad8_arr[8] = {x0, x1, x2, x3, x4, x5, x6, x7};

        // get y coordinates of points around element
        double y0 = u_gq8d.point_pos_y_vec[n0];
        double y1 = u_gq8d.point_pos_y_vec[n1];
        double y2 = u_gq8d.point_pos_y_vec[n2];
        double y3 = u_gq8d.point_pos_y_vec[n3];
        double y4 = u_gq8d.point_pos_y_vec[n4];
        double y5 = u_gq8d.point_pos_y_vec[n5];
        double y6 = u_gq8d.point_pos_y_vec[n6];
        double y7 = u_gq8d.point_pos_y_vec[n7];
        double y_quad8_arr[8] = {y0, y1, y2, y3, y4, y5, y6, y7};

        // get points where the boundary is applied
        int a = u_bq8d.element_flux_pa_loc_vec[k];  // 0 to 7
        int b = u_bq8d.element_flux_pb_loc_vec[k];  // 0 to 7
        int c = u_bq8d.element_flux_pb_loc_vec[k];  // 0 to 7

        // identify boundary type
        int config_id = u_bq8d.element_flux_config_id_vec[k];
        BoundaryConfigData bcd = u_bq8d.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcd.boundary_type_str == "neumann")
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
            b_vec.coeffRef(start_id_u + n_quad8_arr[a]) += -0.5 * bcd.boundary_parameter_vec[0] * dist_pa_pb_pc;
            b_vec.coeffRef(start_id_u + n_quad8_arr[b]) += -0.5 * bcd.boundary_parameter_vec[0] * dist_pa_pb_pc;
            b_vec.coeffRef(start_id_u + n_quad8_arr[c]) += -0.5 * bcd.boundary_parameter_vec[0] * dist_pa_pb_pc;

        }

    }

    // add contents of helper matrix
    b_vec += a_mat_helper*x_ts_vec;

    // iterate for each flux boundary element (y-component of NSE)
    for (int k = 0; v_bq8d.num_element_flux; k++)
    {
        
        // get id of element
        int m = v_bq8d.element_value_id_vec[k];

        // get id of points around quad 8 element
        int n0 = u_gq8d.element_p0_id_vec[m];
        int n1 = u_gq8d.element_p1_id_vec[m];
        int n2 = u_gq8d.element_p2_id_vec[m];
        int n3 = u_gq8d.element_p3_id_vec[m];
        int n4 = u_gq8d.element_p4_id_vec[m];
        int n5 = u_gq8d.element_p5_id_vec[m];
        int n6 = u_gq8d.element_p6_id_vec[m];
        int n7 = u_gq8d.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0, n1, n2, n3, n4, n5, n6, n7};

        // get x coordinates of points around element
        double x0 = u_gq8d.point_pos_x_vec[n0];
        double x1 = u_gq8d.point_pos_x_vec[n1];
        double x2 = u_gq8d.point_pos_x_vec[n2];
        double x3 = u_gq8d.point_pos_x_vec[n3];
        double x4 = u_gq8d.point_pos_x_vec[n4];
        double x5 = u_gq8d.point_pos_x_vec[n5];
        double x6 = u_gq8d.point_pos_x_vec[n6];
        double x7 = u_gq8d.point_pos_x_vec[n7];
        double x_quad8_arr[8] = {x0, x1, x2, x3, x4, x5, x6, x7};

        // get y coordinates of points around element
        double y0 = u_gq8d.point_pos_y_vec[n0];
        double y1 = u_gq8d.point_pos_y_vec[n1];
        double y2 = u_gq8d.point_pos_y_vec[n2];
        double y3 = u_gq8d.point_pos_y_vec[n3];
        double y4 = u_gq8d.point_pos_y_vec[n4];
        double y5 = u_gq8d.point_pos_y_vec[n5];
        double y6 = u_gq8d.point_pos_y_vec[n6];
        double y7 = u_gq8d.point_pos_y_vec[n7];
        double y_quad8_arr[8] = {y0, y1, y2, y3, y4, y5, y6, y7};

        // get points where the boundary is applied
        int a = v_bq8d.element_flux_pa_loc_vec[k];  // 0 to 7
        int b = v_bq8d.element_flux_pb_loc_vec[k];  // 0 to 7
        int c = v_bq8d.element_flux_pb_loc_vec[k];  // 0 to 7

        // identify boundary type
        int config_id = v_bq8d.element_flux_config_id_vec[k];
        BoundaryConfigData bcd = v_bq8d.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcd.boundary_type_str == "neumann")
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
            b_vec.coeffRef(start_id_v + n_quad8_arr[a]) += -0.5 * bcd.boundary_parameter_vec[0] * dist_pa_pb_pc;
            b_vec.coeffRef(start_id_v + n_quad8_arr[b]) += -0.5 * bcd.boundary_parameter_vec[0] * dist_pa_pb_pc;
            b_vec.coeffRef(start_id_v + n_quad8_arr[c]) += -0.5 * bcd.boundary_parameter_vec[0] * dist_pa_pb_pc;

        }

    }

    // iterate for each flux boundary element (continuity)
    for (int k = 0; p_bq4d.num_element_flux; k++)
    {

        // get id of element
        int m = p_bq4d.element_value_id_vec[k];

        // get id of points around quad4 element
        int n0 = p_gq4d.element_p0_id_vec[m];
        int n1 = p_gq4d.element_p1_id_vec[m];
        int n2 = p_gq4d.element_p2_id_vec[m];
        int n3 = p_gq4d.element_p3_id_vec[m];
        int n_quad4_arr[4] = {n0, n1, n2, n3};

        // get x coordinates of points around element
        double x0 = p_gq4d.point_pos_x_vec[n0];
        double x1 = p_gq4d.point_pos_x_vec[n1];
        double x2 = p_gq4d.point_pos_x_vec[n2];
        double x3 = p_gq4d.point_pos_x_vec[n3];
        double x_quad4_arr[4] = {x0, x1, x2, x3};

        // get y coordinates of points around element
        double y0 = p_gq4d.point_pos_y_vec[n0];
        double y1 = p_gq4d.point_pos_y_vec[n1];
        double y2 = p_gq4d.point_pos_y_vec[n2];
        double y3 = p_gq4d.point_pos_y_vec[n3];
        double y_quad4_arr[4] = {y0, y1, y2, y3};

        // get points where the boundary is applied
        int a = p_bq4d.element_flux_pa_loc_vec[k];  // 0 to 3
        int b = p_bq4d.element_flux_pb_loc_vec[k];  // 0 to 3

        // identify boundary type
        int config_id = p_bq4d.element_flux_config_id_vec[k];
        BoundaryConfigData bcd = p_bq4d.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcd.boundary_type_str == "neumann")
        {

            // calculate distance from point a to b
            double xa = x_quad4_arr[a];
            double xb = x_quad4_arr[b];
            double ya = y_quad4_arr[a];
            double yb = y_quad4_arr[b];
            double dist_pa_pb = sqrt((xa - xb)*(xa - xb) + (ya - yb)*(ya - yb));

            // add to b_vec
            b_vec.coeffRef(start_id_p + n_quad4_arr[a]) += -0.5 * bcd.boundary_parameter_vec[0] * dist_pa_pb;
            b_vec.coeffRef(start_id_p + n_quad4_arr[b]) += -0.5 * bcd.boundary_parameter_vec[0] * dist_pa_pb;

        }

    }

    // clear rows with value boundary elements (x-component of NSE)
    for (int k = 0; k < u_bq8d.num_element_value; k++)
    {

        // get id of element
        int m = u_bq8d.element_value_id_vec[k];

        // get id of points around quad8 element
        int n0 = u_gq8d.element_p0_id_vec[m];
        int n1 = u_gq8d.element_p1_id_vec[m];
        int n2 = u_gq8d.element_p2_id_vec[m];
        int n3 = u_gq8d.element_p3_id_vec[m];
        int n4 = u_gq8d.element_p4_id_vec[m];
        int n5 = u_gq8d.element_p5_id_vec[m];
        int n6 = u_gq8d.element_p6_id_vec[m];
        int n7 = u_gq8d.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0, n1, n2, n3, n4, n5, n6, n7};

        // get points where the boundary is applied
        int loc = u_bq8d.element_value_loc_vec[k];  // 0 to 7

        // erase entire row
        a_mat.row(start_id_u + n_quad8_arr[loc]) *= 0.;
        b_vec.coeffRef(start_id_u + n_quad8_arr[loc]) = 0.;

    }

    // clear rows with value boundary elements (y-component of NSE)
    for (int k = 0; k < v_bq8d.num_element_value; k++)
    {

        // get id of element
        int m = v_bq8d.element_value_id_vec[k];

        // get id of points around quad8 element
        int n0 = u_gq8d.element_p0_id_vec[m];
        int n1 = u_gq8d.element_p1_id_vec[m];
        int n2 = u_gq8d.element_p2_id_vec[m];
        int n3 = u_gq8d.element_p3_id_vec[m];
        int n4 = u_gq8d.element_p4_id_vec[m];
        int n5 = u_gq8d.element_p5_id_vec[m];
        int n6 = u_gq8d.element_p6_id_vec[m];
        int n7 = u_gq8d.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0, n1, n2, n3, n4, n5, n6, n7};

        // get points where the boundary is applied
        int loc = v_bq8d.element_value_loc_vec[k];  // 0 to 7

        // erase entire row
        a_mat.row(start_id_v + n_quad8_arr[loc]) *= 0.;
        b_vec.coeffRef(start_id_v + n_quad8_arr[loc]) = 0.;

    }

    // clear rows with value boundary elements (continuity)
    for (int k = 0; k < p_bq4d.num_element_value; k++)
    {

        // get id of element
        int m = p_bq4d.element_value_id_vec[k];

        // get id of points around quad4 element
        int n0 = p_gq4d.element_p0_id_vec[m];
        int n1 = p_gq4d.element_p1_id_vec[m];
        int n2 = p_gq4d.element_p2_id_vec[m];
        int n3 = p_gq4d.element_p3_id_vec[m];
        int n_quad4_arr[4] = {n0, n1, n2, n3};

        // get points where the boundary is applied
        int loc = p_bq4d.element_value_loc_vec[k];  // 0 to 3

        // erase entire row
        a_mat.row(start_id_p + n_quad4_arr[loc]) *= 0.;
        b_vec.coeffRef(start_id_p + n_quad4_arr[loc]) = 0.;

    }

    // iterate for each value boundary element (x-component of NSE)
    for (int k = 0; k < u_bq8d.num_element_value; k++)
    {

        // get id of element
        int m = u_bq8d.element_value_id_vec[k];

        // get id of points around quad8 element
        int n0 = u_gq8d.element_p0_id_vec[m];
        int n1 = u_gq8d.element_p1_id_vec[m];
        int n2 = u_gq8d.element_p2_id_vec[m];
        int n3 = u_gq8d.element_p3_id_vec[m];
        int n4 = u_gq8d.element_p4_id_vec[m];
        int n5 = u_gq8d.element_p5_id_vec[m];
        int n6 = u_gq8d.element_p6_id_vec[m];
        int n7 = u_gq8d.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0, n1, n2, n3, n4, n5, n6, n7};

        // get points where the boundary is applied
        int loc = u_bq8d.element_value_loc_vec[k];  // 0 to 7

        // identify boundary type
        int config_id = u_bq8d.element_value_config_id_vec[k];
        BoundaryConfigData bcd = u_bq8d.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcd.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            a_mat.coeffRef(start_id_u + n_quad8_arr[loc], start_id_u + n_quad8_arr[loc]) += 1.;
            b_vec.coeffRef(start_id_u + n_quad8_arr[loc]) += bcd.boundary_parameter_vec[0];

        }

    }

    // iterate for each value boundary element (y-component of NSE)
    for (int k = 0; k < v_bq8d.num_element_value; k++)
    {

        // get id of element
        int m = v_bq8d.element_value_id_vec[k];

        // get id of points around quad8 element
        int n0 = u_gq8d.element_p0_id_vec[m];
        int n1 = u_gq8d.element_p1_id_vec[m];
        int n2 = u_gq8d.element_p2_id_vec[m];
        int n3 = u_gq8d.element_p3_id_vec[m];
        int n4 = u_gq8d.element_p4_id_vec[m];
        int n5 = u_gq8d.element_p5_id_vec[m];
        int n6 = u_gq8d.element_p6_id_vec[m];
        int n7 = u_gq8d.element_p7_id_vec[m];
        int n_quad8_arr[8] = {n0, n1, n2, n3, n4, n5, n6, n7};

        // get points where the boundary is applied
        int loc = v_bq8d.element_value_loc_vec[k];  // 0 to 7

        // identify boundary type
        int config_id = v_bq8d.element_value_config_id_vec[k];
        BoundaryConfigData bcd = v_bq8d.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcd.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            a_mat.coeffRef(start_id_v + n_quad8_arr[loc], start_id_v + n_quad8_arr[loc]) += 1.;
            b_vec.coeffRef(start_id_v + n_quad8_arr[loc]) += bcd.boundary_parameter_vec[0];

        }

    }

    // iterate for each value boundary element (continuity)
    for (int k = 0; k < p_bq4d.num_element_value; k++)
    {

        // get id of element
        int m = p_bq4d.element_value_id_vec[k];

        // get id of points around quad4 element
        int n0 = p_gq4d.element_p0_id_vec[m];
        int n1 = p_gq4d.element_p1_id_vec[m];
        int n2 = p_gq4d.element_p2_id_vec[m];
        int n3 = p_gq4d.element_p3_id_vec[m];
        int n_quad4_arr[4] = {n0, n1, n2, n3};

        // get points where the boundary is applied
        int loc = p_bq4d.element_value_loc_vec[k];  // 0 to 3

        // identify boundary type
        int config_id = p_bq4d.element_value_config_id_vec[k];
        BoundaryConfigData bcd = p_bq4d.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcd.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            a_mat.coeffRef(start_id_p + n_quad4_arr[loc], start_id_p + n_quad4_arr[loc]) += 1.;
            b_vec.coeffRef(start_id_p + n_quad4_arr[loc]) += bcd.boundary_parameter_vec[0];

        }

    }

}



#endif
