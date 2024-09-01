#ifndef INTEGRAL_QUAD4
#define INTEGRAL_QUAD4
#include <vector>
#include "Eigen/Eigen"

double integral_div_Ni_quad4_dot_div_Nj_quad4(int i, int j, double x_quad4_arr[4], double y_quad4_arr[4])
{

    // unpack x values
    double x0 = x_quad4_arr[0];
    double x1 = x_quad4_arr[1];
    double x2 = x_quad4_arr[2];
    double x3 = x_quad4_arr[3];

    // unpack y values
    double y0 = y_quad4_arr[0];
    double y1 = y_quad4_arr[1];
    double y2 = y_quad4_arr[2];
    double y3 = y_quad4_arr[3];

    // initialize points and weights for gaussian integration
    // weights are 1 for 4-point integration; ignore weights
    const double M_1_SQRT_3 = 1./sqrt(3);
    double a_arr[4] = {+M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3};
    double b_arr[4] = {+M_1_SQRT_3, -M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 4; m++)
    {

        // get a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];

        // get derivatives of test function Ni
        double derivative_Ni_a = 0.;
        double derivative_Ni_b = 0.;
        switch (i)
        {
            case 0: derivative_Ni_a = +0.25*(b - 1.); derivative_Ni_b = +0.25*(a - 1.); break;
            case 1: derivative_Ni_a = -0.25*(b + 1.); derivative_Ni_b = -0.25*(a - 1.); break;
            case 2: derivative_Ni_a = +0.25*(b + 1.); derivative_Ni_b = +0.25*(a + 1.); break;
            case 3: derivative_Ni_a = -0.25*(b - 1.); derivative_Ni_b = -0.25*(a + 1.); break;
        }

        // get derivatives of test function Nj
        double derivative_Nj_a = 0.;
        double derivative_Nj_b = 0.;
        switch (j)
        {
            case 0: derivative_Nj_a = +0.25*(b - 1.); derivative_Nj_b = +0.25*(a - 1.); break;
            case 1: derivative_Nj_a = -0.25*(b + 1.); derivative_Nj_b = -0.25*(a - 1.); break;
            case 2: derivative_Nj_a = +0.25*(b + 1.); derivative_Nj_b = +0.25*(a + 1.); break;
            case 3: derivative_Nj_a = -0.25*(b - 1.); derivative_Nj_b = -0.25*(a + 1.); break;
        }

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

        // get vectors with derivatives of test functions wrt a and b
        Eigen::RowVector2d derivative_Ni_ab_vec;
        Eigen::RowVector2d derivative_Nj_ab_vec;
        derivative_Ni_ab_vec << derivative_Ni_a, derivative_Ni_b;
        derivative_Nj_ab_vec << derivative_Nj_a, derivative_Nj_b;

        // evaluate part of integral value
        integral_value += jacobian_determinant * (derivative_Ni_ab_vec*jacobian_inverse_mat).dot(derivative_Nj_ab_vec*jacobian_inverse_mat);

    }

    return integral_value;

}

double integral_Ni_quad4(int i, double x_quad4_arr[4], double y_quad4_arr[4])
{

    // unpack x values
    double x0 = x_quad4_arr[0];
    double x1 = x_quad4_arr[1];
    double x2 = x_quad4_arr[2];
    double x3 = x_quad4_arr[3];

    // unpack y values
    double y0 = y_quad4_arr[0];
    double y1 = y_quad4_arr[1];
    double y2 = y_quad4_arr[2];
    double y3 = y_quad4_arr[3];

    // initialize points and weights for gaussian integration
    // weights are 1 for 4-point integration; ignore weights
    const double M_1_SQRT_3 = 1./sqrt(3);
    double a_arr[4] = {+M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3};
    double b_arr[4] = {+M_1_SQRT_3, -M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 4; m++)
    {

        // get a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];

        // get test function Ni
        double Ni = 0.;
        switch (i)
        {
            case 0: Ni = 0.25*(1. - a)*(1. - b); break;
            case 1: Ni = 0.25*(1. - a)*(1. + b); break;
            case 2: Ni = 0.25*(1. + a)*(1. + b); break;
            case 3: Ni = 0.25*(1. + a)*(1. - b); break;
        }

        // get derivatives of x and y with respect to a and b
        double derivative_x_a = 0.25*(b*x0 - b*x1 + b*x2 - b*x3 - x0 - x1 + x2 + x3);
        double derivative_x_b = 0.25*(a*x0 - a*x1 + a*x2 - a*x3 - x0 + x1 + x2 - x3);
        double derivative_y_a = 0.25*(b*y0 - b*y1 + b*y2 - b*y3 - y0 - y1 + y2 + y3);
        double derivative_y_b = 0.25*(a*y0 - a*y1 + a*y2 - a*y3 - y0 + y1 + y2 - y3);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        double jacobian_determinant = jacobian_mat.determinant();

        // evaluate part of integral value
        integral_value += jacobian_determinant * Ni;

    }

    return integral_value;

}

double integral_Ni_quad4_Nj_quad4(int i, int j, double x_quad4_arr[4], double y_quad4_arr[4])
{

    // unpack x values
    double x0 = x_quad4_arr[0];
    double x1 = x_quad4_arr[1];
    double x2 = x_quad4_arr[2];
    double x3 = x_quad4_arr[3];

    // unpack y values
    double y0 = y_quad4_arr[0];
    double y1 = y_quad4_arr[1];
    double y2 = y_quad4_arr[2];
    double y3 = y_quad4_arr[3];

    // initialize points and weights for gaussian integration
    // weights are 1 for 4-point integration; ignore weights
    const double M_1_SQRT_3 = 1./sqrt(3);
    double a_arr[4] = {+M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3};
    double b_arr[4] = {+M_1_SQRT_3, -M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3};

    // initialize integral value
    double integral_value = 0.;

    // calcualte for each point
    for (int m = 0; m < 4; m++)
    {

        // get a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];

        // get test function Ni
        double Ni = 0.;
        switch (i)
        {
            case 0: Ni = 0.25*(1. - a)*(1. - b); break;
            case 1: Ni = 0.25*(1. - a)*(1. + b); break;
            case 2: Ni = 0.25*(1. + a)*(1. + b); break;
            case 3: Ni = 0.25*(1. + a)*(1. - b); break;
        }

        // get test function Nj
        double Nj = 0.;
        switch (j)
        {
            case 0: Nj = 0.25*(1. - a)*(1. - b); break;
            case 1: Nj = 0.25*(1. - a)*(1. + b); break;
            case 2: Nj = 0.25*(1. + a)*(1. + b); break;
            case 3: Nj = 0.25*(1. + a)*(1. - b); break;
        }

        // get derivatives of x and y with respect to a and b
        double derivative_x_a = 0.25*(b*x0 - b*x1 + b*x2 - b*x3 - x0 - x1 + x2 + x3);
        double derivative_x_b = 0.25*(a*x0 - a*x1 + a*x2 - a*x3 - x0 + x1 + x2 - x3);
        double derivative_y_a = 0.25*(b*y0 - b*y1 + b*y2 - b*y3 - y0 - y1 + y2 + y3);
        double derivative_y_b = 0.25*(a*y0 - a*y1 + a*y2 - a*y3 - y0 + y1 + y2 - y3);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        double jacobian_determinant = jacobian_mat.determinant();

        // evaluate part of integral value
        integral_value += jacobian_determinant * Ni * Nj;

    }

    return integral_value;

}

#endif
