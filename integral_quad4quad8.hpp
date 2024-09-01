#ifndef INTEGRAL_QUAD4QUAD8
#define INTEGRAL_QUAD4QUAD8
#include <vector>
#include "Eigen/Eigen"

double integral_Ni_quad4_derivative_Mj_quad8_x(int i, int j, double x_quad8_arr[8], double y_quad8_arr[8])
{

    // unpack quad8 x values 
    double x0 = x_quad8_arr[0];
    double x1 = x_quad8_arr[1];
    double x2 = x_quad8_arr[2];
    double x3 = x_quad8_arr[3];
    double x4 = x_quad8_arr[4];
    double x5 = x_quad8_arr[5];
    double x6 = x_quad8_arr[6];
    double x7 = x_quad8_arr[7];

    // unpack quad8 y values
    double y0 = y_quad8_arr[0];
    double y1 = y_quad8_arr[1];
    double y2 = y_quad8_arr[2];
    double y3 = y_quad8_arr[3];
    double y4 = y_quad8_arr[4];
    double y5 = y_quad8_arr[5];
    double y6 = y_quad8_arr[6];
    double y7 = y_quad8_arr[7];

    // initialize points and weights for gaussian integration
    const double M_SQRT_3_5 = sqrt(0.6);
    double a_arr[9] = {+M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., 0.};
    double b_arr[9] = {+M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., +M_SQRT_3_5, -M_SQRT_3_5, 0.};
    double w_arr[9] = {25./81., 25./81., 25./81., 25./81., 40./81., 40./81., 40./81., 40./81., 64./81.};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 9; m++)
    {

        // get weight and a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];
        double w = w_arr[m];

        // get quad4 test function Ni
        double Ni = 0.;
        switch (i)
        {
            case 0: Ni = 0.25*(1. - a)*(1. - b); break;
            case 1: Ni = 0.25*(1. - a)*(1. + b); break;
            case 2: Ni = 0.25*(1. + a)*(1. + b); break;
            case 3: Ni = 0.25*(1. + a)*(1. - b); break;
        }

        // get derivatives of quad8 test function Mj
        double derivative_Mj_a = 0.;
        double derivative_Mj_b = 0.;
        switch (j)
        {
            case 0: derivative_Mj_a = -0.25*(2.*a + b)*(b - 1.); derivative_Mj_b = -0.25*(a - 1.)*(a + 2.*b); break;
            case 1: derivative_Mj_a = +0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = +b*(a - 1.);               break;
            case 2: derivative_Mj_a = +0.25*(2.*a - b)*(b + 1.); derivative_Mj_b = +0.25*(a - 1.)*(a - 2.*b); break;
            case 3: derivative_Mj_a = -a*(b + 1.);               derivative_Mj_b = -0.50*(a - 1.)*(a + 1.);   break;
            case 4: derivative_Mj_a = +0.25*(2.*a + b)*(b + 1.); derivative_Mj_b = +0.25*(a + 1.)*(a + 2.*b); break;
            case 5: derivative_Mj_a = -0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = -b*(a + 1.);               break;
            case 6: derivative_Mj_a = -0.25*(2.*a - b)*(b - 1.); derivative_Mj_b = -0.25*(a + 1.)*(a - 2.*b); break;
            case 7: derivative_Mj_a = +a*(b - 1.);               derivative_Mj_b = +0.50*(a - 1.)*(a + 1.);   break;
        }

        // get derivatives of x and y with respect to a and b
        double derivative_x_a = -0.25*(2.*a*b*x0 - 2.*a*b*x2 + 4.*a*b*x3 - 2.*a*b*x4 + 2.*a*b*x6 - 4.*a*b*x7 - 2.*a*x0 - 2.*a*x2 + 4.*a*x3 - 2.*a*x4 - 2.*a*x6 + 4.*a*x7 + b*b*x0 - 2.*b*b*x1 + b*b*x2 - b*b*x4 + 2.*b*b*x5 - b*b*x6 - b*x0 + b*x2 - b*x4 + b*x6 + 2.*x1 - 2.*x5);
        double derivative_x_b = -0.25*(a*a*x0 - a*a*x2 + 2.*a*a*x3 - a*a*x4 + a*a*x6 - 2.*a*a*x7 + 2.*a*b*x0 - 4.*a*b*x1 + 2.*a*b*x2 - 2.*a*b*x4 + 4.*a*b*x5 - 2.*a*b*x6 - a*x0 + a*x2 - a*x4 + a*x6 - 2.*b*x0 + 4.*b*x1 - 2.*b*x2 - 2.*b*x4 + 4.*b*x5 - 2.*b*x6 - 2.*x3 + 2.*x7);
        double derivative_y_a = -0.25*(2.*a*b*y0 - 2.*a*b*y2 + 4.*a*b*y3 - 2.*a*b*y4 + 2.*a*b*y6 - 4.*a*b*y7 - 2.*a*y0 - 2.*a*y2 + 4.*a*y3 - 2.*a*y4 - 2.*a*y6 + 4.*a*y7 + b*b*y0 - 2.*b*b*y1 + b*b*y2 - b*b*y4 + 2.*b*b*y5 - b*b*y6 - b*y0 + b*y2 - b*y4 + b*y6 + 2.*y1 - 2.*y5);
        double derivative_y_b = -0.25*(a*a*y0 - a*a*y2 + 2.*a*a*y3 - a*a*y4 + a*a*y6 - 2.*a*a*y7 + 2.*a*b*y0 - 4.*a*b*y1 + 2.*a*b*y2 - 2.*a*b*y4 + 4.*a*b*y5 - 2.*a*b*y6 - a*y0 + a*y2 - a*y4 + a*y6 - 2.*b*y0 + 4.*b*y1 - 2.*b*y2 - 2.*b*y4 + 4.*b*y5 - 2.*b*y6 - 2.*y3 + 2.*y7);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
        double jacobian_determinant = jacobian_mat.determinant();

        // get vectors with derivatives of test functions wrt a and b
        Eigen::RowVector2d derivative_Mj_ab_vec;
        derivative_Mj_ab_vec << derivative_Mj_a, derivative_Mj_b;

        // initialize vector with [1 0]
        Eigen::RowVector2d x_component_vec;
        x_component_vec << 1, 0;

        // evaluate part of integral value
        integral_value += w * jacobian_determinant * Ni * x_component_vec.dot(derivative_Mj_ab_vec * jacobian_inverse_mat);

    }

    return integral_value;

}

double integral_Ni_quad4_derivative_Mj_quad8_y(int i, int j, double x_quad8_arr[8], double y_quad8_arr[8])
{

    // unpack quad8 x values 
    double x0 = x_quad8_arr[0];
    double x1 = x_quad8_arr[1];
    double x2 = x_quad8_arr[2];
    double x3 = x_quad8_arr[3];
    double x4 = x_quad8_arr[4];
    double x5 = x_quad8_arr[5];
    double x6 = x_quad8_arr[6];
    double x7 = x_quad8_arr[7];

    // unpack quad8 y values
    double y0 = y_quad8_arr[0];
    double y1 = y_quad8_arr[1];
    double y2 = y_quad8_arr[2];
    double y3 = y_quad8_arr[3];
    double y4 = y_quad8_arr[4];
    double y5 = y_quad8_arr[5];
    double y6 = y_quad8_arr[6];
    double y7 = y_quad8_arr[7];

    // initialize points and weights for gaussian integration
    const double M_SQRT_3_5 = sqrt(0.6);
    double a_arr[9] = {+M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., 0.};
    double b_arr[9] = {+M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., +M_SQRT_3_5, -M_SQRT_3_5, 0.};
    double w_arr[9] = {25./81., 25./81., 25./81., 25./81., 40./81., 40./81., 40./81., 40./81., 64./81.};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 9; m++)
    {

        // get weight and a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];
        double w = w_arr[m];

        // get quad4 test function Ni
        double Ni = 0.;
        switch (i)
        {
            case 0: Ni = 0.25*(1. - a)*(1. - b); break;
            case 1: Ni = 0.25*(1. - a)*(1. + b); break;
            case 2: Ni = 0.25*(1. + a)*(1. + b); break;
            case 3: Ni = 0.25*(1. + a)*(1. - b); break;
        }

        // get derivatives of quad8 test function Mj
        double derivative_Mj_a = 0.;
        double derivative_Mj_b = 0.;
        switch (j)
        {
            case 0: derivative_Mj_a = -0.25*(2.*a + b)*(b - 1.); derivative_Mj_b = -0.25*(a - 1.)*(a + 2.*b); break;
            case 1: derivative_Mj_a = +0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = +b*(a - 1.);               break;
            case 2: derivative_Mj_a = +0.25*(2.*a - b)*(b + 1.); derivative_Mj_b = +0.25*(a - 1.)*(a - 2.*b); break;
            case 3: derivative_Mj_a = -a*(b + 1.);               derivative_Mj_b = -0.50*(a - 1.)*(a + 1.);   break;
            case 4: derivative_Mj_a = +0.25*(2.*a + b)*(b + 1.); derivative_Mj_b = +0.25*(a + 1.)*(a + 2.*b); break;
            case 5: derivative_Mj_a = -0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = -b*(a + 1.);               break;
            case 6: derivative_Mj_a = -0.25*(2.*a - b)*(b - 1.); derivative_Mj_b = -0.25*(a + 1.)*(a - 2.*b); break;
            case 7: derivative_Mj_a = +a*(b - 1.);               derivative_Mj_b = +0.50*(a - 1.)*(a + 1.);   break;
        }

        // get derivatives of x and y with respect to a and b
        double derivative_x_a = -0.25*(2.*a*b*x0 - 2.*a*b*x2 + 4.*a*b*x3 - 2.*a*b*x4 + 2.*a*b*x6 - 4.*a*b*x7 - 2.*a*x0 - 2.*a*x2 + 4.*a*x3 - 2.*a*x4 - 2.*a*x6 + 4.*a*x7 + b*b*x0 - 2.*b*b*x1 + b*b*x2 - b*b*x4 + 2.*b*b*x5 - b*b*x6 - b*x0 + b*x2 - b*x4 + b*x6 + 2.*x1 - 2.*x5);
        double derivative_x_b = -0.25*(a*a*x0 - a*a*x2 + 2.*a*a*x3 - a*a*x4 + a*a*x6 - 2.*a*a*x7 + 2.*a*b*x0 - 4.*a*b*x1 + 2.*a*b*x2 - 2.*a*b*x4 + 4.*a*b*x5 - 2.*a*b*x6 - a*x0 + a*x2 - a*x4 + a*x6 - 2.*b*x0 + 4.*b*x1 - 2.*b*x2 - 2.*b*x4 + 4.*b*x5 - 2.*b*x6 - 2.*x3 + 2.*x7);
        double derivative_y_a = -0.25*(2.*a*b*y0 - 2.*a*b*y2 + 4.*a*b*y3 - 2.*a*b*y4 + 2.*a*b*y6 - 4.*a*b*y7 - 2.*a*y0 - 2.*a*y2 + 4.*a*y3 - 2.*a*y4 - 2.*a*y6 + 4.*a*y7 + b*b*y0 - 2.*b*b*y1 + b*b*y2 - b*b*y4 + 2.*b*b*y5 - b*b*y6 - b*y0 + b*y2 - b*y4 + b*y6 + 2.*y1 - 2.*y5);
        double derivative_y_b = -0.25*(a*a*y0 - a*a*y2 + 2.*a*a*y3 - a*a*y4 + a*a*y6 - 2.*a*a*y7 + 2.*a*b*y0 - 4.*a*b*y1 + 2.*a*b*y2 - 2.*a*b*y4 + 4.*a*b*y5 - 2.*a*b*y6 - a*y0 + a*y2 - a*y4 + a*y6 - 2.*b*y0 + 4.*b*y1 - 2.*b*y2 - 2.*b*y4 + 4.*b*y5 - 2.*b*y6 - 2.*y3 + 2.*y7);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
        double jacobian_determinant = jacobian_mat.determinant();

        // get vectors with derivatives of test functions wrt a and b
        Eigen::RowVector2d derivative_Mj_ab_vec;
        derivative_Mj_ab_vec << derivative_Mj_a, derivative_Mj_b;

        // initialize vector with [0 1]
        Eigen::RowVector2d y_component_vec;
        y_component_vec << 0, 1;

        // evaluate part of integral value
        integral_value += w * jacobian_determinant * Ni * y_component_vec.dot(derivative_Mj_ab_vec * jacobian_inverse_mat);

    }

    return integral_value;

}

double integral_Mi_quad8_u_derivative_Mj_quad8_x(int i, int j, double u_quad8_arr[8], double x_quad8_arr[8], double y_quad8_arr[8])
{
    
    // unpack quad8 x values 
    double x0 = x_quad8_arr[0];
    double x1 = x_quad8_arr[1];
    double x2 = x_quad8_arr[2];
    double x3 = x_quad8_arr[3];
    double x4 = x_quad8_arr[4];
    double x5 = x_quad8_arr[5];
    double x6 = x_quad8_arr[6];
    double x7 = x_quad8_arr[7];

    // unpack quad8 y values
    double y0 = y_quad8_arr[0];
    double y1 = y_quad8_arr[1];
    double y2 = y_quad8_arr[2];
    double y3 = y_quad8_arr[3];
    double y4 = y_quad8_arr[4];
    double y5 = y_quad8_arr[5];
    double y6 = y_quad8_arr[6];
    double y7 = y_quad8_arr[7];

    // initialize points and weights for gaussian integration
    const double M_SQRT_3_5 = sqrt(0.6);
    double a_arr[9] = {+M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., 0.};
    double b_arr[9] = {+M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., +M_SQRT_3_5, -M_SQRT_3_5, 0.};
    double w_arr[9] = {25./81., 25./81., 25./81., 25./81., 40./81., 40./81., 40./81., 40./81., 64./81.};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 9; m++)
    {

        // get weight and a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];
        double w = w_arr[m];

        // initialize u
        double u = 0;

        // iterate from M0 to M7 and solve for u
        for (int k = 0; k < 8; k++)
        {

            // get quad8 test function Mk
            double Mk = 0.;
            switch (k)
            {
                case 0: Mk = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
                case 1: Mk = 0.50*(1. - a)*(1. - b*b);             break;
                case 2: Mk = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
                case 3: Mk = 0.50*(1. - a*a)*(1. + b);             break;
                case 4: Mk = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
                case 5: Mk = 0.50*(1. + a)*(1. - b*b);             break;
                case 6: Mk = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
                case 7: Mk = 0.50*(1. - a*a)*(1 - b);              break;
            }

            // add contributions to u and dv/dx
            u += u_quad8_arr[k]*Mk;

        }

        // get quad8 test function Mi
        double Mi = 0.;
        switch (i)
        {
            case 0: Mi = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
            case 1: Mi = 0.50*(1. - a)*(1. - b*b);             break;
            case 2: Mi = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
            case 3: Mi = 0.50*(1. - a*a)*(1. + b);             break;
            case 4: Mi = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
            case 5: Mi = 0.50*(1. + a)*(1. - b*b);             break;
            case 6: Mi = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
            case 7: Mi = 0.50*(1. - a*a)*(1 - b);              break;
        }

        // get derivatives of quad8 test function Mj
        double derivative_Mj_a = 0.;
        double derivative_Mj_b = 0.;
        switch (j)
        {
            case 0: derivative_Mj_a = -0.25*(2.*a + b)*(b - 1.); derivative_Mj_b = -0.25*(a - 1.)*(a + 2.*b); break;
            case 1: derivative_Mj_a = +0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = +b*(a - 1.);               break;
            case 2: derivative_Mj_a = +0.25*(2.*a - b)*(b + 1.); derivative_Mj_b = +0.25*(a - 1.)*(a - 2.*b); break;
            case 3: derivative_Mj_a = -a*(b + 1.);               derivative_Mj_b = -0.50*(a - 1.)*(a + 1.);   break;
            case 4: derivative_Mj_a = +0.25*(2.*a + b)*(b + 1.); derivative_Mj_b = +0.25*(a + 1.)*(a + 2.*b); break;
            case 5: derivative_Mj_a = -0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = -b*(a + 1.);               break;
            case 6: derivative_Mj_a = -0.25*(2.*a - b)*(b - 1.); derivative_Mj_b = -0.25*(a + 1.)*(a - 2.*b); break;
            case 7: derivative_Mj_a = +a*(b - 1.);               derivative_Mj_b = +0.50*(a - 1.)*(a + 1.);   break;
        }

        // get derivatives of x and y with respect to a and b
        double derivative_x_a = -0.25*(2.*a*b*x0 - 2.*a*b*x2 + 4.*a*b*x3 - 2.*a*b*x4 + 2.*a*b*x6 - 4.*a*b*x7 - 2.*a*x0 - 2.*a*x2 + 4.*a*x3 - 2.*a*x4 - 2.*a*x6 + 4.*a*x7 + b*b*x0 - 2.*b*b*x1 + b*b*x2 - b*b*x4 + 2.*b*b*x5 - b*b*x6 - b*x0 + b*x2 - b*x4 + b*x6 + 2.*x1 - 2.*x5);
        double derivative_x_b = -0.25*(a*a*x0 - a*a*x2 + 2.*a*a*x3 - a*a*x4 + a*a*x6 - 2.*a*a*x7 + 2.*a*b*x0 - 4.*a*b*x1 + 2.*a*b*x2 - 2.*a*b*x4 + 4.*a*b*x5 - 2.*a*b*x6 - a*x0 + a*x2 - a*x4 + a*x6 - 2.*b*x0 + 4.*b*x1 - 2.*b*x2 - 2.*b*x4 + 4.*b*x5 - 2.*b*x6 - 2.*x3 + 2.*x7);
        double derivative_y_a = -0.25*(2.*a*b*y0 - 2.*a*b*y2 + 4.*a*b*y3 - 2.*a*b*y4 + 2.*a*b*y6 - 4.*a*b*y7 - 2.*a*y0 - 2.*a*y2 + 4.*a*y3 - 2.*a*y4 - 2.*a*y6 + 4.*a*y7 + b*b*y0 - 2.*b*b*y1 + b*b*y2 - b*b*y4 + 2.*b*b*y5 - b*b*y6 - b*y0 + b*y2 - b*y4 + b*y6 + 2.*y1 - 2.*y5);
        double derivative_y_b = -0.25*(a*a*y0 - a*a*y2 + 2.*a*a*y3 - a*a*y4 + a*a*y6 - 2.*a*a*y7 + 2.*a*b*y0 - 4.*a*b*y1 + 2.*a*b*y2 - 2.*a*b*y4 + 4.*a*b*y5 - 2.*a*b*y6 - a*y0 + a*y2 - a*y4 + a*y6 - 2.*b*y0 + 4.*b*y1 - 2.*b*y2 - 2.*b*y4 + 4.*b*y5 - 2.*b*y6 - 2.*y3 + 2.*y7);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
        double jacobian_determinant = jacobian_mat.determinant();

        // get vectors with derivatives of test functions wrt a and b
        Eigen::RowVector2d derivative_Mj_ab_vec;
        derivative_Mj_ab_vec << derivative_Mj_a, derivative_Mj_b;

        // initialize vector with [1 0]
        Eigen::RowVector2d x_component_vec;
        x_component_vec << 1, 0;

        // evaluate part of integral value
        integral_value += w * jacobian_determinant * Mi * u * x_component_vec.dot(derivative_Mj_ab_vec * jacobian_inverse_mat);

    }

    return integral_value;

}

double integral_Mi_quad8_u_derivative_Mj_quad8_y(int i, int j, double u_quad8_arr[8], double x_quad8_arr[8], double y_quad8_arr[8])
{
    
    // unpack quad8 x values 
    double x0 = x_quad8_arr[0];
    double x1 = x_quad8_arr[1];
    double x2 = x_quad8_arr[2];
    double x3 = x_quad8_arr[3];
    double x4 = x_quad8_arr[4];
    double x5 = x_quad8_arr[5];
    double x6 = x_quad8_arr[6];
    double x7 = x_quad8_arr[7];

    // unpack quad8 y values
    double y0 = y_quad8_arr[0];
    double y1 = y_quad8_arr[1];
    double y2 = y_quad8_arr[2];
    double y3 = y_quad8_arr[3];
    double y4 = y_quad8_arr[4];
    double y5 = y_quad8_arr[5];
    double y6 = y_quad8_arr[6];
    double y7 = y_quad8_arr[7];

    // initialize points and weights for gaussian integration
    const double M_SQRT_3_5 = sqrt(0.6);
    double a_arr[9] = {+M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., 0.};
    double b_arr[9] = {+M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., +M_SQRT_3_5, -M_SQRT_3_5, 0.};
    double w_arr[9] = {25./81., 25./81., 25./81., 25./81., 40./81., 40./81., 40./81., 40./81., 64./81.};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 9; m++)
    {

        // get weight and a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];
        double w = w_arr[m];

        // initialize u
        double u = 0;

        // iterate from M0 to M7 and solve for u
        for (int k = 0; k < 8; k++)
        {

            // get quad8 test function Mk
            double Mk = 0.;
            switch (k)
            {
                case 0: Mk = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
                case 1: Mk = 0.50*(1. - a)*(1. - b*b);             break;
                case 2: Mk = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
                case 3: Mk = 0.50*(1. - a*a)*(1. + b);             break;
                case 4: Mk = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
                case 5: Mk = 0.50*(1. + a)*(1. - b*b);             break;
                case 6: Mk = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
                case 7: Mk = 0.50*(1. - a*a)*(1 - b);              break;
            }

            // add contributions to u and dv/dx
            u += u_quad8_arr[k]*Mk;

        }

        // get quad8 test function Mi
        double Mi = 0.;
        switch (i)
        {
            case 0: Mi = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
            case 1: Mi = 0.50*(1. - a)*(1. - b*b);             break;
            case 2: Mi = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
            case 3: Mi = 0.50*(1. - a*a)*(1. + b);             break;
            case 4: Mi = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
            case 5: Mi = 0.50*(1. + a)*(1. - b*b);             break;
            case 6: Mi = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
            case 7: Mi = 0.50*(1. - a*a)*(1 - b);              break;
        }

        // get derivatives of quad8 test function Mj
        double derivative_Mj_a = 0.;
        double derivative_Mj_b = 0.;
        switch (j)
        {
            case 0: derivative_Mj_a = -0.25*(2.*a + b)*(b - 1.); derivative_Mj_b = -0.25*(a - 1.)*(a + 2.*b); break;
            case 1: derivative_Mj_a = +0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = +b*(a - 1.);               break;
            case 2: derivative_Mj_a = +0.25*(2.*a - b)*(b + 1.); derivative_Mj_b = +0.25*(a - 1.)*(a - 2.*b); break;
            case 3: derivative_Mj_a = -a*(b + 1.);               derivative_Mj_b = -0.50*(a - 1.)*(a + 1.);   break;
            case 4: derivative_Mj_a = +0.25*(2.*a + b)*(b + 1.); derivative_Mj_b = +0.25*(a + 1.)*(a + 2.*b); break;
            case 5: derivative_Mj_a = -0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = -b*(a + 1.);               break;
            case 6: derivative_Mj_a = -0.25*(2.*a - b)*(b - 1.); derivative_Mj_b = -0.25*(a + 1.)*(a - 2.*b); break;
            case 7: derivative_Mj_a = +a*(b - 1.);               derivative_Mj_b = +0.50*(a - 1.)*(a + 1.);   break;
        }

        // get derivatives of x and y with respect to a and b
        double derivative_x_a = -0.25*(2.*a*b*x0 - 2.*a*b*x2 + 4.*a*b*x3 - 2.*a*b*x4 + 2.*a*b*x6 - 4.*a*b*x7 - 2.*a*x0 - 2.*a*x2 + 4.*a*x3 - 2.*a*x4 - 2.*a*x6 + 4.*a*x7 + b*b*x0 - 2.*b*b*x1 + b*b*x2 - b*b*x4 + 2.*b*b*x5 - b*b*x6 - b*x0 + b*x2 - b*x4 + b*x6 + 2.*x1 - 2.*x5);
        double derivative_x_b = -0.25*(a*a*x0 - a*a*x2 + 2.*a*a*x3 - a*a*x4 + a*a*x6 - 2.*a*a*x7 + 2.*a*b*x0 - 4.*a*b*x1 + 2.*a*b*x2 - 2.*a*b*x4 + 4.*a*b*x5 - 2.*a*b*x6 - a*x0 + a*x2 - a*x4 + a*x6 - 2.*b*x0 + 4.*b*x1 - 2.*b*x2 - 2.*b*x4 + 4.*b*x5 - 2.*b*x6 - 2.*x3 + 2.*x7);
        double derivative_y_a = -0.25*(2.*a*b*y0 - 2.*a*b*y2 + 4.*a*b*y3 - 2.*a*b*y4 + 2.*a*b*y6 - 4.*a*b*y7 - 2.*a*y0 - 2.*a*y2 + 4.*a*y3 - 2.*a*y4 - 2.*a*y6 + 4.*a*y7 + b*b*y0 - 2.*b*b*y1 + b*b*y2 - b*b*y4 + 2.*b*b*y5 - b*b*y6 - b*y0 + b*y2 - b*y4 + b*y6 + 2.*y1 - 2.*y5);
        double derivative_y_b = -0.25*(a*a*y0 - a*a*y2 + 2.*a*a*y3 - a*a*y4 + a*a*y6 - 2.*a*a*y7 + 2.*a*b*y0 - 4.*a*b*y1 + 2.*a*b*y2 - 2.*a*b*y4 + 4.*a*b*y5 - 2.*a*b*y6 - a*y0 + a*y2 - a*y4 + a*y6 - 2.*b*y0 + 4.*b*y1 - 2.*b*y2 - 2.*b*y4 + 4.*b*y5 - 2.*b*y6 - 2.*y3 + 2.*y7);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
        double jacobian_determinant = jacobian_mat.determinant();

        // get vectors with derivatives of test functions wrt a and b
        Eigen::RowVector2d derivative_Mj_ab_vec;
        derivative_Mj_ab_vec << derivative_Mj_a, derivative_Mj_b;

        // initialize vector with [0 1]
        Eigen::RowVector2d y_component_vec;
        y_component_vec << 0, 1;

        // evaluate part of integral value
        integral_value += w * jacobian_determinant * Mi * u * y_component_vec.dot(derivative_Mj_ab_vec * jacobian_inverse_mat);

    }

    return integral_value;

}

double integral_Mi_quad8_derivative_Nj_quad4_x(int i, int j, double x_quad8_arr[8], double y_quad8_arr[8])
{

    // unpack quad8 x values 
    double x0 = x_quad8_arr[0];
    double x1 = x_quad8_arr[1];
    double x2 = x_quad8_arr[2];
    double x3 = x_quad8_arr[3];
    double x4 = x_quad8_arr[4];
    double x5 = x_quad8_arr[5];
    double x6 = x_quad8_arr[6];
    double x7 = x_quad8_arr[7];

    // unpack quad8 y values
    double y0 = y_quad8_arr[0];
    double y1 = y_quad8_arr[1];
    double y2 = y_quad8_arr[2];
    double y3 = y_quad8_arr[3];
    double y4 = y_quad8_arr[4];
    double y5 = y_quad8_arr[5];
    double y6 = y_quad8_arr[6];
    double y7 = y_quad8_arr[7];

    // initialize points and weights for gaussian integration
    const double M_SQRT_3_5 = sqrt(0.6);
    double a_arr[9] = {+M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., 0.};
    double b_arr[9] = {+M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., +M_SQRT_3_5, -M_SQRT_3_5, 0.};
    double w_arr[9] = {25./81., 25./81., 25./81., 25./81., 40./81., 40./81., 40./81., 40./81., 64./81.};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 9; m++)
    {

        // get weight and a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];
        double w = w_arr[m];

        // get quad8 test function Mi
        double Mi = 0.;
        switch (i)
        {
            case 0: Mi = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
            case 1: Mi = 0.50*(1. - a)*(1. - b*b);             break;
            case 2: Mi = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
            case 3: Mi = 0.50*(1. - a*a)*(1. + b);             break;
            case 4: Mi = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
            case 5: Mi = 0.50*(1. + a)*(1. - b*b);             break;
            case 6: Mi = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
            case 7: Mi = 0.50*(1. - a*a)*(1 - b);              break;
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
        double derivative_x_a = -0.25*(2.*a*b*x0 - 2.*a*b*x2 + 4.*a*b*x3 - 2.*a*b*x4 + 2.*a*b*x6 - 4.*a*b*x7 - 2.*a*x0 - 2.*a*x2 + 4.*a*x3 - 2.*a*x4 - 2.*a*x6 + 4.*a*x7 + b*b*x0 - 2.*b*b*x1 + b*b*x2 - b*b*x4 + 2.*b*b*x5 - b*b*x6 - b*x0 + b*x2 - b*x4 + b*x6 + 2.*x1 - 2.*x5);
        double derivative_x_b = -0.25*(a*a*x0 - a*a*x2 + 2.*a*a*x3 - a*a*x4 + a*a*x6 - 2.*a*a*x7 + 2.*a*b*x0 - 4.*a*b*x1 + 2.*a*b*x2 - 2.*a*b*x4 + 4.*a*b*x5 - 2.*a*b*x6 - a*x0 + a*x2 - a*x4 + a*x6 - 2.*b*x0 + 4.*b*x1 - 2.*b*x2 - 2.*b*x4 + 4.*b*x5 - 2.*b*x6 - 2.*x3 + 2.*x7);
        double derivative_y_a = -0.25*(2.*a*b*y0 - 2.*a*b*y2 + 4.*a*b*y3 - 2.*a*b*y4 + 2.*a*b*y6 - 4.*a*b*y7 - 2.*a*y0 - 2.*a*y2 + 4.*a*y3 - 2.*a*y4 - 2.*a*y6 + 4.*a*y7 + b*b*y0 - 2.*b*b*y1 + b*b*y2 - b*b*y4 + 2.*b*b*y5 - b*b*y6 - b*y0 + b*y2 - b*y4 + b*y6 + 2.*y1 - 2.*y5);
        double derivative_y_b = -0.25*(a*a*y0 - a*a*y2 + 2.*a*a*y3 - a*a*y4 + a*a*y6 - 2.*a*a*y7 + 2.*a*b*y0 - 4.*a*b*y1 + 2.*a*b*y2 - 2.*a*b*y4 + 4.*a*b*y5 - 2.*a*b*y6 - a*y0 + a*y2 - a*y4 + a*y6 - 2.*b*y0 + 4.*b*y1 - 2.*b*y2 - 2.*b*y4 + 4.*b*y5 - 2.*b*y6 - 2.*y3 + 2.*y7);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
        double jacobian_determinant = jacobian_mat.determinant();

        // get vectors with derivatives of test functions wrt a and b
        Eigen::RowVector2d derivative_Nj_ab_vec;
        derivative_Nj_ab_vec << derivative_Nj_a, derivative_Nj_b;

        // initialize vector with [1 0]
        Eigen::RowVector2d x_component_vec;
        x_component_vec << 1, 0;

        // evaluate part of integral value
        integral_value += w * jacobian_determinant * Mi * x_component_vec.dot(derivative_Nj_ab_vec * jacobian_inverse_mat);

    }

    return integral_value;

}

double integral_Mi_quad8_derivative_Nj_quad4_y(int i, int j, double x_quad8_arr[8], double y_quad8_arr[8])
{

    // unpack quad8 x values 
    double x0 = x_quad8_arr[0];
    double x1 = x_quad8_arr[1];
    double x2 = x_quad8_arr[2];
    double x3 = x_quad8_arr[3];
    double x4 = x_quad8_arr[4];
    double x5 = x_quad8_arr[5];
    double x6 = x_quad8_arr[6];
    double x7 = x_quad8_arr[7];

    // unpack quad8 y values
    double y0 = y_quad8_arr[0];
    double y1 = y_quad8_arr[1];
    double y2 = y_quad8_arr[2];
    double y3 = y_quad8_arr[3];
    double y4 = y_quad8_arr[4];
    double y5 = y_quad8_arr[5];
    double y6 = y_quad8_arr[6];
    double y7 = y_quad8_arr[7];

    // initialize points and weights for gaussian integration
    const double M_SQRT_3_5 = sqrt(0.6);
    double a_arr[9] = {+M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., 0.};
    double b_arr[9] = {+M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., +M_SQRT_3_5, -M_SQRT_3_5, 0.};
    double w_arr[9] = {25./81., 25./81., 25./81., 25./81., 40./81., 40./81., 40./81., 40./81., 64./81.};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 9; m++)
    {

        // get weight and a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];
        double w = w_arr[m];

        // get quad8 test function Mi
        double Mi = 0.;
        switch (i)
        {
            case 0: Mi = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
            case 1: Mi = 0.50*(1. - a)*(1. - b*b);             break;
            case 2: Mi = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
            case 3: Mi = 0.50*(1. - a*a)*(1. + b);             break;
            case 4: Mi = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
            case 5: Mi = 0.50*(1. + a)*(1. - b*b);             break;
            case 6: Mi = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
            case 7: Mi = 0.50*(1. - a*a)*(1 - b);              break;
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
        double derivative_x_a = -0.25*(2.*a*b*x0 - 2.*a*b*x2 + 4.*a*b*x3 - 2.*a*b*x4 + 2.*a*b*x6 - 4.*a*b*x7 - 2.*a*x0 - 2.*a*x2 + 4.*a*x3 - 2.*a*x4 - 2.*a*x6 + 4.*a*x7 + b*b*x0 - 2.*b*b*x1 + b*b*x2 - b*b*x4 + 2.*b*b*x5 - b*b*x6 - b*x0 + b*x2 - b*x4 + b*x6 + 2.*x1 - 2.*x5);
        double derivative_x_b = -0.25*(a*a*x0 - a*a*x2 + 2.*a*a*x3 - a*a*x4 + a*a*x6 - 2.*a*a*x7 + 2.*a*b*x0 - 4.*a*b*x1 + 2.*a*b*x2 - 2.*a*b*x4 + 4.*a*b*x5 - 2.*a*b*x6 - a*x0 + a*x2 - a*x4 + a*x6 - 2.*b*x0 + 4.*b*x1 - 2.*b*x2 - 2.*b*x4 + 4.*b*x5 - 2.*b*x6 - 2.*x3 + 2.*x7);
        double derivative_y_a = -0.25*(2.*a*b*y0 - 2.*a*b*y2 + 4.*a*b*y3 - 2.*a*b*y4 + 2.*a*b*y6 - 4.*a*b*y7 - 2.*a*y0 - 2.*a*y2 + 4.*a*y3 - 2.*a*y4 - 2.*a*y6 + 4.*a*y7 + b*b*y0 - 2.*b*b*y1 + b*b*y2 - b*b*y4 + 2.*b*b*y5 - b*b*y6 - b*y0 + b*y2 - b*y4 + b*y6 + 2.*y1 - 2.*y5);
        double derivative_y_b = -0.25*(a*a*y0 - a*a*y2 + 2.*a*a*y3 - a*a*y4 + a*a*y6 - 2.*a*a*y7 + 2.*a*b*y0 - 4.*a*b*y1 + 2.*a*b*y2 - 2.*a*b*y4 + 4.*a*b*y5 - 2.*a*b*y6 - a*y0 + a*y2 - a*y4 + a*y6 - 2.*b*y0 + 4.*b*y1 - 2.*b*y2 - 2.*b*y4 + 4.*b*y5 - 2.*b*y6 - 2.*y3 + 2.*y7);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
        double jacobian_determinant = jacobian_mat.determinant();

        // get vectors with derivatives of test functions wrt a and b
        Eigen::RowVector2d derivative_Nj_ab_vec;
        derivative_Nj_ab_vec << derivative_Nj_a, derivative_Nj_b;

        // initialize vector with [0 1]
        Eigen::RowVector2d y_component_vec;
        y_component_vec << 0, 1;

        // evaluate part of integral value
        integral_value += w * jacobian_determinant * Mi * y_component_vec.dot(derivative_Nj_ab_vec * jacobian_inverse_mat);

    }

    return integral_value;

}

double integral_div_Mi_quad8_dot_div_Mj_quad8(int i, int j, double x_quad8_arr[8], double y_quad8_arr[8])
{

    // unpack quad8 x values 
    double x0 = x_quad8_arr[0];
    double x1 = x_quad8_arr[1];
    double x2 = x_quad8_arr[2];
    double x3 = x_quad8_arr[3];
    double x4 = x_quad8_arr[4];
    double x5 = x_quad8_arr[5];
    double x6 = x_quad8_arr[6];
    double x7 = x_quad8_arr[7];

    // unpack quad8 y values
    double y0 = y_quad8_arr[0];
    double y1 = y_quad8_arr[1];
    double y2 = y_quad8_arr[2];
    double y3 = y_quad8_arr[3];
    double y4 = y_quad8_arr[4];
    double y5 = y_quad8_arr[5];
    double y6 = y_quad8_arr[6];
    double y7 = y_quad8_arr[7];

    // initialize points and weights for gaussian integration
    const double M_SQRT_3_5 = sqrt(0.6);
    double a_arr[9] = {+M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., 0.};
    double b_arr[9] = {+M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., +M_SQRT_3_5, -M_SQRT_3_5, 0.};
    double w_arr[9] = {25./81., 25./81., 25./81., 25./81., 40./81., 40./81., 40./81., 40./81., 64./81.};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 9; m++)
    {

        // get weight and a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];
        double w = w_arr[m];

        // get derivatives of quad8 test function Mi
        double derivative_Mi_a = 0.;
        double derivative_Mi_b = 0.;
        switch (i)
        {
            case 0: derivative_Mi_a = -0.25*(2.*a + b)*(b - 1.); derivative_Mi_b = -0.25*(a - 1.)*(a + 2.*b); break;
            case 1: derivative_Mi_a = +0.50*(b - 1.)*(b + 1.);   derivative_Mi_b = +b*(a - 1.);               break;
            case 2: derivative_Mi_a = +0.25*(2.*a - b)*(b + 1.); derivative_Mi_b = +0.25*(a - 1.)*(a - 2.*b); break;
            case 3: derivative_Mi_a = -a*(b + 1.);               derivative_Mi_b = -0.50*(a - 1.)*(a + 1.);   break;
            case 4: derivative_Mi_a = +0.25*(2.*a + b)*(b + 1.); derivative_Mi_b = +0.25*(a + 1.)*(a + 2.*b); break;
            case 5: derivative_Mi_a = -0.50*(b - 1.)*(b + 1.);   derivative_Mi_b = -b*(a + 1.);               break;
            case 6: derivative_Mi_a = -0.25*(2.*a - b)*(b - 1.); derivative_Mi_b = -0.25*(a + 1.)*(a - 2.*b); break;
            case 7: derivative_Mi_a = +a*(b - 1.);               derivative_Mi_b = +0.50*(a - 1.)*(a + 1.);   break;
        }

        // get derivatives of quad8 test function Mj
        double derivative_Mj_a = 0.;
        double derivative_Mj_b = 0.;
        switch (j)
        {
            case 0: derivative_Mj_a = -0.25*(2.*a + b)*(b - 1.); derivative_Mj_b = -0.25*(a - 1.)*(a + 2.*b); break;
            case 1: derivative_Mj_a = +0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = +b*(a - 1.);               break;
            case 2: derivative_Mj_a = +0.25*(2.*a - b)*(b + 1.); derivative_Mj_b = +0.25*(a - 1.)*(a - 2.*b); break;
            case 3: derivative_Mj_a = -a*(b + 1.);               derivative_Mj_b = -0.50*(a - 1.)*(a + 1.);   break;
            case 4: derivative_Mj_a = +0.25*(2.*a + b)*(b + 1.); derivative_Mj_b = +0.25*(a + 1.)*(a + 2.*b); break;
            case 5: derivative_Mj_a = -0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = -b*(a + 1.);               break;
            case 6: derivative_Mj_a = -0.25*(2.*a - b)*(b - 1.); derivative_Mj_b = -0.25*(a + 1.)*(a - 2.*b); break;
            case 7: derivative_Mj_a = +a*(b - 1.);               derivative_Mj_b = +0.50*(a - 1.)*(a + 1.);   break;
        }

        // get derivatives of x and y with respect to a and b
        double derivative_x_a = -0.25*(2.*a*b*x0 - 2.*a*b*x2 + 4.*a*b*x3 - 2.*a*b*x4 + 2.*a*b*x6 - 4.*a*b*x7 - 2.*a*x0 - 2.*a*x2 + 4.*a*x3 - 2.*a*x4 - 2.*a*x6 + 4.*a*x7 + b*b*x0 - 2.*b*b*x1 + b*b*x2 - b*b*x4 + 2.*b*b*x5 - b*b*x6 - b*x0 + b*x2 - b*x4 + b*x6 + 2.*x1 - 2.*x5);
        double derivative_x_b = -0.25*(a*a*x0 - a*a*x2 + 2.*a*a*x3 - a*a*x4 + a*a*x6 - 2.*a*a*x7 + 2.*a*b*x0 - 4.*a*b*x1 + 2.*a*b*x2 - 2.*a*b*x4 + 4.*a*b*x5 - 2.*a*b*x6 - a*x0 + a*x2 - a*x4 + a*x6 - 2.*b*x0 + 4.*b*x1 - 2.*b*x2 - 2.*b*x4 + 4.*b*x5 - 2.*b*x6 - 2.*x3 + 2.*x7);
        double derivative_y_a = -0.25*(2.*a*b*y0 - 2.*a*b*y2 + 4.*a*b*y3 - 2.*a*b*y4 + 2.*a*b*y6 - 4.*a*b*y7 - 2.*a*y0 - 2.*a*y2 + 4.*a*y3 - 2.*a*y4 - 2.*a*y6 + 4.*a*y7 + b*b*y0 - 2.*b*b*y1 + b*b*y2 - b*b*y4 + 2.*b*b*y5 - b*b*y6 - b*y0 + b*y2 - b*y4 + b*y6 + 2.*y1 - 2.*y5);
        double derivative_y_b = -0.25*(a*a*y0 - a*a*y2 + 2.*a*a*y3 - a*a*y4 + a*a*y6 - 2.*a*a*y7 + 2.*a*b*y0 - 4.*a*b*y1 + 2.*a*b*y2 - 2.*a*b*y4 + 4.*a*b*y5 - 2.*a*b*y6 - a*y0 + a*y2 - a*y4 + a*y6 - 2.*b*y0 + 4.*b*y1 - 2.*b*y2 - 2.*b*y4 + 4.*b*y5 - 2.*b*y6 - 2.*y3 + 2.*y7);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
        double jacobian_determinant = jacobian_mat.determinant();

        // get vectors with derivatives of test functions wrt a and b
        Eigen::RowVector2d derivative_Mi_ab_vec;
        Eigen::RowVector2d derivative_Mj_ab_vec;
        derivative_Mi_ab_vec << derivative_Mi_a, derivative_Mi_b;
        derivative_Mj_ab_vec << derivative_Mj_a, derivative_Mj_b;

        // evaluate part of integral value
        integral_value += w * jacobian_determinant * (derivative_Mi_ab_vec * jacobian_inverse_mat).dot(derivative_Mj_ab_vec * jacobian_inverse_mat);

    }

    return integral_value;

}

double integral_Mi_quad8(int i, double x_quad8_arr[8], double y_quad8_arr[8])
{

    // unpack quad8 x values 
    double x0 = x_quad8_arr[0];
    double x1 = x_quad8_arr[1];
    double x2 = x_quad8_arr[2];
    double x3 = x_quad8_arr[3];
    double x4 = x_quad8_arr[4];
    double x5 = x_quad8_arr[5];
    double x6 = x_quad8_arr[6];
    double x7 = x_quad8_arr[7];

    // unpack quad8 y values
    double y0 = y_quad8_arr[0];
    double y1 = y_quad8_arr[1];
    double y2 = y_quad8_arr[2];
    double y3 = y_quad8_arr[3];
    double y4 = y_quad8_arr[4];
    double y5 = y_quad8_arr[5];
    double y6 = y_quad8_arr[6];
    double y7 = y_quad8_arr[7];

    // initialize points and weights for gaussian integration
    const double M_SQRT_3_5 = sqrt(0.6);
    double a_arr[9] = {+M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., 0.};
    double b_arr[9] = {+M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., +M_SQRT_3_5, -M_SQRT_3_5, 0.};
    double w_arr[9] = {25./81., 25./81., 25./81., 25./81., 40./81., 40./81., 40./81., 40./81., 64./81.};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 9; m++)
    {

        // get weight and a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];
        double w = w_arr[m];

        // get quad8 test function Mi
        double Mi = 0.;
        switch (i)
        {
            case 0: Mi = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
            case 1: Mi = 0.50*(1. - a)*(1. - b*b);             break;
            case 2: Mi = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
            case 3: Mi = 0.50*(1. - a*a)*(1. + b);             break;
            case 4: Mi = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
            case 5: Mi = 0.50*(1. + a)*(1. - b*b);             break;
            case 6: Mi = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
            case 7: Mi = 0.50*(1. - a*a)*(1 - b);              break;
        }        

        // get derivatives of x and y with respect to a and b
        double derivative_x_a = -0.25*(2.*a*b*x0 - 2.*a*b*x2 + 4.*a*b*x3 - 2.*a*b*x4 + 2.*a*b*x6 - 4.*a*b*x7 - 2.*a*x0 - 2.*a*x2 + 4.*a*x3 - 2.*a*x4 - 2.*a*x6 + 4.*a*x7 + b*b*x0 - 2.*b*b*x1 + b*b*x2 - b*b*x4 + 2.*b*b*x5 - b*b*x6 - b*x0 + b*x2 - b*x4 + b*x6 + 2.*x1 - 2.*x5);
        double derivative_x_b = -0.25*(a*a*x0 - a*a*x2 + 2.*a*a*x3 - a*a*x4 + a*a*x6 - 2.*a*a*x7 + 2.*a*b*x0 - 4.*a*b*x1 + 2.*a*b*x2 - 2.*a*b*x4 + 4.*a*b*x5 - 2.*a*b*x6 - a*x0 + a*x2 - a*x4 + a*x6 - 2.*b*x0 + 4.*b*x1 - 2.*b*x2 - 2.*b*x4 + 4.*b*x5 - 2.*b*x6 - 2.*x3 + 2.*x7);
        double derivative_y_a = -0.25*(2.*a*b*y0 - 2.*a*b*y2 + 4.*a*b*y3 - 2.*a*b*y4 + 2.*a*b*y6 - 4.*a*b*y7 - 2.*a*y0 - 2.*a*y2 + 4.*a*y3 - 2.*a*y4 - 2.*a*y6 + 4.*a*y7 + b*b*y0 - 2.*b*b*y1 + b*b*y2 - b*b*y4 + 2.*b*b*y5 - b*b*y6 - b*y0 + b*y2 - b*y4 + b*y6 + 2.*y1 - 2.*y5);
        double derivative_y_b = -0.25*(a*a*y0 - a*a*y2 + 2.*a*a*y3 - a*a*y4 + a*a*y6 - 2.*a*a*y7 + 2.*a*b*y0 - 4.*a*b*y1 + 2.*a*b*y2 - 2.*a*b*y4 + 4.*a*b*y5 - 2.*a*b*y6 - a*y0 + a*y2 - a*y4 + a*y6 - 2.*b*y0 + 4.*b*y1 - 2.*b*y2 - 2.*b*y4 + 4.*b*y5 - 2.*b*y6 - 2.*y3 + 2.*y7);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        double jacobian_determinant = jacobian_mat.determinant();

        // evaluate part of integral value
        integral_value += w * jacobian_determinant * Mi;

    }

    return integral_value;

}

double integral_Mi_quad8_Mj_quad8(int i, int j, double x_quad8_arr[8], double y_quad8_arr[8])
{

    // unpack quad8 x values 
    double x0 = x_quad8_arr[0];
    double x1 = x_quad8_arr[1];
    double x2 = x_quad8_arr[2];
    double x3 = x_quad8_arr[3];
    double x4 = x_quad8_arr[4];
    double x5 = x_quad8_arr[5];
    double x6 = x_quad8_arr[6];
    double x7 = x_quad8_arr[7];

    // unpack quad8 y values
    double y0 = y_quad8_arr[0];
    double y1 = y_quad8_arr[1];
    double y2 = y_quad8_arr[2];
    double y3 = y_quad8_arr[3];
    double y4 = y_quad8_arr[4];
    double y5 = y_quad8_arr[5];
    double y6 = y_quad8_arr[6];
    double y7 = y_quad8_arr[7];

    // initialize points and weights for gaussian integration
    const double M_SQRT_3_5 = sqrt(0.6);
    double a_arr[9] = {+M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., 0.};
    double b_arr[9] = {+M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., +M_SQRT_3_5, -M_SQRT_3_5, 0.};
    double w_arr[9] = {25./81., 25./81., 25./81., 25./81., 40./81., 40./81., 40./81., 40./81., 64./81.};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 9; m++)
    {

        // get weight and a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];
        double w = w_arr[m];

        // get quad8 test function Mi
        double Mi = 0.;
        switch (i)
        {
            case 0: Mi = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
            case 1: Mi = 0.50*(1. - a)*(1. - b*b);             break;
            case 2: Mi = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
            case 3: Mi = 0.50*(1. - a*a)*(1. + b);             break;
            case 4: Mi = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
            case 5: Mi = 0.50*(1. + a)*(1. - b*b);             break;
            case 6: Mi = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
            case 7: Mi = 0.50*(1. - a*a)*(1 - b);              break;
        }        

        // get quad8 test function Mj
        double Mj = 0.;
        switch (j)
        {
            case 0: Mj = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
            case 1: Mj = 0.50*(1. - a)*(1. - b*b);             break;
            case 2: Mj = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
            case 3: Mj = 0.50*(1. - a*a)*(1. + b);             break;
            case 4: Mj = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
            case 5: Mj = 0.50*(1. + a)*(1. - b*b);             break;
            case 6: Mj = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
            case 7: Mj = 0.50*(1. - a*a)*(1 - b);              break;
        }      

        // get derivatives of x and y with respect to a and b
        double derivative_x_a = -0.25*(2.*a*b*x0 - 2.*a*b*x2 + 4.*a*b*x3 - 2.*a*b*x4 + 2.*a*b*x6 - 4.*a*b*x7 - 2.*a*x0 - 2.*a*x2 + 4.*a*x3 - 2.*a*x4 - 2.*a*x6 + 4.*a*x7 + b*b*x0 - 2.*b*b*x1 + b*b*x2 - b*b*x4 + 2.*b*b*x5 - b*b*x6 - b*x0 + b*x2 - b*x4 + b*x6 + 2.*x1 - 2.*x5);
        double derivative_x_b = -0.25*(a*a*x0 - a*a*x2 + 2.*a*a*x3 - a*a*x4 + a*a*x6 - 2.*a*a*x7 + 2.*a*b*x0 - 4.*a*b*x1 + 2.*a*b*x2 - 2.*a*b*x4 + 4.*a*b*x5 - 2.*a*b*x6 - a*x0 + a*x2 - a*x4 + a*x6 - 2.*b*x0 + 4.*b*x1 - 2.*b*x2 - 2.*b*x4 + 4.*b*x5 - 2.*b*x6 - 2.*x3 + 2.*x7);
        double derivative_y_a = -0.25*(2.*a*b*y0 - 2.*a*b*y2 + 4.*a*b*y3 - 2.*a*b*y4 + 2.*a*b*y6 - 4.*a*b*y7 - 2.*a*y0 - 2.*a*y2 + 4.*a*y3 - 2.*a*y4 - 2.*a*y6 + 4.*a*y7 + b*b*y0 - 2.*b*b*y1 + b*b*y2 - b*b*y4 + 2.*b*b*y5 - b*b*y6 - b*y0 + b*y2 - b*y4 + b*y6 + 2.*y1 - 2.*y5);
        double derivative_y_b = -0.25*(a*a*y0 - a*a*y2 + 2.*a*a*y3 - a*a*y4 + a*a*y6 - 2.*a*a*y7 + 2.*a*b*y0 - 4.*a*b*y1 + 2.*a*b*y2 - 2.*a*b*y4 + 4.*a*b*y5 - 2.*a*b*y6 - a*y0 + a*y2 - a*y4 + a*y6 - 2.*b*y0 + 4.*b*y1 - 2.*b*y2 - 2.*b*y4 + 4.*b*y5 - 2.*b*y6 - 2.*y3 + 2.*y7);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        double jacobian_determinant = jacobian_mat.determinant();

        // evaluate part of integral value
        integral_value += w * jacobian_determinant * Mi * Mj;

    }

    return integral_value;

}

double integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_x(int i, int j, int k, double x_quad8_arr[8], double y_quad8_arr[8])
{

    // unpack quad8 x values 
    double x0 = x_quad8_arr[0];
    double x1 = x_quad8_arr[1];
    double x2 = x_quad8_arr[2];
    double x3 = x_quad8_arr[3];
    double x4 = x_quad8_arr[4];
    double x5 = x_quad8_arr[5];
    double x6 = x_quad8_arr[6];
    double x7 = x_quad8_arr[7];

    // unpack quad8 y values
    double y0 = y_quad8_arr[0];
    double y1 = y_quad8_arr[1];
    double y2 = y_quad8_arr[2];
    double y3 = y_quad8_arr[3];
    double y4 = y_quad8_arr[4];
    double y5 = y_quad8_arr[5];
    double y6 = y_quad8_arr[6];
    double y7 = y_quad8_arr[7];

    // initialize points and weights for gaussian integration
    const double M_SQRT_3_5 = sqrt(0.6);
    double a_arr[9] = {+M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., 0.};
    double b_arr[9] = {+M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., +M_SQRT_3_5, -M_SQRT_3_5, 0.};
    double w_arr[9] = {25./81., 25./81., 25./81., 25./81., 40./81., 40./81., 40./81., 40./81., 64./81.};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 9; m++)
    {

        // get weight and a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];
        double w = w_arr[m];

        // get quad8 test function Mi
        double Mi = 0.;
        switch (i)
        {
            case 0: Mi = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
            case 1: Mi = 0.50*(1. - a)*(1. - b*b);             break;
            case 2: Mi = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
            case 3: Mi = 0.50*(1. - a*a)*(1. + b);             break;
            case 4: Mi = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
            case 5: Mi = 0.50*(1. + a)*(1. - b*b);             break;
            case 6: Mi = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
            case 7: Mi = 0.50*(1. - a*a)*(1 - b);              break;
        }

        // get quad8 test function Mk
        double Mk = 0.;
        switch (k)
        {
            case 0: Mk = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
            case 1: Mk = 0.50*(1. - a)*(1. - b*b);             break;
            case 2: Mk = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
            case 3: Mk = 0.50*(1. - a*a)*(1. + b);             break;
            case 4: Mk = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
            case 5: Mk = 0.50*(1. + a)*(1. - b*b);             break;
            case 6: Mk = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
            case 7: Mk = 0.50*(1. - a*a)*(1 - b);              break;
        }

        // get derivatives of quad8 test function Mj
        double derivative_Mj_a = 0.;
        double derivative_Mj_b = 0.;
        switch (j)
        {
            case 0: derivative_Mj_a = -0.25*(2.*a + b)*(b - 1.); derivative_Mj_b = -0.25*(a - 1.)*(a + 2.*b); break;
            case 1: derivative_Mj_a = +0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = +b*(a - 1.);               break;
            case 2: derivative_Mj_a = +0.25*(2.*a - b)*(b + 1.); derivative_Mj_b = +0.25*(a - 1.)*(a - 2.*b); break;
            case 3: derivative_Mj_a = -a*(b + 1.);               derivative_Mj_b = -0.50*(a - 1.)*(a + 1.);   break;
            case 4: derivative_Mj_a = +0.25*(2.*a + b)*(b + 1.); derivative_Mj_b = +0.25*(a + 1.)*(a + 2.*b); break;
            case 5: derivative_Mj_a = -0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = -b*(a + 1.);               break;
            case 6: derivative_Mj_a = -0.25*(2.*a - b)*(b - 1.); derivative_Mj_b = -0.25*(a + 1.)*(a - 2.*b); break;
            case 7: derivative_Mj_a = +a*(b - 1.);               derivative_Mj_b = +0.50*(a - 1.)*(a + 1.);   break;
        }

        // get derivatives of x and y with respect to a and b
        double derivative_x_a = -0.25*(2.*a*b*x0 - 2.*a*b*x2 + 4.*a*b*x3 - 2.*a*b*x4 + 2.*a*b*x6 - 4.*a*b*x7 - 2.*a*x0 - 2.*a*x2 + 4.*a*x3 - 2.*a*x4 - 2.*a*x6 + 4.*a*x7 + b*b*x0 - 2.*b*b*x1 + b*b*x2 - b*b*x4 + 2.*b*b*x5 - b*b*x6 - b*x0 + b*x2 - b*x4 + b*x6 + 2.*x1 - 2.*x5);
        double derivative_x_b = -0.25*(a*a*x0 - a*a*x2 + 2.*a*a*x3 - a*a*x4 + a*a*x6 - 2.*a*a*x7 + 2.*a*b*x0 - 4.*a*b*x1 + 2.*a*b*x2 - 2.*a*b*x4 + 4.*a*b*x5 - 2.*a*b*x6 - a*x0 + a*x2 - a*x4 + a*x6 - 2.*b*x0 + 4.*b*x1 - 2.*b*x2 - 2.*b*x4 + 4.*b*x5 - 2.*b*x6 - 2.*x3 + 2.*x7);
        double derivative_y_a = -0.25*(2.*a*b*y0 - 2.*a*b*y2 + 4.*a*b*y3 - 2.*a*b*y4 + 2.*a*b*y6 - 4.*a*b*y7 - 2.*a*y0 - 2.*a*y2 + 4.*a*y3 - 2.*a*y4 - 2.*a*y6 + 4.*a*y7 + b*b*y0 - 2.*b*b*y1 + b*b*y2 - b*b*y4 + 2.*b*b*y5 - b*b*y6 - b*y0 + b*y2 - b*y4 + b*y6 + 2.*y1 - 2.*y5);
        double derivative_y_b = -0.25*(a*a*y0 - a*a*y2 + 2.*a*a*y3 - a*a*y4 + a*a*y6 - 2.*a*a*y7 + 2.*a*b*y0 - 4.*a*b*y1 + 2.*a*b*y2 - 2.*a*b*y4 + 4.*a*b*y5 - 2.*a*b*y6 - a*y0 + a*y2 - a*y4 + a*y6 - 2.*b*y0 + 4.*b*y1 - 2.*b*y2 - 2.*b*y4 + 4.*b*y5 - 2.*b*y6 - 2.*y3 + 2.*y7);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
        double jacobian_determinant = jacobian_mat.determinant();

        // get vectors with derivatives of test functions wrt a and b
        Eigen::RowVector2d derivative_Mj_ab_vec;
        derivative_Mj_ab_vec << derivative_Mj_a, derivative_Mj_b;

        // initialize vector with [1 0]
        Eigen::RowVector2d x_component_vec;
        x_component_vec << 1, 0;

        // evaluate part of integral value
        integral_value += w * jacobian_determinant * Mi * Mk * x_component_vec.dot(derivative_Mj_ab_vec * jacobian_inverse_mat);

    }

    return integral_value;

}

double integral_Mi_quad8_Mk_quad8_derivative_Mj_quad8_y(int i, int j, int k, double x_quad8_arr[8], double y_quad8_arr[8])
{

    // unpack quad8 x values 
    double x0 = x_quad8_arr[0];
    double x1 = x_quad8_arr[1];
    double x2 = x_quad8_arr[2];
    double x3 = x_quad8_arr[3];
    double x4 = x_quad8_arr[4];
    double x5 = x_quad8_arr[5];
    double x6 = x_quad8_arr[6];
    double x7 = x_quad8_arr[7];

    // unpack quad8 y values
    double y0 = y_quad8_arr[0];
    double y1 = y_quad8_arr[1];
    double y2 = y_quad8_arr[2];
    double y3 = y_quad8_arr[3];
    double y4 = y_quad8_arr[4];
    double y5 = y_quad8_arr[5];
    double y6 = y_quad8_arr[6];
    double y7 = y_quad8_arr[7];

    // initialize points and weights for gaussian integration
    const double M_SQRT_3_5 = sqrt(0.6);
    double a_arr[9] = {+M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., 0.};
    double b_arr[9] = {+M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, 0., 0., +M_SQRT_3_5, -M_SQRT_3_5, 0.};
    double w_arr[9] = {25./81., 25./81., 25./81., 25./81., 40./81., 40./81., 40./81., 40./81., 64./81.};

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 9; m++)
    {

        // get weight and a and b values where function is evaluated
        double a = a_arr[m];
        double b = b_arr[m];
        double w = w_arr[m];

        // get quad8 test function Mi
        double Mi = 0.;
        switch (i)
        {
            case 0: Mi = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
            case 1: Mi = 0.50*(1. - a)*(1. - b*b);             break;
            case 2: Mi = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
            case 3: Mi = 0.50*(1. - a*a)*(1. + b);             break;
            case 4: Mi = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
            case 5: Mi = 0.50*(1. + a)*(1. - b*b);             break;
            case 6: Mi = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
            case 7: Mi = 0.50*(1. - a*a)*(1 - b);              break;
        }

        // get quad8 test function Mk
        double Mk = 0.;
        switch (k)
        {
            case 0: Mk = 0.25*(1. - a)*(1. - b)*(-a - b - 1.); break;
            case 1: Mk = 0.50*(1. - a)*(1. - b*b);             break;
            case 2: Mk = 0.25*(1. - a)*(1. + b)*(-a + b - 1.); break;
            case 3: Mk = 0.50*(1. - a*a)*(1. + b);             break;
            case 4: Mk = 0.25*(1. + a)*(1. + b)*(+a + b - 1.); break;
            case 5: Mk = 0.50*(1. + a)*(1. - b*b);             break;
            case 6: Mk = 0.25*(1. + a)*(1. - b)*(+a - b - 1.); break;
            case 7: Mk = 0.50*(1. - a*a)*(1 - b);              break;
        }

        // get derivatives of quad8 test function Mj
        double derivative_Mj_a = 0.;
        double derivative_Mj_b = 0.;
        switch (j)
        {
            case 0: derivative_Mj_a = -0.25*(2.*a + b)*(b - 1.); derivative_Mj_b = -0.25*(a - 1.)*(a + 2.*b); break;
            case 1: derivative_Mj_a = +0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = +b*(a - 1.);               break;
            case 2: derivative_Mj_a = +0.25*(2.*a - b)*(b + 1.); derivative_Mj_b = +0.25*(a - 1.)*(a - 2.*b); break;
            case 3: derivative_Mj_a = -a*(b + 1.);               derivative_Mj_b = -0.50*(a - 1.)*(a + 1.);   break;
            case 4: derivative_Mj_a = +0.25*(2.*a + b)*(b + 1.); derivative_Mj_b = +0.25*(a + 1.)*(a + 2.*b); break;
            case 5: derivative_Mj_a = -0.50*(b - 1.)*(b + 1.);   derivative_Mj_b = -b*(a + 1.);               break;
            case 6: derivative_Mj_a = -0.25*(2.*a - b)*(b - 1.); derivative_Mj_b = -0.25*(a + 1.)*(a - 2.*b); break;
            case 7: derivative_Mj_a = +a*(b - 1.);               derivative_Mj_b = +0.50*(a - 1.)*(a + 1.);   break;
        }

        // get derivatives of x and y with respect to a and b
        double derivative_x_a = -0.25*(2.*a*b*x0 - 2.*a*b*x2 + 4.*a*b*x3 - 2.*a*b*x4 + 2.*a*b*x6 - 4.*a*b*x7 - 2.*a*x0 - 2.*a*x2 + 4.*a*x3 - 2.*a*x4 - 2.*a*x6 + 4.*a*x7 + b*b*x0 - 2.*b*b*x1 + b*b*x2 - b*b*x4 + 2.*b*b*x5 - b*b*x6 - b*x0 + b*x2 - b*x4 + b*x6 + 2.*x1 - 2.*x5);
        double derivative_x_b = -0.25*(a*a*x0 - a*a*x2 + 2.*a*a*x3 - a*a*x4 + a*a*x6 - 2.*a*a*x7 + 2.*a*b*x0 - 4.*a*b*x1 + 2.*a*b*x2 - 2.*a*b*x4 + 4.*a*b*x5 - 2.*a*b*x6 - a*x0 + a*x2 - a*x4 + a*x6 - 2.*b*x0 + 4.*b*x1 - 2.*b*x2 - 2.*b*x4 + 4.*b*x5 - 2.*b*x6 - 2.*x3 + 2.*x7);
        double derivative_y_a = -0.25*(2.*a*b*y0 - 2.*a*b*y2 + 4.*a*b*y3 - 2.*a*b*y4 + 2.*a*b*y6 - 4.*a*b*y7 - 2.*a*y0 - 2.*a*y2 + 4.*a*y3 - 2.*a*y4 - 2.*a*y6 + 4.*a*y7 + b*b*y0 - 2.*b*b*y1 + b*b*y2 - b*b*y4 + 2.*b*b*y5 - b*b*y6 - b*y0 + b*y2 - b*y4 + b*y6 + 2.*y1 - 2.*y5);
        double derivative_y_b = -0.25*(a*a*y0 - a*a*y2 + 2.*a*a*y3 - a*a*y4 + a*a*y6 - 2.*a*a*y7 + 2.*a*b*y0 - 4.*a*b*y1 + 2.*a*b*y2 - 2.*a*b*y4 + 4.*a*b*y5 - 2.*a*b*y6 - a*y0 + a*y2 - a*y4 + a*y6 - 2.*b*y0 + 4.*b*y1 - 2.*b*y2 - 2.*b*y4 + 4.*b*y5 - 2.*b*y6 - 2.*y3 + 2.*y7);

        // get jacobian and its inverse and determinant
        Eigen::Matrix2d jacobian_mat;
        jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
        Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
        double jacobian_determinant = jacobian_mat.determinant();

        // get vectors with derivatives of test functions wrt a and b
        Eigen::RowVector2d derivative_Mj_ab_vec;
        derivative_Mj_ab_vec << derivative_Mj_a, derivative_Mj_b;

        // initialize vector with [0 1]
        Eigen::RowVector2d y_component_vec;
        y_component_vec << 0, 1;

        // evaluate part of integral value
        integral_value += w * jacobian_determinant * Mi * Mk * y_component_vec.dot(derivative_Mj_ab_vec * jacobian_inverse_mat);

    }

    return integral_value;

}

#endif
