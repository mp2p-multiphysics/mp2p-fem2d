#ifndef MODULE_HEAT_TRANSIENT_QUAD4
#define MODULE_HEAT_TRANSIENT_QUAD4
#include <vector>
#include "Eigen/Eigen"
#include "integral_quad4.hpp"
#include "grid_quad4.hpp"
#include "scalar_quad4.hpp"

struct HeatTransientQuad4IntegralStruct
{
    
    // integrals with one shape function
    std::vector<std::vector<double>> integral_Ni_quad4_vec;
    
    // integrals with two shape functions
    std::vector<std::vector<std::vector<double>>> integral_div_Ni_quad4_dot_div_Nj_quad4_vec;
    std::vector<std::vector<std::vector<double>>> integral_Ni_quad4_Nj_quad4_vec;

};

class HeatTransientQuad4Class
{

    public:

    // variables
    GridQuad4Struct gq4s;
    BoundaryQuad4Struct bq4s;

    // functions
    void integral_evaluate();
    void matrix_fill(
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec, Eigen::VectorXd &x_last_timestep_vec,
        ScalarQuad4Class &density_sq4c, ScalarQuad4Class &heat_capacity_sq4c, ScalarQuad4Class &thermal_conductivity_sq4c, ScalarQuad4Class &heat_generation_sq4c,
        double dt, int start_id
    );

    // struct with integrals
    HeatTransientQuad4IntegralStruct hsq4is;

    // constructor
    HeatTransientQuad4Class(GridQuad4Struct &gq4s_in, BoundaryQuad4Struct &bq4s_in)
    {
        gq4s = gq4s_in;
        bq4s = bq4s_in;
    }

};

void HeatTransientQuad4Class::integral_evaluate()
{
    
    // iterate for each element
    for (int m = 0; m < gq4s.num_element; m++)
    {

        // get id of points around element
        int n0 = gq4s.element_p0_id_vec[m];
        int n1 = gq4s.element_p1_id_vec[m];
        int n2 = gq4s.element_p2_id_vec[m];
        int n3 = gq4s.element_p3_id_vec[m];

        // get x coordinates of points around element
        double x0 = gq4s.point_pos_x_vec[n0];
        double x1 = gq4s.point_pos_x_vec[n1];
        double x2 = gq4s.point_pos_x_vec[n2];
        double x3 = gq4s.point_pos_x_vec[n3];
        double x_arr[4] = {x0, x1, x2, x3};

        // get y coordinates of points around element
        double y0 = gq4s.point_pos_y_vec[n0];
        double y1 = gq4s.point_pos_y_vec[n1];
        double y2 = gq4s.point_pos_y_vec[n2];
        double y3 = gq4s.point_pos_y_vec[n3];
        double y_arr[4] = {y0, y1, y2, y3};

        // temporarily stores values
        double part_val = 0;

        // evaluate integral_Ni_quad4_vec
        std::vector<double> integral_Ni_quad4_part_i_vec;
        for (int i = 0; i < 4; i++)
        {
            part_val = integral_Ni_quad4(i, x_arr, y_arr);
            integral_Ni_quad4_part_i_vec.push_back(part_val);
        }

        // evaluate integral_div_Ni_quad4_dot_div_Nj_quad4_vec
        std::vector<std::vector<double>> integral_div_Ni_quad4_dot_div_Nj_quad4_part_i_vec;
        for (int i = 0; i < 4; i++)
        {
            std::vector<double> integral_div_Ni_quad4_dot_div_Nj_quad4_part_ij_vec;
            for (int j = 0; j < 4; j++)
            {
                part_val = integral_div_Ni_quad4_dot_div_Nj_quad4(i, j, x_arr, y_arr);
                integral_div_Ni_quad4_dot_div_Nj_quad4_part_ij_vec.push_back(part_val);
            }
            integral_div_Ni_quad4_dot_div_Nj_quad4_part_i_vec.push_back(integral_div_Ni_quad4_dot_div_Nj_quad4_part_ij_vec);
        }

        // evaluate integral_Ni_quad4_Nj_quad4_vec
        std::vector<std::vector<double>> integral_Ni_quad4_Nj_quad4_part_i_vec;
        for (int i = 0; i < 4; i++)
        {
            std::vector<double> integral_Ni_quad4_Nj_quad4_part_ij_vec;
            for (int j = 0; j < 4; j++)
            {
                part_val = integral_Ni_quad4_Nj_quad4(i, j, x_arr, y_arr);
                integral_Ni_quad4_Nj_quad4_part_ij_vec.push_back(part_val);
            }
            integral_Ni_quad4_Nj_quad4_part_i_vec.push_back(integral_Ni_quad4_Nj_quad4_part_ij_vec);
        }

        // store integrals to struct
        hsq4is.integral_Ni_quad4_vec.push_back(integral_Ni_quad4_part_i_vec);
        hsq4is.integral_div_Ni_quad4_dot_div_Nj_quad4_vec.push_back(integral_div_Ni_quad4_dot_div_Nj_quad4_part_i_vec);
        hsq4is.integral_Ni_quad4_Nj_quad4_vec.push_back(integral_Ni_quad4_Nj_quad4_part_i_vec);

    }

}

void HeatTransientQuad4Class::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec, Eigen::VectorXd &x_last_timestep_vec,
    ScalarQuad4Class &density_sq4c, ScalarQuad4Class &heat_capacity_sq4c, ScalarQuad4Class &thermal_conductivity_sq4c, ScalarQuad4Class &heat_generation_sq4c,
    double dt, int start_id
)
{

    // initialize
    double inverse_dt = 1./dt;
    Eigen::SparseMatrix<double> help_mat(gq4s.num_point, gq4s.num_point);  // helper matrix used to make b_vec

    // iterate for each grid element
    for (int m = 0; m < gq4s.num_element; m++)
    {

        // get id of points around element
        int n0 = gq4s.element_p0_id_vec[m];
        int n1 = gq4s.element_p1_id_vec[m];
        int n2 = gq4s.element_p2_id_vec[m];
        int n3 = gq4s.element_p3_id_vec[m];
        int n_arr[4] = {n0, n1, n2, n3};

        // get density of points around element
        int density0 = density_sq4c.scalar_vec[n0];
        int density1 = density_sq4c.scalar_vec[n1];
        int density2 = density_sq4c.scalar_vec[n2];
        int density3 = density_sq4c.scalar_vec[n3];
        int density_arr[4] = {density0, density1, density2, density3};

        // get heat capacity of points around element
        int heatcap0 = heat_capacity_sq4c.scalar_vec[n0];
        int heatcap1 = heat_capacity_sq4c.scalar_vec[n1];
        int heatcap2 = heat_capacity_sq4c.scalar_vec[n2];
        int heatcap3 = heat_capacity_sq4c.scalar_vec[n3];
        int heatcap_arr[4] = {heatcap0, heatcap1, heatcap2, heatcap3};

        // get thermal conductivity of points around element
        int thermcond0 = thermal_conductivity_sq4c.scalar_vec[n0];
        int thermcond1 = thermal_conductivity_sq4c.scalar_vec[n1];
        int thermcond2 = thermal_conductivity_sq4c.scalar_vec[n2];
        int thermcond3 = thermal_conductivity_sq4c.scalar_vec[n3];
        int thermcond_arr[4] = {thermcond0, thermcond1, thermcond2, thermcond3};

        // get heat generation of points around element
        int heatgen0 = heat_generation_sq4c.scalar_vec[n0];
        int heatgen1 = heat_generation_sq4c.scalar_vec[n1];
        int heatgen2 = heat_generation_sq4c.scalar_vec[n2];
        int heatgen3 = heat_generation_sq4c.scalar_vec[n3];
        int heatgen_arr[4] = {heatgen0, heatgen1, heatgen2, heatgen3};

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable [u0 ... uN v0 ... vN p0 ... pN]

        // calculate a_mat coefficients
        for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            a_mat.coeffRef(start_id + n_arr[i], start_id + n_arr[j]) += inverse_dt*density_arr[i]*heatcap_arr[i]*hsq4is.integral_Ni_quad4_Nj_quad4_vec[m][i][j] + 0.5*thermcond_arr[i]*hsq4is.integral_div_Ni_quad4_dot_div_Nj_quad4_vec[m][i][j];
        }}

        // calculate help_mat coefficients
        for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            help_mat.coeffRef(n_arr[i], n_arr[j]) += inverse_dt*density_arr[i]*heatcap_arr[i]*hsq4is.integral_Ni_quad4_Nj_quad4_vec[m][i][j] - 0.5*thermcond_arr[i]*hsq4is.integral_div_Ni_quad4_dot_div_Nj_quad4_vec[m][i][j];
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 4; i++)
        {
            b_vec.coeffRef(start_id + n_arr[i]) += heatgen_arr[i]*hsq4is.integral_Ni_quad4_vec[m][i];
        }

    }

    // add additional terms for b_vec
    b_vec += help_mat*x_last_timestep_vec(Eigen::seqN(start_id, gq4s.num_point));

    // iterate for each flux boundary element
    for (int k = 0; k < bq4s.num_element_flux; k++)
    {

        // get id of element
        int m = bq4s.element_flux_id_vec[k];

        // get id of points around element
        int n0 = gq4s.element_p0_id_vec[m];
        int n1 = gq4s.element_p1_id_vec[m];
        int n2 = gq4s.element_p2_id_vec[m];
        int n3 = gq4s.element_p3_id_vec[m];
        int n_arr[4] = {n0, n1, n2, n3};

        // get x coordinates of points around element
        double x0 = gq4s.point_pos_x_vec[n0];
        double x1 = gq4s.point_pos_x_vec[n1];
        double x2 = gq4s.point_pos_x_vec[n2];
        double x3 = gq4s.point_pos_x_vec[n3];
        double x_arr[4] = {x0, x1, x2, x3};

        // get y coordinates of points around element
        double y0 = gq4s.point_pos_y_vec[n0];
        double y1 = gq4s.point_pos_y_vec[n1];
        double y2 = gq4s.point_pos_y_vec[n2];
        double y3 = gq4s.point_pos_y_vec[n3];
        double y_arr[4] = {y0, y1, y2, y3};

        // get points where the boundary is applied
        int a = bq4s.element_flux_pa_loc_vec[k];  // 0, 1, 2, or 3
        int b = bq4s.element_flux_pb_loc_vec[k];  // 0, 1, 2, or 3

        // identify boundary type
        int config_id = bq4s.element_flux_config_id_vec[k];
        BoundaryConfigQuad4Struct bcq4s = bq4s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcq4s.boundary_type_str == "neumann")
        {

            // calculate distance from point a to b
            double xa = x_arr[a];
            double xb = x_arr[b];
            double ya = y_arr[a];
            double yb = y_arr[b];
            double dist_pa_pb = sqrt((xa - xb)*(xa - xb) + (ya - yb)*(ya - yb));

            // add to b_vec
            b_vec.coeffRef(start_id + n_arr[a]) += -0.5 * bcq4s.boundary_parameter_vec[0] * dist_pa_pb;
            b_vec.coeffRef(start_id + n_arr[b]) += -0.5 * bcq4s.boundary_parameter_vec[0] * dist_pa_pb;

        }

    }

    // clear rows with value boundary elements
    for (int k = 0; k < bq4s.num_element_value; k++)
    {

        // get id of element
        int m = bq4s.element_value_id_vec[k];

        // get id of points around element
        int n0 = gq4s.element_p0_id_vec[m];
        int n1 = gq4s.element_p1_id_vec[m];
        int n2 = gq4s.element_p2_id_vec[m];
        int n3 = gq4s.element_p3_id_vec[m];
        int n_arr[4] = {n0, n1, n2, n3};

        // get points where the boundary is applied
        int a = bq4s.element_value_pa_loc_vec[k];  // 1, 2, 3, or 4
        int b = bq4s.element_value_pb_loc_vec[k];  // 1, 2, 3, or 4

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            a_mat.row(start_id + n_arr[a]) *= 0.;
            b_vec.coeffRef(start_id + n_arr[a]) = 0.;
        }
        if (b != -1)
        {
            a_mat.row(start_id + n_arr[b]) *= 0.;
            b_vec.coeffRef(start_id + n_arr[b]) = 0.;
        }

    }

    // iterate for each value boundary element
    for (int k = 0; k < bq4s.num_element_value; k++)
    {

        // get id of element
        int m = bq4s.element_value_id_vec[k];

        // get id of points around element
        int n0 = gq4s.element_p0_id_vec[m];
        int n1 = gq4s.element_p1_id_vec[m];
        int n2 = gq4s.element_p2_id_vec[m];
        int n3 = gq4s.element_p3_id_vec[m];
        int n_arr[4] = {n0, n1, n2, n3};

        // get points where the boundary is applied
        int a = bq4s.element_value_pa_loc_vec[k];  // 1, 2, 3, or 4
        int b = bq4s.element_value_pb_loc_vec[k];  // 1, 2, 3, or 4
        
        // identify boundary type
        int config_id = bq4s.element_value_config_id_vec[k];
        BoundaryConfigQuad4Struct bcq4s = bq4s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bcq4s.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            if (a != -1)
            {
                a_mat.coeffRef(start_id + n_arr[a], start_id + n_arr[a]) += 1.;
                b_vec.coeffRef(start_id + n_arr[a]) += bcq4s.boundary_parameter_vec[0];
            }
            if (b != -1)
            {
                a_mat.coeffRef(start_id + n_arr[b], start_id + n_arr[b]) += 1.;
                b_vec.coeffRef(start_id + n_arr[b]) += bcq4s.boundary_parameter_vec[0];
            }

        }

    }

}

#endif
