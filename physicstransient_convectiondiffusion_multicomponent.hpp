#ifndef PHYSICSTRANSIENT_CONVECTIONDIFFUSION_MULTICOMPONENT
#define PHYSICSTRANSIENT_CONVECTIONDIFFUSION_MULTICOMPONENT
#include <unordered_map>
#include <vector>
#include "Eigen/Eigen"
#include "boundary_physicsgroup_vector.hpp"
#include "boundary_physicsgroup.hpp"
#include "container_typedef.hpp"
#include "integral_physicsgroup_vector.hpp"
#include "integral_physicsgroup.hpp"
#include "mesh_physicsgroup.hpp"
#include "physicstransient_base.hpp"
#include "scalar_fieldgroup_matrix.hpp"
#include "scalar_fieldgroup_vector.hpp"
#include "scalar_fieldgroup.hpp"
#include "variable_fieldgroup_vector.hpp"
#include "variable_fieldgroup.hpp"

class PhysicsTransientConvectionDiffusionMulticomponent : public PhysicsTransientBase
{
    /*

    Multi-component transient convection-diffusion equation.   

    a * du/dt = -div(-b * grad(u) + u * v) + c wherein:

    a1 * d(u1)/dt = -div(-b11 * grad(u1) - b12 * grad(u2) - b13 * grad(u3) - ... + u1 * v) + c1
    a2 * d(u2)/dt = -div(-b21 * grad(u1) - b22 * grad(u2) - b23 * grad(u3) - ... + u2 * v) + c2
    a3 * d(u3)/dt = -div(-b31 * grad(u1) - b32 * grad(u2) - b33 * grad(u3) - ... + u3 * v) + c3
    ...

    Variables
    =========
    mesh_physics_in : MeshPhysicsGroup
        Meshes where this physics is applied to.
    boundary_physics_ptr_vec_in : BoundaryPhysicsGroupVector
        Boundary conditions pertinent to this physics.
    integral_physics_ptr_vec_in : IntegralPhysicsGroupVector
        Test function integrals of the meshes.
    value_field_ptr_vec_in : VariableFieldGroupVector
        u in a * du/dt = -div(-b * grad(u) + u * v) + c.
        This will be solved for by the matrix equation.
    derivativecoefficient_field_ptr_vec_in : ScalarFieldGroupVector
        a in a * du/dt = -div(-b * grad(u) + u * v) + c.
    diffusioncoefficient_field_ptr_mat_in : ScalarFieldGroupMatrix
        b in a * du/dt = -div(-b * grad(u) + u * v) + c.
    velocity_x_field_in : ScalarFieldGroup
        x-component of v in a * du/dt = -div(-b * grad(u) + u * v) + c.
    velocity_x_field_in : ScalarFieldGroup
        y-component of v in a * du/dt = -div(-b * grad(u) + u * v) + c.
    generationcoefficient_field_ptr_vec_in : ScalarFieldGroupVector
        c in a * du/dt = -div(-b * grad(u) + u * v) + c.

    Functions
    =========
    matrix_fill : void
        Fill up the matrix equation Ax = b with entries as dictated by the physics. 
    set_start_row : void
        Sets the starting row in A and b where entries are filled up.
    get_start_row : int
        Returns the starting row.
    get_variable_field_ptr_vec() : vector<VariableFieldGroup*>
        Returns the vector containing pointers to VariableFieldGroup objects tied to this physics.

    */

    public:

    // physics groups
    MeshPhysicsGroup *mesh_physics_ptr;
    BoundaryPhysicsGroupVector boundary_physics_ptr_vec;
    IntegralPhysicsGroupVector integral_physics_ptr_vec;

    // field groups
    VariableFieldGroupVector value_field_ptr_vec;
    ScalarFieldGroupVector derivativecoefficient_field_ptr_vec;
    ScalarFieldGroupMatrix diffusioncoefficient_field_ptr_mat;
    ScalarFieldGroup *velocity_x_field_ptr;
    ScalarFieldGroup *velocity_y_field_ptr;
    ScalarFieldGroupVector generationcoefficient_field_ptr_vec;

    // vector of variable fields
    std::vector<VariableFieldGroup*> variable_field_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(
        Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
        Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt
    );
    void set_start_row(int start_row_in);
    int get_start_row();
    std::vector<VariableFieldGroup*> get_variable_field_ptr_vec();

    // default constructor
    PhysicsTransientConvectionDiffusionMulticomponent()
    {

    }

    // constructor
    PhysicsTransientConvectionDiffusionMulticomponent
    (
        MeshPhysicsGroup &mesh_physics_in, BoundaryPhysicsGroupVector &boundary_physics_ptr_vec_in, IntegralPhysicsGroupVector &integral_physics_ptr_vec_in,
        VariableFieldGroupVector &value_field_ptr_vec_in,
        ScalarFieldGroupVector &derivativecoefficient_field_ptr_vec_in, ScalarFieldGroupMatrix &diffusioncoefficient_field_ptr_mat_in, ScalarFieldGroup &velocity_x_field_in, ScalarFieldGroup &velocity_y_field_in, ScalarFieldGroupVector &generationcoefficient_field_ptr_vec_in
    )
    {
        
        // store physics groups
        mesh_physics_ptr = &mesh_physics_in;
        boundary_physics_ptr_vec = boundary_physics_ptr_vec_in;
        integral_physics_ptr_vec = integral_physics_ptr_vec_in;

        // store field 
        value_field_ptr_vec = value_field_ptr_vec_in;
        derivativecoefficient_field_ptr_vec = derivativecoefficient_field_ptr_vec_in;
        diffusioncoefficient_field_ptr_mat = diffusioncoefficient_field_ptr_mat_in;
        velocity_x_field_ptr = &velocity_x_field_in;
        velocity_y_field_ptr = &velocity_y_field_in;
        generationcoefficient_field_ptr_vec = generationcoefficient_field_ptr_vec_in;

        // vector of variable fields 
        variable_field_ptr_vec = value_field_ptr_vec.get_vector();

        // iterate through integrals and evaluate
        for (auto &integral_physics_ptr : integral_physics_ptr_vec.integral_physicsgroup_ptr_vec)
        {

            // calculate integrals
            integral_physics_ptr->evaluate_Ni_derivative();
            integral_physics_ptr->evaluate_integral_div_Ni_dot_div_Nj();
            integral_physics_ptr->evaluate_integral_Ni_derivative_Nj_x();
            integral_physics_ptr->evaluate_integral_Ni_derivative_Nj_y();
            integral_physics_ptr->evaluate_integral_Ni_Nj();
            integral_physics_ptr->evaluate_integral_Ni();
            integral_physics_ptr->evaluate_integral_Ni_Nj_derivative_Nk_x();
            integral_physics_ptr->evaluate_integral_Ni_Nj_derivative_Nk_y();

            // calculate boundary integrals
            integral_physics_ptr->evaluate_boundary_Ni_derivative();
            integral_physics_ptr->evaluate_boundary_normal();
            integral_physics_ptr->evaluate_boundary_integral_Ni();
            integral_physics_ptr->evaluate_boundary_integral_Ni_Nj();

        }

    }

    private:
    void matrix_fill_domain_t3
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
        Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
        MeshTri3 *mesh_ptr, BoundaryTri3 *boundary_ptr, IntegralTri3 *integral_ptr,
        VariableFieldGroup *value_field_ptr,
        ScalarTri3 *derivativecoefficient_ptr, ScalarTri3 *diffusioncoefficient_ptr, ScalarTri3 *velocity_x_ptr, ScalarTri3 *velocity_y_ptr, ScalarTri3 *generationcoefficient_ptr
    );
    void matrix_fill_domain_q4
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
        Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
        MeshQuad4 *mesh_ptr, BoundaryQuad4 *boundary_ptr, IntegralQuad4 *integral_ptr,
        VariableFieldGroup *value_field_ptr,
        ScalarQuad4 *derivativecoefficient_ptr, ScalarQuad4 *diffusioncoefficient_ptr, ScalarQuad4 *velocity_x_ptr, ScalarQuad4 *velocity_y_ptr, ScalarQuad4 *generationcoefficient_ptr
    );

};

void PhysicsTransientConvectionDiffusionMulticomponent::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
    Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt
)
{
    /*

    Fill up the matrix equation Ax(t+1) = Cx(t) + d with entries as dictated by the physics. 

    Arguments
    =========
    a_mat : Eigen::SparseMatrix<double>
        A in Ax(t+1) = Cx(t) + d.
    c_mat : Eigen::SparseMatrix<double>
        C in Ax(t+1) = Cx(t) + d.
    d_vec : Eigen::VectorXd
        d in Ax(t+1) = Cx(t) + d.
    x_vec : Eigen::VectorXd
        x(t+1) in Ax(t+1) = Cx(t) + d.
    x_last_timestep_vec : Eigen::VectorXd
        x(t) in Ax(t+1) = Cx(t) + d.
    dt : double
        Length of the timestep.

    Returns
    =======
    (none)

    */

    // iterate through each domain covered by a tri3 mesh
    for (int indx_d = 0; indx_d < mesh_physics_ptr->mesh_t3_ptr_vec.size(); indx_d++)
    {

        // subset the mesh
        MeshTri3 *mesh_ptr = mesh_physics_ptr->mesh_t3_ptr_vec[indx_d];

        // subset the velocity
        ScalarTri3 *velocity_x_ptr = velocity_x_field_ptr->mesh_to_scalar_t3_ptr_map[mesh_ptr];
        ScalarTri3 *velocity_y_ptr = velocity_y_field_ptr->mesh_to_scalar_t3_ptr_map[mesh_ptr];

        // indx_r - row in diffusion matrix (diffusion equation)
        // indx_c - column in diffusion matrix (variable)

        // iterate through each diffusion equation
        for (int indx_r = 0; indx_r < boundary_physics_ptr_vec.num_entry; indx_r++)
        {
            
            // subset the boundary conditions and integrals
            BoundaryTri3 *boundary_ptr = boundary_physics_ptr_vec.get_entry(indx_r)->boundary_t3_ptr_vec[indx_d];
            IntegralTri3 *integral_ptr = integral_physics_ptr_vec.get_entry(indx_r)->integral_t3_ptr_vec[indx_d];

            // subset value and derivative and generation coefficient
            VariableFieldGroup *value_field_ptr = value_field_ptr_vec.get_entry(indx_r);
            ScalarTri3 *derivativecoefficient_ptr = derivativecoefficient_field_ptr_vec.get_entry(indx_r)->mesh_to_scalar_t3_ptr_map[mesh_ptr];
            ScalarTri3 *generationcoefficient_ptr = generationcoefficient_field_ptr_vec.get_entry(indx_r)->mesh_to_scalar_t3_ptr_map[mesh_ptr];

            // iterate through each variable
            for (auto &key_value : diffusioncoefficient_field_ptr_mat.get_row(indx_r))
            {

                // subset diffusion coefficient
                int indx_c = key_value.first;
                ScalarTri3 *diffusioncoefficient_ptr = diffusioncoefficient_field_ptr_mat.get_entry(indx_r, indx_c)->mesh_to_scalar_t3_ptr_map[mesh_ptr];

                // determine matrix coefficients for the domain and equation
                matrix_fill_domain_t3(
                    a_mat, c_mat, d_vec,
                    x_vec, x_last_timestep_vec, dt,
                    mesh_ptr, boundary_ptr, integral_ptr,
                    value_field_ptr,
                    derivativecoefficient_ptr, diffusioncoefficient_ptr, velocity_x_ptr, velocity_y_ptr, generationcoefficient_ptr
                );

            }

        }

    }


    // iterate through each domain covered by a quad4 mesh
    for (int indx_d = 0; indx_d < mesh_physics_ptr->mesh_q4_ptr_vec.size(); indx_d++)
    {

        // subset the mesh
        MeshQuad4 *mesh_ptr = mesh_physics_ptr->mesh_q4_ptr_vec[indx_d];

        // subset the velocity
        ScalarQuad4 *velocity_x_ptr = velocity_x_field_ptr->mesh_to_scalar_q4_ptr_map[mesh_ptr];
        ScalarQuad4 *velocity_y_ptr = velocity_y_field_ptr->mesh_to_scalar_q4_ptr_map[mesh_ptr];

        // indx_r - row in diffusion matrix (diffusion equation)
        // indx_c - column in diffusion matrix (variable)

        // iterate through each diffusion equation
        for (int indx_r = 0; indx_r < boundary_physics_ptr_vec.num_entry; indx_r++)
        {
            
            // subset the boundary conditions and integrals
            BoundaryQuad4 *boundary_ptr = boundary_physics_ptr_vec.get_entry(indx_r)->boundary_q4_ptr_vec[indx_d];
            IntegralQuad4 *integral_ptr = integral_physics_ptr_vec.get_entry(indx_r)->integral_q4_ptr_vec[indx_d];

            // subset value and derivative and generation coefficient
            VariableFieldGroup *value_field_ptr = value_field_ptr_vec.get_entry(indx_r);
            ScalarQuad4 *derivativecoefficient_ptr = derivativecoefficient_field_ptr_vec.get_entry(indx_r)->mesh_to_scalar_q4_ptr_map[mesh_ptr];
            ScalarQuad4 *generationcoefficient_ptr = generationcoefficient_field_ptr_vec.get_entry(indx_r)->mesh_to_scalar_q4_ptr_map[mesh_ptr];

            // iterate through each variable
            for (auto &key_value : diffusioncoefficient_field_ptr_mat.get_row(indx_r))
            {

                // subset diffusion coefficient
                int indx_c = key_value.first;
                ScalarQuad4 *diffusioncoefficient_ptr = diffusioncoefficient_field_ptr_mat.get_entry(indx_r, indx_c)->mesh_to_scalar_q4_ptr_map[mesh_ptr];

                // determine matrix coefficients for the domain and equation
                matrix_fill_domain_q4(
                    a_mat, c_mat, d_vec,
                    x_vec, x_last_timestep_vec, dt,
                    mesh_ptr, boundary_ptr, integral_ptr,
                    value_field_ptr,
                    derivativecoefficient_ptr, diffusioncoefficient_ptr, velocity_x_ptr, velocity_y_ptr, generationcoefficient_ptr
                );

            }

        }

    }

}

void PhysicsTransientConvectionDiffusionMulticomponent::matrix_fill_domain_t3
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
    Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
    MeshTri3 *mesh_ptr, BoundaryTri3 *boundary_ptr, IntegralTri3 *integral_ptr,
    VariableFieldGroup *value_field_ptr,
    ScalarTri3 *derivativecoefficient_ptr, ScalarTri3 *diffusioncoefficient_ptr, ScalarTri3 *velocity_x_ptr, ScalarTri3 *velocity_y_ptr, ScalarTri3 *generationcoefficient_ptr
)
{

    // calculate adjustment in starting row
    // e.g., 1st variable starts at start_row, 2nd starts after 1st, etc.
    int adjust_start_row = 0;
    for (auto &value_field_ptr_sub : value_field_ptr_vec.get_vector())
    {
        
        // stop incrementing if dealing with current variable
        if (value_field_ptr_sub == value_field_ptr)
        {
            break;
        }

        // add number of rows occupied by previous variables
        adjust_start_row += value_field_ptr->num_point_field;

    }

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_ptr->num_element_domain; element_did++)
    {

        // get global ID of points around element
        int p0_gid = mesh_ptr->element_p0_gid_vec[element_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[element_did];
        int p2_gid = mesh_ptr->element_p2_gid_vec[element_did];

        // get domain ID of points
        // used for getting properties and integrals
        int p0_did = mesh_ptr->point_gid_to_did_map[p0_gid];
        int p1_did = mesh_ptr->point_gid_to_did_map[p1_gid];
        int p2_did = mesh_ptr->point_gid_to_did_map[p2_gid];
        int did_arr[3] = {p0_did, p1_did, p2_did};

        // get derivative coefficient of points around element
        double dervcoeff_p0 = derivativecoefficient_ptr->point_value_vec[p0_did];
        double dervcoeff_p1 = derivativecoefficient_ptr->point_value_vec[p1_did];
        double dervcoeff_p2 = derivativecoefficient_ptr->point_value_vec[p2_did];
        double dervcoeff_arr[3] = {dervcoeff_p0, dervcoeff_p1, dervcoeff_p2};

        // get diffusion coefficient of points around element
        double diffcoeff_p0 = diffusioncoefficient_ptr->point_value_vec[p0_did];
        double diffcoeff_p1 = diffusioncoefficient_ptr->point_value_vec[p1_did];
        double diffcoeff_p2 = diffusioncoefficient_ptr->point_value_vec[p2_did];
        double diffcoeff_arr[3] = {diffcoeff_p0, diffcoeff_p1, diffcoeff_p2};

        // get x velocity of points around element
        double velx_p0 = velocity_x_ptr->point_value_vec[p0_did];
        double velx_p1 = velocity_x_ptr->point_value_vec[p1_did];
        double velx_p2 = velocity_x_ptr->point_value_vec[p2_did];
        double velx_arr[3] = {velx_p0, velx_p1, velx_p2};

        // get y velocity of points around element
        double vely_p0 = velocity_y_ptr->point_value_vec[p0_did];
        double vely_p1 = velocity_y_ptr->point_value_vec[p1_did];
        double vely_p2 = velocity_y_ptr->point_value_vec[p2_did];
        double vely_arr[3] = {vely_p0, vely_p1, vely_p2};

        // get generation coefficient of points around element
        double gencoeff_p0 = generationcoefficient_ptr->point_value_vec[p0_did];
        double gencoeff_p1 = generationcoefficient_ptr->point_value_vec[p1_did];
        double gencoeff_p2 = generationcoefficient_ptr->point_value_vec[p2_did];
        double gencoeff_arr[3] = {gencoeff_p0, gencoeff_p1, gencoeff_p2};

        // calculate a_mat coefficients
        // matrix row = start_row of test function (physics) + field ID of variable
        // matrix column = start_column of variable + field ID of variable

        // get field ID of temperature points
        // used for getting matrix rows and columns
        int p0_fid = value_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_field_ptr->point_gid_to_fid_map[p1_gid];
        int p2_fid = value_field_ptr->point_gid_to_fid_map[p2_gid];
        int fid_arr[3] = {p0_fid, p1_fid, p2_fid};

        // calculate a_mat and c_mat coefficients
        for (int indx_i = 0; indx_i < 3; indx_i++){
        for (int indx_j = 0; indx_j < 3; indx_j++){
            
            // calculate matrix indices
            int mat_row = start_row + adjust_start_row + fid_arr[indx_i];
            int mat_col = value_field_ptr->start_col + fid_arr[indx_j];

            // calculate velocity derivative
            double integral_Ni_Nj_dvelx_dx = 0;
            double integral_Ni_Nj_dvely_dy = 0;
            for (int indx_k = 0; indx_k < 3; indx_k++){
                integral_Ni_Nj_dvelx_dx += velx_arr[indx_k]*integral_ptr->integral_Ni_Nj_derivative_Nk_x_vec[element_did][indx_i][indx_j][indx_k];
                integral_Ni_Nj_dvely_dy += vely_arr[indx_k]*integral_ptr->integral_Ni_Nj_derivative_Nk_y_vec[element_did][indx_i][indx_j][indx_k];
            }

            // calculate a_mat coefficients
            a_mat.coeffRef(mat_row, mat_col) += (
                (dervcoeff_arr[indx_i]/dt)*integral_ptr->integral_Ni_Nj_vec[element_did][indx_i][indx_j] +
                diffcoeff_arr[indx_i]*integral_ptr->integral_div_Ni_dot_div_Nj_vec[element_did][indx_i][indx_j] +
                velx_arr[indx_i]*integral_ptr->integral_Ni_derivative_Nj_x_vec[element_did][indx_i][indx_j] +
                vely_arr[indx_i]*integral_ptr->integral_Ni_derivative_Nj_y_vec[element_did][indx_i][indx_j] +
                integral_Ni_Nj_dvelx_dx + integral_Ni_Nj_dvely_dy
            );

            // calculate c_mat coefficients
            c_mat.coeffRef(mat_row, mat_col) += (dervcoeff_arr[indx_i]/dt)*integral_ptr->integral_Ni_Nj_vec[element_did][indx_i][indx_j];

        }}

        // calculate d_vec coefficients
        for (int indx_i = 0; indx_i < 3; indx_i++)
        {
            int mat_row = start_row + adjust_start_row + fid_arr[indx_i];
            d_vec.coeffRef(mat_row) += gencoeff_arr[indx_i]*integral_ptr->integral_Ni_vec[element_did][indx_i];
        }

    }

    // iterate for each flux boundary element
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_element_flux_domain; boundary_id++)
    {

        // get global ID of element
        int ea_gid = boundary_ptr->element_flux_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int ea_did = mesh_ptr->element_gid_to_did_map[ea_gid];

        // get global ID of points
        int p0_gid = mesh_ptr->element_p0_gid_vec[ea_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[ea_did];
        int p2_gid = mesh_ptr->element_p2_gid_vec[ea_did];

        // get local ID of point where boundary is applied
        int pa_lid = boundary_ptr->element_flux_pa_lid_vec[boundary_id];  // 0 to 2
        int pb_lid = boundary_ptr->element_flux_pb_lid_vec[boundary_id];  // 0 to 2

        // get edge where boundary is applied
        int helper_num = pa_lid + pb_lid + 1;
        int boundary_key = (helper_num*helper_num - helper_num % 2)/4 + std::min(pa_lid, pb_lid);

        // identify boundary type
        int config_id = boundary_ptr->element_flux_boundaryconfig_id_vec[boundary_id];
        BoundaryConfigStruct boundaryconfig = boundary_ptr->boundaryconfig_vec[config_id];

        // get field ID of value points
        // used for getting matrix rows and columns
        int p0_fid = value_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_field_ptr->point_gid_to_fid_map[p1_gid];
        int p2_fid = value_field_ptr->point_gid_to_fid_map[p2_gid];
        int fid_arr[3] = {p0_fid, p1_fid, p2_fid};

        // get domain ID of points
        // used for getting properties and integrals
        int p0_did = mesh_ptr->point_gid_to_did_map[p0_gid];
        int p1_did = mesh_ptr->point_gid_to_did_map[p1_gid];
        int p2_did = mesh_ptr->point_gid_to_did_map[p2_gid];

        // get x velocity of points around element
        double velx_p0 = velocity_x_ptr->point_value_vec[p0_did];
        double velx_p1 = velocity_x_ptr->point_value_vec[p1_did];
        double velx_p2 = velocity_x_ptr->point_value_vec[p2_did];

        // get y velocity of points around element
        double vely_p0 = velocity_y_ptr->point_value_vec[p0_did];
        double vely_p1 = velocity_y_ptr->point_value_vec[p1_did];
        double vely_p2 = velocity_y_ptr->point_value_vec[p2_did];

        // get dot product of velocity and normal
        double norm_x = integral_ptr->boundary_normal_x_map[ea_did][boundary_key];
        double norm_y = integral_ptr->boundary_normal_y_map[ea_did][boundary_key];
        double vel_dot_normal_p0 = velx_p0*norm_x + vely_p0*norm_y;
        double vel_dot_normal_p1 = velx_p1*norm_x + vely_p1*norm_y;
        double vel_dot_normal_p2 = velx_p2*norm_x + vely_p2*norm_y;             
        double vel_dot_normal_arr[3] = {vel_dot_normal_p0, vel_dot_normal_p1, vel_dot_normal_p2};

        // apply boundary condition
        if (boundaryconfig.type_str == "neumann")
        {
            
            // set a_mat and d_vec
            // -1 values indicate invalid points
            if (pa_lid != -1)
            {
                
                // constant part - add terms to b vector
                int mat_row = start_row + adjust_start_row + fid_arr[pa_lid];
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0] * integral_ptr->boundary_integral_Ni_map[ea_did][boundary_key][pa_lid];

                // linear part - iterate over all test functions in element
                for (int indx_j = 0; indx_j < 3; indx_j++){
                    int mat_col = value_field_ptr->start_col + fid_arr[indx_j];  
                    a_mat.coeffRef(mat_row, mat_col) += -vel_dot_normal_arr[pa_lid] * integral_ptr->boundary_integral_Ni_Nj_map[ea_did][boundary_key][pa_lid][indx_j];
                }
                
            }
            if (pb_lid != -1)
            {

                // constant part - add terms to b vector
                int mat_row = start_row + adjust_start_row + fid_arr[pb_lid];
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0] * integral_ptr->boundary_integral_Ni_map[ea_did][boundary_key][pb_lid];

                // linear part - iterate over all test functions in element
                for (int indx_j = 0; indx_j < 3; indx_j++){
                    int mat_col = value_field_ptr->start_col + fid_arr[indx_j];  
                    a_mat.coeffRef(mat_row, mat_col) += -vel_dot_normal_arr[pb_lid] * integral_ptr->boundary_integral_Ni_Nj_map[ea_did][boundary_key][pb_lid][indx_j];
                }

            }

        }
        else if (boundaryconfig.type_str == "robin")
        {

            // set a_mat and d_vec
            // -1 values indicate invalid points
            if (pa_lid != -1)
            {
                
                // constant part - add terms to b vector
                int mat_row = start_row + adjust_start_row + fid_arr[pa_lid];
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0] * integral_ptr->boundary_integral_Ni_map[ea_did][boundary_key][pa_lid];

                // linear part - iterate over all test functions in element
                for (int indx_j = 0; indx_j < 3; indx_j++){
                    int mat_col = value_field_ptr->start_col + fid_arr[indx_j];  
                    a_mat.coeffRef(mat_row, mat_col) += -(boundaryconfig.parameter_vec[1] + vel_dot_normal_arr[pa_lid]) * integral_ptr->boundary_integral_Ni_Nj_map[ea_did][boundary_key][pa_lid][indx_j];
                }
                
            }
            if (pb_lid != -1)
            {

                // constant part - add terms to b vector
                int mat_row = start_row + adjust_start_row + fid_arr[pb_lid];
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0] * integral_ptr->boundary_integral_Ni_map[ea_did][boundary_key][pb_lid];

                // linear part - iterate over all test functions in element
                for (int indx_j = 0; indx_j < 3; indx_j++){
                    int mat_col = value_field_ptr->start_col + fid_arr[indx_j];  
                    a_mat.coeffRef(mat_row, mat_col) += -(boundaryconfig.parameter_vec[1] + vel_dot_normal_arr[pb_lid]) * integral_ptr->boundary_integral_Ni_Nj_map[ea_did][boundary_key][pb_lid][indx_j];
                }

            }

        }

    }

    // clear rows with value boundary elements
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_element_value_domain; boundary_id++)
    {

        // get global ID of element
        int ea_gid = boundary_ptr->element_value_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int ea_did = mesh_ptr->element_gid_to_did_map[ea_gid];

        // get global ID of points
        int p0_gid = mesh_ptr->element_p0_gid_vec[ea_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[ea_did];
        int p2_gid = mesh_ptr->element_p2_gid_vec[ea_did];

        // get local ID of point where boundary is applied
        int pa_lid = boundary_ptr->element_value_pa_lid_vec[boundary_id];  // 0 to 2
        int pb_lid = boundary_ptr->element_value_pb_lid_vec[boundary_id];  // 0 to 2

        // get field ID of value points
        // used for getting matrix rows and columns
        int p0_fid = value_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_field_ptr->point_gid_to_fid_map[p1_gid];
        int p2_fid = value_field_ptr->point_gid_to_fid_map[p2_gid];
        int fid_arr[3] = {p0_fid, p1_fid, p2_fid};

        // erase entire row
        // -1 values indicate invalid points
        if (pa_lid != -1)
        {
            int mat_row = start_row + adjust_start_row + fid_arr[pa_lid];
            a_mat.row(mat_row) *= 0.;
            c_mat.row(mat_row) *= 0.;
            d_vec.coeffRef(mat_row) = 0.;
        }
        if (pb_lid != -1)
        {
            int mat_row = start_row + adjust_start_row + fid_arr[pa_lid];
            a_mat.row(mat_row) *= 0.;
            c_mat.row(mat_row) *= 0.;
            d_vec.coeffRef(mat_row) = 0.;
        }

    }

    // iterate for each value boundary element
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_element_value_domain; boundary_id++)
    {

        // get global ID of element
        int ea_gid = boundary_ptr->element_value_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int ea_did = mesh_ptr->element_gid_to_did_map[ea_gid];

        // get global ID of points
        int p0_gid = mesh_ptr->element_p0_gid_vec[ea_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[ea_did];
        int p2_gid = mesh_ptr->element_p2_gid_vec[ea_did];

        // get local ID of point where boundary is applied
        int pa_lid = boundary_ptr->element_value_pa_lid_vec[boundary_id];  // 0 to 2
        int pb_lid = boundary_ptr->element_value_pb_lid_vec[boundary_id];  // 0 to 2
        
        // identify boundary type
        int config_id = boundary_ptr->element_value_boundaryconfig_id_vec[boundary_id];
        BoundaryConfigStruct boundaryconfig = boundary_ptr->boundaryconfig_vec[config_id];

        // get field ID of value points
        // used for getting matrix rows and columns
        int p0_fid = value_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_field_ptr->point_gid_to_fid_map[p1_gid];
        int p2_fid = value_field_ptr->point_gid_to_fid_map[p2_gid];
        int fid_arr[3] = {p0_fid, p1_fid, p2_fid};

        // apply boundary condition
        if (boundaryconfig.type_str == "dirichlet")
        {

            // set a_mat and d_vec
            // -1 values indicate invalid points
            if (pa_lid != -1)
            {
                int mat_row = start_row + adjust_start_row + fid_arr[pa_lid];
                int mat_col = value_field_ptr->start_col + fid_arr[pa_lid];
                a_mat.coeffRef(mat_row, mat_col) += 1.;
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0];
            }
            if (pb_lid != -1)
            {
                int mat_row = start_row + adjust_start_row + fid_arr[pb_lid];
                int mat_col = value_field_ptr->start_col + fid_arr[pb_lid];
                a_mat.coeffRef(mat_row, mat_col) += 1.;
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0];
            }

        }

    }

}

void PhysicsTransientConvectionDiffusionMulticomponent::matrix_fill_domain_q4
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
    Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
    MeshQuad4 *mesh_ptr, BoundaryQuad4 *boundary_ptr, IntegralQuad4 *integral_ptr,
    VariableFieldGroup *value_field_ptr,
    ScalarQuad4 *derivativecoefficient_ptr, ScalarQuad4 *diffusioncoefficient_ptr, ScalarQuad4 *velocity_x_ptr, ScalarQuad4 *velocity_y_ptr, ScalarQuad4 *generationcoefficient_ptr
)
{

    // calculate adjustment in starting row
    // e.g., 1st variable starts at start_row, 2nd starts after 1st, etc.
    int adjust_start_row = 0;
    for (auto &value_field_ptr_sub : value_field_ptr_vec.get_vector())
    {
        
        // stop incrementing if dealing with current variable
        if (value_field_ptr_sub == value_field_ptr)
        {
            break;
        }

        // add number of rows occupied by previous variables
        adjust_start_row += value_field_ptr->num_point_field;

    }

    // iterate for each domain element
    for (int element_did = 0; element_did < mesh_ptr->num_element_domain; element_did++)
    {

        // get global ID of points around element
        int p0_gid = mesh_ptr->element_p0_gid_vec[element_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[element_did];
        int p2_gid = mesh_ptr->element_p2_gid_vec[element_did];
        int p3_gid = mesh_ptr->element_p3_gid_vec[element_did];

        // get domain ID of points
        // used for getting properties and integrals
        int p0_did = mesh_ptr->point_gid_to_did_map[p0_gid];
        int p1_did = mesh_ptr->point_gid_to_did_map[p1_gid];
        int p2_did = mesh_ptr->point_gid_to_did_map[p2_gid];
        int p3_did = mesh_ptr->point_gid_to_did_map[p3_gid];
        int did_arr[4] = {p0_did, p1_did, p2_did, p3_did};

        // get derivative coefficient of points around element
        double dervcoeff_p0 = derivativecoefficient_ptr->point_value_vec[p0_did];
        double dervcoeff_p1 = derivativecoefficient_ptr->point_value_vec[p1_did];
        double dervcoeff_p2 = derivativecoefficient_ptr->point_value_vec[p2_did];
        double dervcoeff_p3 = derivativecoefficient_ptr->point_value_vec[p3_did];
        double dervcoeff_arr[4] = {dervcoeff_p0, dervcoeff_p1, dervcoeff_p2, dervcoeff_p3};

        // get diffusion coefficient of points around element
        double diffcoeff_p0 = diffusioncoefficient_ptr->point_value_vec[p0_did];
        double diffcoeff_p1 = diffusioncoefficient_ptr->point_value_vec[p1_did];
        double diffcoeff_p2 = diffusioncoefficient_ptr->point_value_vec[p2_did];
        double diffcoeff_p3 = diffusioncoefficient_ptr->point_value_vec[p3_did];
        double diffcoeff_arr[4] = {diffcoeff_p0, diffcoeff_p1, diffcoeff_p2, diffcoeff_p3};

        // get x velocity of points around element
        double velx_p0 = velocity_x_ptr->point_value_vec[p0_did];
        double velx_p1 = velocity_x_ptr->point_value_vec[p1_did];
        double velx_p2 = velocity_x_ptr->point_value_vec[p2_did];
        double velx_p3 = velocity_x_ptr->point_value_vec[p3_did];
        double velx_arr[4] = {velx_p0, velx_p1, velx_p2};

        // get y velocity of points around element
        double vely_p0 = velocity_y_ptr->point_value_vec[p0_did];
        double vely_p1 = velocity_y_ptr->point_value_vec[p1_did];
        double vely_p2 = velocity_y_ptr->point_value_vec[p2_did];
        double vely_p3 = velocity_y_ptr->point_value_vec[p3_did];
        double vely_arr[4] = {vely_p0, vely_p1, vely_p2, vely_p3};

        // get generation coefficient of points around element
        double gencoeff_p0 = generationcoefficient_ptr->point_value_vec[p0_did];
        double gencoeff_p1 = generationcoefficient_ptr->point_value_vec[p1_did];
        double gencoeff_p2 = generationcoefficient_ptr->point_value_vec[p2_did];
        double gencoeff_p3 = generationcoefficient_ptr->point_value_vec[p3_did];
        double gencoeff_arr[4] = {gencoeff_p0, gencoeff_p1, gencoeff_p2, gencoeff_p3};

        // calculate a_mat coefficients
        // matrix row = start_row of test function (physics) + field ID of variable
        // matrix column = start_column of variable + field ID of variable

        // get field ID of value points
        // used for getting matrix rows and columns
        int p0_fid = value_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_field_ptr->point_gid_to_fid_map[p1_gid];
        int p2_fid = value_field_ptr->point_gid_to_fid_map[p2_gid];
        int p3_fid = value_field_ptr->point_gid_to_fid_map[p3_gid];
        int fid_arr[4] = {p0_fid, p1_fid, p2_fid, p3_fid};

        // calculate a_mat and c_mat coefficients
        for (int indx_i = 0; indx_i < 4; indx_i++){
        for (int indx_j = 0; indx_j < 4; indx_j++){
            
            // calculate matrix indices
            int mat_row = start_row + adjust_start_row + fid_arr[indx_i];
            int mat_col = value_field_ptr->start_col + fid_arr[indx_j];

            // calculate velocity derivative
            double integral_Ni_Nj_dvelx_dx = 0;
            double integral_Ni_Nj_dvely_dy = 0;
            for (int indx_k = 0; indx_k < 4; indx_k++){
                integral_Ni_Nj_dvelx_dx += velx_arr[indx_k]*integral_ptr->integral_Ni_Nj_derivative_Nk_x_vec[element_did][indx_i][indx_j][indx_k];
                integral_Ni_Nj_dvely_dy += vely_arr[indx_k]*integral_ptr->integral_Ni_Nj_derivative_Nk_y_vec[element_did][indx_i][indx_j][indx_k];
            }

            // calculate a_mat coefficients
            a_mat.coeffRef(mat_row, mat_col) += (
                (dervcoeff_arr[indx_i]/dt)*integral_ptr->integral_Ni_Nj_vec[element_did][indx_i][indx_j] +
                diffcoeff_arr[indx_i]*integral_ptr->integral_div_Ni_dot_div_Nj_vec[element_did][indx_i][indx_j] +
                velx_arr[indx_i]*integral_ptr->integral_Ni_derivative_Nj_x_vec[element_did][indx_i][indx_j] +
                vely_arr[indx_i]*integral_ptr->integral_Ni_derivative_Nj_y_vec[element_did][indx_i][indx_j] +
                integral_Ni_Nj_dvelx_dx + integral_Ni_Nj_dvely_dy
            );

            // calculate c_mat coefficients
            c_mat.coeffRef(mat_row, mat_col) += (dervcoeff_arr[indx_i]/dt)*integral_ptr->integral_Ni_Nj_vec[element_did][indx_i][indx_j];

        }}

        // calculate d_vec coefficients
        for (int indx_i = 0; indx_i < 4; indx_i++)
        {
            int mat_row = start_row + adjust_start_row + fid_arr[indx_i];
            d_vec.coeffRef(mat_row) += gencoeff_arr[indx_i]*integral_ptr->integral_Ni_vec[element_did][indx_i];
        }

    }

    // iterate for each flux boundary element
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_element_flux_domain; boundary_id++)
    {

        // get global ID of element
        int ea_gid = boundary_ptr->element_flux_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int ea_did = mesh_ptr->element_gid_to_did_map[ea_gid];

        // get global ID of points
        int p0_gid = mesh_ptr->element_p0_gid_vec[ea_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[ea_did];
        int p2_gid = mesh_ptr->element_p2_gid_vec[ea_did];
        int p3_gid = mesh_ptr->element_p3_gid_vec[ea_did];

        // get local ID of point where boundary is applied
        int pa_lid = boundary_ptr->element_flux_pa_lid_vec[boundary_id];  // 0 to 3
        int pb_lid = boundary_ptr->element_flux_pb_lid_vec[boundary_id];  // 0 to 3

        // get edge where boundary is applied
        int helper_num = pa_lid + pb_lid + 1;
        int boundary_key = (helper_num*helper_num - helper_num % 2)/4 + std::min(pa_lid, pb_lid);

        // identify boundary type
        int config_id = boundary_ptr->element_flux_boundaryconfig_id_vec[boundary_id];
        BoundaryConfigStruct boundaryconfig = boundary_ptr->boundaryconfig_vec[config_id];

        // get field ID of value points
        // used for getting matrix rows and columns
        int p0_fid = value_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_field_ptr->point_gid_to_fid_map[p1_gid];
        int p2_fid = value_field_ptr->point_gid_to_fid_map[p2_gid];
        int p3_fid = value_field_ptr->point_gid_to_fid_map[p3_gid];
        int fid_arr[4] = {p0_fid, p1_fid, p2_fid, p3_fid};

        // get domain ID of points
        // used for getting properties and integrals
        int p0_did = mesh_ptr->point_gid_to_did_map[p0_gid];
        int p1_did = mesh_ptr->point_gid_to_did_map[p1_gid];
        int p2_did = mesh_ptr->point_gid_to_did_map[p2_gid];
        int p3_did = mesh_ptr->point_gid_to_did_map[p3_gid];

        // get x velocity of points around element
        double velx_p0 = velocity_x_ptr->point_value_vec[p0_did];
        double velx_p1 = velocity_x_ptr->point_value_vec[p1_did];
        double velx_p2 = velocity_x_ptr->point_value_vec[p2_did];
        double velx_p3 = velocity_x_ptr->point_value_vec[p3_did];

        // get y velocity of points around element
        double vely_p0 = velocity_y_ptr->point_value_vec[p0_did];
        double vely_p1 = velocity_y_ptr->point_value_vec[p1_did];
        double vely_p2 = velocity_y_ptr->point_value_vec[p2_did];
        double vely_p3 = velocity_y_ptr->point_value_vec[p3_did];

        // get dot product of velocity and normal
        double norm_x = integral_ptr->boundary_normal_x_map[ea_did][boundary_key];
        double norm_y = integral_ptr->boundary_normal_y_map[ea_did][boundary_key];
        double vel_dot_normal_p0 = velx_p0*norm_x + vely_p0*norm_y;
        double vel_dot_normal_p1 = velx_p1*norm_x + vely_p1*norm_y;
        double vel_dot_normal_p2 = velx_p2*norm_x + vely_p2*norm_y;
        double vel_dot_normal_p3 = velx_p3*norm_x + vely_p3*norm_y;
        double vel_dot_normal_arr[4] = {vel_dot_normal_p0, vel_dot_normal_p1, vel_dot_normal_p2, vel_dot_normal_p3};

        // apply boundary condition
        if (boundaryconfig.type_str == "neumann")
        {
            
            // set a_mat and d_vec
            // -1 values indicate invalid points
            if (pa_lid != -1)
            {
                
                // constant part - add terms to b vector
                int mat_row = start_row + adjust_start_row + fid_arr[pa_lid];
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0] * integral_ptr->boundary_integral_Ni_map[ea_did][boundary_key][pa_lid];

                // linear part - iterate over all test functions in element
                for (int indx_j = 0; indx_j < 4; indx_j++){
                    int mat_col = value_field_ptr->start_col + fid_arr[indx_j];  
                    a_mat.coeffRef(mat_row, mat_col) += -vel_dot_normal_arr[pa_lid] * integral_ptr->boundary_integral_Ni_Nj_map[ea_did][boundary_key][pa_lid][indx_j];
                }

            }
            if (pb_lid != -1)
            {

                // constant part - add terms to b vector
                int mat_row = start_row + adjust_start_row + fid_arr[pb_lid];
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0] * integral_ptr->boundary_integral_Ni_map[ea_did][boundary_key][pb_lid];

                // linear part - iterate over all test functions in element
                for (int indx_j = 0; indx_j < 4; indx_j++){
                    int mat_col = value_field_ptr->start_col + fid_arr[indx_j];  
                    a_mat.coeffRef(mat_row, mat_col) += -vel_dot_normal_arr[pb_lid] * integral_ptr->boundary_integral_Ni_Nj_map[ea_did][boundary_key][pb_lid][indx_j];
                }

            }

        }
        else if (boundaryconfig.type_str == "robin")
        {

            // set a_mat and d_vec
            // -1 values indicate invalid points
            if (pa_lid != -1)
            {
                
                // constant part - add terms to b vector
                int mat_row = start_row + adjust_start_row + fid_arr[pa_lid];
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0] * integral_ptr->boundary_integral_Ni_map[ea_did][boundary_key][pa_lid];

                // linear part - iterate over all test functions in element
                for (int indx_j = 0; indx_j < 4; indx_j++){
                    int mat_col = value_field_ptr->start_col + fid_arr[indx_j];  
                    a_mat.coeffRef(mat_row, mat_col) += -(vel_dot_normal_arr[pa_lid] + boundaryconfig.parameter_vec[1]) * integral_ptr->boundary_integral_Ni_Nj_map[ea_did][boundary_key][pa_lid][indx_j];
                }
                
            }
            if (pb_lid != -1)
            {

                // constant part - add terms to b vector
                int mat_row = start_row + adjust_start_row + fid_arr[pb_lid];
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0] * integral_ptr->boundary_integral_Ni_map[ea_did][boundary_key][pb_lid];

                // linear part - iterate over all test functions in element
                for (int indx_j = 0; indx_j < 4; indx_j++){
                    int mat_col = value_field_ptr->start_col + fid_arr[indx_j];  
                    a_mat.coeffRef(mat_row, mat_col) += -(vel_dot_normal_arr[pb_lid] + boundaryconfig.parameter_vec[1]) * integral_ptr->boundary_integral_Ni_Nj_map[ea_did][boundary_key][pb_lid][indx_j];
                }

            }

        }

    }

    // clear rows with value boundary elements
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_element_value_domain; boundary_id++)
    {

        // get global ID of element
        int ea_gid = boundary_ptr->element_value_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int ea_did = mesh_ptr->element_gid_to_did_map[ea_gid];

        // get global ID of points
        int p0_gid = mesh_ptr->element_p0_gid_vec[ea_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[ea_did];
        int p2_gid = mesh_ptr->element_p2_gid_vec[ea_did];
        int p3_gid = mesh_ptr->element_p3_gid_vec[ea_did];

        // get local ID of point where boundary is applied
        int pa_lid = boundary_ptr->element_value_pa_lid_vec[boundary_id];  // 0 to 3
        int pb_lid = boundary_ptr->element_value_pb_lid_vec[boundary_id];  // 0 to 3

        // get field ID of value points
        // used for getting matrix rows and columns
        int p0_fid = value_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_field_ptr->point_gid_to_fid_map[p1_gid];
        int p2_fid = value_field_ptr->point_gid_to_fid_map[p2_gid];
        int p3_fid = value_field_ptr->point_gid_to_fid_map[p3_gid];
        int fid_arr[4] = {p0_fid, p1_fid, p2_fid, p3_fid};

        // erase entire row
        // -1 values indicate invalid points
        if (pa_lid != -1)
        {
            int mat_row = start_row + adjust_start_row + fid_arr[pa_lid];
            a_mat.row(mat_row) *= 0.;
            c_mat.row(mat_row) *= 0.;
            d_vec.coeffRef(mat_row) = 0.;
        }
        if (pb_lid != -1)
        {
            int mat_row = start_row + adjust_start_row + fid_arr[pa_lid];
            a_mat.row(mat_row) *= 0.;
            c_mat.row(mat_row) *= 0.;
            d_vec.coeffRef(mat_row) = 0.;
        }

    }

    // iterate for each value boundary element
    for (int boundary_id = 0; boundary_id < boundary_ptr->num_element_value_domain; boundary_id++)
    {

        // get global ID of element
        int ea_gid = boundary_ptr->element_value_gid_vec[boundary_id];

        // get domain ID of element
        // used for getting global ID of points
        int ea_did = mesh_ptr->element_gid_to_did_map[ea_gid];

        // get global ID of points
        int p0_gid = mesh_ptr->element_p0_gid_vec[ea_did];
        int p1_gid = mesh_ptr->element_p1_gid_vec[ea_did];
        int p2_gid = mesh_ptr->element_p2_gid_vec[ea_did];
        int p3_gid = mesh_ptr->element_p3_gid_vec[ea_did];

        // get local ID of point where boundary is applied
        int pa_lid = boundary_ptr->element_value_pa_lid_vec[boundary_id];  // 0 to 3
        int pb_lid = boundary_ptr->element_value_pb_lid_vec[boundary_id];  // 0 to 3
        
        // identify boundary type
        int config_id = boundary_ptr->element_value_boundaryconfig_id_vec[boundary_id];
        BoundaryConfigStruct boundaryconfig = boundary_ptr->boundaryconfig_vec[config_id];

        // get field ID of value points
        // used for getting matrix rows and columns
        int p0_fid = value_field_ptr->point_gid_to_fid_map[p0_gid];
        int p1_fid = value_field_ptr->point_gid_to_fid_map[p1_gid];
        int p2_fid = value_field_ptr->point_gid_to_fid_map[p2_gid];
        int p3_fid = value_field_ptr->point_gid_to_fid_map[p3_gid];
        int fid_arr[4] = {p0_fid, p1_fid, p2_fid, p3_fid};

        // apply boundary condition
        if (boundaryconfig.type_str == "dirichlet")
        {

            // set a_mat and d_vec
            // -1 values indicate invalid points
            if (pa_lid != -1)
            {
                int mat_row = start_row + adjust_start_row + fid_arr[pa_lid];
                int mat_col = value_field_ptr->start_col + fid_arr[pa_lid];
                a_mat.coeffRef(mat_row, mat_col) += 1.;
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0];
            }
            if (pb_lid != -1)
            {
                int mat_row = start_row + adjust_start_row + fid_arr[pb_lid];
                int mat_col = value_field_ptr->start_col + fid_arr[pb_lid];
                a_mat.coeffRef(mat_row, mat_col) += 1.;
                d_vec.coeffRef(mat_row) += boundaryconfig.parameter_vec[0];
            }

        }

    }

}

void PhysicsTransientConvectionDiffusionMulticomponent::set_start_row(int start_row_in)
{
    /*

    Sets the starting row in A and b where entries are filled up.

    Arguments
    =========
    start_row_in : int
        Starting row in A and b.

    Returns
    =======
    (none)

    */

    start_row = start_row_in;

}

int PhysicsTransientConvectionDiffusionMulticomponent::get_start_row()
{
    /*

    Returns the starting row.

    Arguments
    =========
    (none)

    Returns
    =======
    start_row : int
        Starting row in A and b.

    */

    return start_row;

}

std::vector<VariableFieldGroup*> PhysicsTransientConvectionDiffusionMulticomponent::get_variable_field_ptr_vec()
{
    /*

    Returns the vector containing pointers to VariableFieldGroup objects tied to this physics.

    Arguments
    =========
    (none)

    Returns
    =======
    variable_field_ptr : vector<VariableFieldGroup*>
        Vector containing pointers to VariableFieldGroup objects.

    */

    return variable_field_ptr_vec;

}

#endif
