/*
####################################################
####################################################
###   __  __ ____ ____  ____                     ###
###  |  \/  |  _ \___ \|  _ \     Multi-purpose  ###
###  | |\/| | |_) |__) | |_) |    Multiphysics   ###
###  | |  | |  __// __/|  __/     Program        ###
###  |_|  |_|_|  |_____|_|        (FEM 2D)       ###
###                                              ###
####################################################
####################################################
*/

#include "boundary_physicsgroup.hpp"
#include "boundary_quad4.hpp"
#include "container_typedef.hpp"
#include "initialize_mesh_quad4_csv.hpp"
#include "integral_physicsgroup.hpp"
#include "integral_quad4.hpp"
#include "matrixequation_steady.hpp"
#include "mesh_physicsgroup.hpp"
#include "mesh_quad4.hpp"
#include "physicssteady_base.hpp"
#include "physicssteady_diffusion.hpp"
#include "scalar_fieldgroup.hpp"
#include "scalar_quad4.hpp"
#include "variable_fieldgroup.hpp"
#include "variable_quad4.hpp"

// #include "boundary_line2.hpp"
// #include "boundary_physicsgroup_vector.hpp"
// #include "boundary_physicsgroup.hpp"
// #include "initialize_mesh_line2_csv.hpp"
// #include "integral_line2.hpp"
// #include "integral_physicsgroup.hpp"
// #include "matrixequation_steady.hpp"
// #include "matrixequation_transient.hpp"
// #include "mesh_line2.hpp"
// #include "mesh_physicsgroup.hpp"
// #include "physicssteady_base.hpp"
// #include "physicssteady_convectiondiffusion_multicomponent.hpp"
// #include "physicssteady_convectiondiffusion.hpp"
// #include "physicssteady_diffusion_multicomponent.hpp"
// #include "physicssteady_diffusion.hpp"
// #include "physicstransient_base.hpp"
// #include "physicstransient_convectiondiffusion_multicomponent.hpp"
// #include "physicstransient_convectiondiffusion.hpp"
// #include "physicstransient_diffusion_multicomponent.hpp"
// #include "physicstransient_diffusion.hpp"
// #include "scalar_fieldgroup_matrix.hpp"
// #include "scalar_fieldgroup_vector.hpp"
// #include "scalar_fieldgroup.hpp"
// #include "scalar_line2.hpp"
// #include "variable_fieldgroup_vector.hpp"
// #include "variable_fieldgroup.hpp"
// #include "variable_line2.hpp"
