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
#include "boundary_tri3.hpp"
#include "container_typedef.hpp"
#include "integral_physicsgroup.hpp"
#include "integral_quad4.hpp"
#include "integral_tri3.hpp"
#include "matrixequation_steady.hpp"
#include "matrixequation_transient.hpp"
#include "mesh_physicsgroup.hpp"
#include "mesh_quad4.hpp"
#include "mesh_tri3.hpp"
#include "physicssteady_base.hpp"
#include "physicssteady_convectiondiffusion.hpp"
#include "physicssteady_diffusion.hpp"
#include "physicstransient_base.hpp"
#include "physicstransient_convectiondiffusion.hpp"
#include "physicstransient_diffusion.hpp"
#include "scalar_fieldgroup.hpp"
#include "scalar_quad4.hpp"
#include "variable_fieldgroup.hpp"
#include "variable_quad4.hpp"
