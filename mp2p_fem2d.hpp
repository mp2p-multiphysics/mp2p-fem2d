/*
###################################################
###################################################
###   __  __ ____ ____  ____                    ###
###  |  \/  |  _ \___ \|  _ \     Multipurpose  ###
###  | |\/| | |_) |__) | |_) |    Multiphysics  ###
###  | |  | |  __// __/|  __/     Package       ###
###  |_|  |_|_|  |_____|_|        (FEM 2D)      ###
###                                             ###
###################################################
###################################################
*/

#include "boundary_group.hpp"
#include "boundary_quad4.hpp"
#include "boundary_tri3.hpp"
#include "container_typedef.hpp"
#include "domain_group.hpp"
#include "domain_quad4.hpp"
#include "domain_tri3.hpp"
#include "integral_group.hpp"
#include "integral_quad4.hpp"
#include "integral_tri3.hpp"
#include "matrixequation_steady.hpp"
#include "matrixequation_transient.hpp"
#include "physicssteady_base.hpp"
#include "physicssteady_diffusion.hpp"
#include "physicstransient_base.hpp"
#include "physicstransient_diffusion.hpp"
#include "scalar_group.hpp"
#include "scalar_quad4.hpp"
#include "scalar_tri3.hpp"
#include "variable_group.hpp"
#include "variable_quad4.hpp"
#include "variable_tri3.hpp"


// #include "boundary_group.hpp"
// #include "boundary_line2.hpp"
// #include "container_typedef.hpp"
// #include "domain_group.hpp"
// #include "domain_line2.hpp"
// #include "integral_group.hpp"
// #include "integral_line2.hpp"
// #include "matrixequation_steady.hpp"
// #include "matrixequation_transient.hpp"
// #include "physicssteady_base.hpp"
// #include "physicssteady_convectiondiffusion.hpp"
// #include "physicssteady_diffusion.hpp"
// #include "physicstransient_base.hpp"
// #include "physicstransient_convectiondiffusion.hpp"
// #include "physicstransient_diffusion.hpp"
// #include "scalar_group.hpp"
// #include "scalar_line2.hpp"
// #include "variable_group.hpp"
// #include "variable_line2.hpp"
// 