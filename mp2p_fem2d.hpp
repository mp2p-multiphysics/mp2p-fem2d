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

#include "container_typedef.hpp"
#include "domain_1d.hpp"
#include "domain_2d.hpp"
#include "integral_1d.hpp"
#include "integral_2d.hpp"
#include "matrixequation_steady.hpp"
#include "matrixequation_transient.hpp"
#include "physicssteady_base.hpp"
#include "physicssteady_diffusion.hpp"
#include "physicstransient_base.hpp"
#include "physicstransient_diffusion.hpp"
#include "scalar_1d.hpp"
#include "scalar_2d.hpp"
#include "variable_2d.hpp"
#include "variable_group.hpp"


// #include "boundary_group.hpp"
// #include "boundary_unit.hpp"
// #include "boundaryintegral_group.hpp"
// #include "boundaryintegral_unit.hpp"
// #include "container_typedef.hpp"
// #include "domain_group.hpp"
// #include "domain_unit.hpp"
// #include "domainintegral_group.hpp"
// #include "domainintegral_unit.hpp"
// #include "matrixequation_steady.hpp"
// #include "matrixequation_transient.hpp"
// #include "physicssteady_base.hpp"
// #include "physicssteady_convectiondiffusion.hpp"
// #include "physicssteady_diffusion.hpp"
// #include "physicstransient_base.hpp"
// #include "physicstransient_convectiondiffusion.hpp"
// #include "physicstransient_diffusion.hpp"
// #include "scalar_group.hpp"
// #include "scalar_unit.hpp"
// #include "variable_group.hpp"
// #include "variable_unit.hpp"
