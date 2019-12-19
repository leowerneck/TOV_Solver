/* .-----------------------------------------.
 * |  Copyright (c) 2019, Leonardo Werneck   |
 * | Licensed under the BSD 2-Clause License |
 * .-----------------------------------------.
 */

/* Program     : TOV Solver
 * File        : Polytropic_EOS__struct_initialization.C
 * Author      : Leo Werneck (werneck@if.usp.br)
 * Date        : October 29, 2019
 *
 * Description : This file is dedicated to the initialization
 *               of the eos_struct array.
 *
 * Dependencies: string.h, math.h, & TOV_headers.h
 *
 * Reference(s): Read et al., PRD 79, 124032 (2009) | (https://arxiv.org/pdf/0812.2163.pdf)
 *  
 */

/* Include dependencies */
#include <string>
#include <math.h>
#include "TOV_headers.h"

/* .---------------------.
 * | Function prototypes |
 * .---------------------.
 */
inline void convert_EOS_struct_to_geometrized_units( eos_struct &eos );
inline void populate_Kpoly_PPEOS( eos_struct &eos );
inline void populate_epsIC_PPEOS( eos_struct &eos );
inline int get_EOS_parameters_from_EOSname( string EOSname, eos_struct &eos );

/* Function   : initialize_EOS_struct()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Initializes the EOS struct from the user input
 *
 * Input(s)   : EOSname - String containing the name of the EOS to be initialized into
 *                        the EOS struct. The names and parameters follow those of
 *                        Read et al. (see references on top of file).
 *            : eos     - Struct containing EOS information
  *
 * Outputs(s) : eos     - Struct containing EOS information, fully initialized
 */
inline void initialize_EOS_struct( string EOSname, eos_struct &eos ) {

  /* First initialize the EOS parameters from table III in Read et al.
   * This sets:
   * .--------------------.
   * | eos.Gamma_PPEOS[4] |
   * | eos.Gamma_PPEOS[5] |
   * | eos.Gamma_PPEOS[6] | Piecewise Polytropic
   * | eos.Press_PPEOS[4] |      EOS case
   * .--------------------.
   * |  Return value: 0   |
   * .--------------------.
   *
   * .--------------------.
   * |  Return value: 1   | Single Polytropic EOS case
   * .--------------------.
   *
   * .--------------------.
   * |  Return value: 2   | Unsupported/unknown EOS
   * .--------------------.
   */
  int EOS_type = get_EOS_parameters_from_EOSname(EOSname,eos);

  /* Check if no error occured */
  if( EOS_type == 2 ) {
    printf("ERROR! Unrecognized EOS name: %s\n",EOSname.c_str());
    exit(2);
  }

  /* Single polytrope case */
  if( EOS_type == 1 ) {
    /* .---------------------------------.
     * | Single Polytrope Initialization |
     * .---------------------------------.
     * For the case of a single polytrope,
     * we initialize the EOS parameters:
     * .--------------------------.
     * | eos.neos           =  1  |
     * | eos.Kpoly_PPEOS[0] = 1.0 |
     * | eos.Gamma_PPEOS[0] = 2.0 |
     * .--------------------------.
     */
    eos.neos           =  1 ;
    eos.Kpoly_PPEOS[0] = 1.0;
    eos.Gamma_PPEOS[0] = 2.0;

    /* We are done with the initialization, no need to continue */
    return;
  }

  /* Read et al. tables II and III case */
  if( EOS_type == 0 ) {
    /* .------------------------------------.
     * | Piecewise Polytrope Initialization |
     * |     Read et al. EOS parameters     |
     * .------------------------------------.
     * For the case of a piecewise polytrope, we initialize the following EOS parameters
     * for *all* EOSnames. These values can be found in table II and Fig 3. of Read et al.
     * (see references at the top of this file)
     * .----------------------------------.
     * | eos.neos = 7                     |
     * .----------------------------------.
     * | eos.rho_b_PPEOS[0] = 2.44034e+07 |
     * | eos.rho_b_PPEOS[1] = 3.78358e+11 |
     * | eos.rho_b_PPEOS[2] = 2.62780e+12 |
     * | eos.rho_b_PPEOS[4] = 5.01187e+14 |
     * | eos.rho_b_PPEOS[5] = 1.00000e+15 |
     * .----------------------------------.
     * | eos.Kpoly_PPEOS[0] = 6.80110e-09 |
     * .----------------------------------.
     * | eos.Gamma_PPEOS[0] = 1.58425     |
     * | eos.Gamma_PPEOS[1] = 1.28733     |
     * | eos.Gamma_PPEOS[2] = 0.62223     |
     * | eos.Gamma_PPEOS[3] = 1.35692     |
     * .----------------------------------.
     */
    eos.neos = 7;
    eos.rho_b_PPEOS[0] = 2.44034e+07;
    eos.rho_b_PPEOS[1] = 3.78358e+11;
    eos.rho_b_PPEOS[2] = 2.62780e+12;
    eos.rho_b_PPEOS[4] = 5.01187e+14;
    eos.rho_b_PPEOS[5] = 1.00000e+15;
    eos.Kpoly_PPEOS[0] = 6.80110e-09;
    eos.Gamma_PPEOS[0] = 1.58425;
    eos.Gamma_PPEOS[1] = 1.28733;
    eos.Gamma_PPEOS[2] = 0.62223;
    eos.Gamma_PPEOS[3] = 1.35692;
  }

  /* At this point we still need to determine:
   *  - rho_{3} = eos.rho_b_PPEOS[3]
   *  - K_{j}   = eos.Kpoly_PPEOS[j], 1 <= j < eos.neos
   *  - P_{j}   = eos.Press_PPEOS[j], j != 4
   *  - C_{j}   = eos.epsIC_PPEOS[j], 0 <= j < eos.neos
   *
   * We start by computing K_{j} and P_{j} using the
   * values of rho_{j} we have available at this point
   * (see functions populate_K_PPEOS and populate_PRESS_PPEOS)
   */
  for(int j=1; j <= 3; j++) {
    eos.Kpoly_PPEOS[j] = eos.Kpoly_PPEOS[j-1]*pow(eos.rho_b_PPEOS[j-1],eos.Gamma_PPEOS[j-1]-eos.Gamma_PPEOS[j]);
  }
  for(int j=0; j <= 2; j++) {
    eos.Press_PPEOS[j] = eos.Kpoly_PPEOS[j] * pow(eos.rho_b_PPEOS[j],eos.Gamma_PPEOS[j]);
  }
  
  /* Computing rho_{3} is then done in a slightly unintuitive way.
   * Remember that we currently know the values of P_{4}, rho_{4},
   * and Gamma_{4}. Therefore, we compute
   * .-------------------------------------.
   * | K_{4} = P_{4} / rho_{4}^(Gamma_{4}) |
   * .-------------------------------------.
   */
  eos.Kpoly_PPEOS[4] = eos.Press_PPEOS[4] / pow(eos.rho_b_PPEOS[4],eos.Gamma_PPEOS[4]);

  /* Then, continuity of the pressure tells us that
   *
   * K_{4} = K_{3} rho_{3}^( Gamma_{3} - Gamma_{4} ) ,
   *
   * which implies
   * .-----------------------------------------------------.
   * | rho_{3) = ( K_{4}/K_{3} )^(1/(Gamma_{3}-Gamma_{4})) |
   * .-----------------------------------------------------.
   */
  eos.rho_b_PPEOS[3] = pow( eos.Kpoly_PPEOS[4]/eos.Kpoly_PPEOS[3], 1.0/(eos.Gamma_PPEOS[3] - eos.Gamma_PPEOS[4]) );
  
  /* Then, we finish populating K_{j} */
  for(int j=5; j < eos.neos; j++) {
    eos.Kpoly_PPEOS[j] = eos.Kpoly_PPEOS[j-1]*pow(eos.rho_b_PPEOS[j-1],eos.Gamma_PPEOS[j-1]-eos.Gamma_PPEOS[j]);
  }
  
  /* And finally we finish populating P_{j} */
  for(int j=3; j < eos.neos-1; j++) {
    eos.Press_PPEOS[j] = eos.Kpoly_PPEOS[j] * pow(eos.rho_b_PPEOS[j],eos.Gamma_PPEOS[j]);
  }
  
  /* At this point we have populated:
   *
   *  - eos.neos
   *  - eos.rho_b_PPEOS
   *  - eos.Kpoly_PPEOS
   *  - eos.Gamma_PPEOS
   *  - eos.Press_PPEOS
   *
   * and would be left with the task of of populating eos.epsIC_PPEOS. However,
   * the values we have populated the struct with thus far are in units such that
   * [P] = [rho], where [A] indicate the units of the quantity A. This is not a desireable
   * feature, since we end up having to work with very small/large numbers throughout the TOV
   * solver. Instead, what we will do it convert the results we have this far to *geometrized*
   * units, where G = 1 = c. In doing so, we will  also populate the eos.epsIC_PPEOS array.
   * For more details, please look below for the implementation of the function
   * convert_EOS_struct_to_geometrized_units().
   */
  convert_EOS_struct_to_geometrized_units(eos);
  
}

/* Function   : convert_EOS_struct_to_geometrized_units()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Compute Kpoly_PPEOS by imposing that the pressure be
 *              everywhere continuous.
 *
 * Input(s)   : eos - Struct containing EOS information in [P] = [rho] units
 *
 * Outputs(s) : eos - Fully populated struct in G = 1 = c units
 */
inline void convert_EOS_struct_to_geometrized_units( eos_struct &eos ) {

  /* The values of the polytropic EOS parameters we have used are in units
   * such that [P] = [rho] = g/cm^{3}. In order to keep quantities of order
   * unity and not have to worry about keeping track of powers of c (speed
   * of light) and G (gravitational constant), we will rewrite all quantities
   * in *geometrized* units, where G = 1 = c. In these units,
   * .------------------------.
   * | [P] = [rho] = 1/cm^{2} |
   * .------------------------.
   * Now we want to rewrite our quantities wuch that [rho_{neos-1}] = 1. We
   * also want this rescaling to preserve the ratios rho_{j}/rho_{j-1} we
   * had before. Thus, we have
   * .-------------------------------------------------------------.
   * | rho_rescaled_{0} = 1                                        |
   * | rho_rescaled_{j} = rho_rescaled_{j-1} * (rho_{j}/rho_{j-1}) |
   * .-------------------------------------------------------------.
   */
  REAL rho_b_rescaled[eos.neos-1];
  rho_b_rescaled[eos.neos-2] = 1.0;
  for(int j=eos.neos-2; j>0; j--) {
    rho_b_rescaled[j-1] = rho_b_rescaled[j]*(eos.rho_b_PPEOS[j-1]/eos.rho_b_PPEOS[j]);
  }
  /* Because [P] = [rho], they must be rescaled by the same factor, namely
   *
   * P_rescaled_{j}/P_{j} = rho_rescaled_{j}/rho_{j} ,
   *
   * implying the relation
   * .-----------------------------------------------------.
   * | P_rescaled_{j} = P_{j} * (rho_rescaled_{j}/rho_{j}) |
   * .-----------------------------------------------------. 
   */
  REAL Press_rescaled[eos.neos-1];
  for(int j=0; j<eos.neos-1; j++) {
    Press_rescaled[j] = eos.Press_PPEOS[j] * (rho_b_rescaled[j]/eos.rho_b_PPEOS[j]);
  }

  /* Then, update our struct with the rescaled values of P and rho_b */
  for(int j=0; j<eos.neos-1; j++) {
    eos.rho_b_PPEOS[j] = rho_b_rescaled[j];
    eos.Press_PPEOS[j] = Press_rescaled[j];
  }

  /* Recompute K_{0} and then update K_{j}, 1 <= j < eos.neos */
  eos.Kpoly_PPEOS[0] = eos.Press_PPEOS[0]/pow(eos.rho_b_PPEOS[0],eos.Gamma_PPEOS[0]);
  populate_Kpoly_PPEOS(eos);

  /* Compute the integration constants of eps */
  populate_epsIC_PPEOS(eos);

}

/* Function   : populate_Kpoly_PPEOS()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Compute Kpoly_PPEOS by imposing that the pressure be
 *              everywhere continuous.
 *
 * Input(s)   : eos             - Struct containing EOS information
 *
 * Outputs(s) : eos.Kpoly_PPEOS - Fully populated array
 */
inline void populate_Kpoly_PPEOS( eos_struct &eos ) {
  
  /* Remember that for a piecewise polytropic EOS, we have
   * .---------------------------------------------------------------------------------.
   * |             / K_{0}   rho_b ^ (Gamma_{0})  , if              rho_b <= rho_{0}   |
   * |             | K_{1}   rho_b ^ (Gamma_{1})  , if   rho_{0} <= rho_b <= rho_{1}   |
   * |             |             ...              ,                  ...               |
   * | P(rho_b) = <  K_{j}   rho_b ^ (Gamma_{j})  , if rho_{j-1} <= rho_b <= rho_{j}   |
   * |             |             ...              ,                  ...               |
   * |             | K_{N-2} rho_b ^ (Gamma_{N-2}), if rho_{N-3} <= rho_b <= rho_{N-2} |
   * |             \ K_{N-1} rho_b ^ (Gamma_{N-1}), if              rho_b >= rho_{N-2} |
   * .---------------------------------------------------------------------------------.
   * where N = eos.neos. Thus, imposing continuity for the first two EOSs, we would have
   *
   * P(rho_{0}) = K_{0} rho_{0} ^ (Gamma_{0}) = K_{1} rho_{0} ^ (Gamma_{1})
   *
   *       => K_{1} = K_{0} rho_{0} ^ (Gamma_{0} - Gamma_{1})
   *
   * Doing this for a generic case we end up with
   * .---------------------------------------------------------.
   * | K_{j} = K_{j-1} * rho_{j-1} ^ (Gamma_{j-1} - Gamma_{j}) |
   * .---------------------------------------------------------.
   */
  for(int j=1; j<eos.neos; j++)
    eos.Kpoly_PPEOS[j] = eos.Kpoly_PPEOS[j-1] * pow(eos.rho_b_PPEOS[j-1],eos.Gamma_PPEOS[j-1] - eos.Gamma_PPEOS[j]);
}

/* Function   : populate_Press_PPEOS()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Compute Press_PPEOS using the appropriate EOSs.
 *
 * Input(s)   : eos             - Struct containing EOS information
 *
 * Outputs(s) : eos.Press_PPEOS - Fully populated array
 */
inline void populate_Press_PPEOS( eos_struct &eos ) {
  /* In this function we simply compute the pressure values
   * associated with rho_{j}, i.e.
   * .---------------------------------------.
   * | P_{j} = K_{j} * rho_{j} ^ (Gamma_{j}) |
   * .---------------------------------------.
   */
  for(int j=0; j<eos.neos-1; j++)
    eos.Press_PPEOS[j] = eos.Kpoly_PPEOS[j] * pow(eos.rho_b_PPEOS[j],eos.Gamma_PPEOS[j]);
}

/* Function   : populate_epsIC_PPEOS()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Compute Kpoly_PPEOS by imposing that the pressure be
 *              everywhere continuous.
 *
 * Input(s)   : eos             - Struct containing EOS information
 *
 * Outputs(s) : eos.Kpoly_PPEOS - Fully populated array
 */
inline void populate_epsIC_PPEOS( eos_struct &eos ) {
  
  /* Remember that for a piecewise polytropic EOS, we have
   * .-------------------------------------------------------------------------------------------------.
   * |        / C_{0}   + K_{0}  rho_b^( Gamma_{0}-1 )/( Gamma_{0}-1 ), if            rho_b<=rho_{0}   |
   * |        | C_{1}   + K_{1}  rho_b^( Gamma_{1}-1 )/( Gamma_{1}-1 ), if   rho_{0}<=rho_b<=rho_{1}   |
   * |        |                           ...                         ,                ...             |
   * | eps = <  C_{j}   + K_{j}  rho_b^( Gamma_{j}-1 )/( Gamma_{j}-1 ), if rho_{j-1}<=rho_b<=rho_{j}   |
   * |        |                           ...                         ,                ...             |
   * |        | C_{N-2} + K_{N-2}rho_b^(Gamma_{N-2}-1)/(Gamma_{N-2}-1), if rho_{N-3}<=rho_b<=rho_{N-2} |
   * |        \ C_{N-1} + K_{N-1}rho_b^(Gamma_{N-1}-1)/(Gamma_{N-1}-1), if            rho_b>=rho_{N-2} |
   * .-------------------------------------------------------------------------------------------------.
   * where the C's are the integration constants which we refer to as epsIC_PPEOS,
   * and N = eos.neos. Our first observation is the demand that
   * .-------------------------------.
   * | eps(rho_b=0) = 0 => C_{0} = 0 |
   * .-------------------------------.
   * Then, imposing that eps be everywhere continuous, we have
   *
   * eps(rho_{0}) = C_{0} + K_{0} rho_{0} ^ (Gamma_{0}-1)/(Gamma_{0}-1)
   *              = C_{1} + K_{1} rho_{0} ^ (Gamma_{1}-1)/(Gamma_{1}-1)
   *
   * => C_{1} = C_{0} 
   *          + K_{0} rho_{0} ^ (Gamma_{0}-1)/(Gamma_{0}-1)
   *          - K_{1} rho_{0} ^ (Gamma_{1}-1)/(Gamma_{1}-1)
   *
   * Doing this for a generic case we end up with
   * .-------------------------------------------------------------.
   * | C_{j} = C_{j-1}                                             |
   * |       + K_{j-1} rho_{j-1} ^ (Gamma_{j-1}-1)/(Gamma_{j-1}-1) |
   * |       - K_{ j } rho_{j-1} ^ (Gamma_{ j }-1)/(Gamma_{ j }-1) |
   * .-------------------------------------------------------------.
   */
  eos.epsIC_PPEOS[0] = 0.0;
  for(int j=1; j<eos.neos; j++) {
    /* Auxiliary variable: computes the second term in
     *                     the boxed equation above
     */
    REAL aux_jm1 = eos.Kpoly_PPEOS[j-1]*pow(eos.rho_b_PPEOS[j-1],eos.Gamma_PPEOS[j-1]-1.0)/(eos.Gamma_PPEOS[j-1]-1.0);
    /* Auxiliary variable: computes the third term in
     *                     the boxed equation above
     */
    REAL aux_j = eos.Kpoly_PPEOS[j]*pow(eos.rho_b_PPEOS[j-1],eos.Gamma_PPEOS[j]-1.0)/(eos.Gamma_PPEOS[j]-1.0);
    /* Update the integration constant according to the boxed equation above */
    eos.epsIC_PPEOS[j] = eos.epsIC_PPEOS[j-1] + aux_jm1 - aux_j;
  }

}

/* Function   : get_EOS_parameters_from_EOSname()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Compute Kpoly_PPEOS by imposing that the pressure be
 *              everywhere continuous.
 *
 * Input(s)   : EOSname            - Name of the EOS. Can be "Single" or any of the
 *                                   EOS names in table III of Read et al.
 *            : eos                - Struct containing EOS information
 *
 * Outputs(s) : return value:
 *                               0 - Piecewise Polytropic EOS case
 *                               1 - Single Polytropic EOS case
 *                               2 - Error case, unknown EOS
 *
 *            : eos.Press_PPEOS[4] - Value from table III of Read et al. (Piecewise Polytropic EOS only)
 *              eos.Gamma_PPEOS[4] - Value from table III of Read et al. (Piecewise Polytropic EOS only)
 *              eos.Gamma_PPEOS[5] - Value from table III of Read et al. (Piecewise Polytropic EOS only)
 *              eos.Gamma_PPEOS[6] - Value from table III of Read et al. (Piecewise Polytropic EOS only)

 */
inline int get_EOS_parameters_from_EOSname( string EOSname, eos_struct &eos ) {

  if( EOSname == "Single" )
  {
    return 1;
  }
  else if( EOSname == "PAL6" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.380 2.227 2.189 2.159
     */
    eos.Press_PPEOS[4] = 34.380; eos.Gamma_PPEOS[4] = 2.227; eos.Gamma_PPEOS[5] = 2.189; eos.Gamma_PPEOS[6] = 2.159;
  }
  else if( EOSname == "SLy" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.384 3.005 2.988 2.851
     */
    eos.Press_PPEOS[4] = 34.384; eos.Gamma_PPEOS[4] = 3.005; eos.Gamma_PPEOS[5] = 2.988; eos.Gamma_PPEOS[6] = 2.851;
  }
    
  else if( EOSname == "APR1" )
  {
    /* Copy & paste from Table III of Read et al.
     * 33.943 2.442 3.256 2.908
     */
    eos.Press_PPEOS[4] = 33.943; eos.Gamma_PPEOS[4] = 2.442; eos.Gamma_PPEOS[5] = 3.256; eos.Gamma_PPEOS[6] = 2.908;
  }
    
  else if( EOSname == "APR2" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.126 2.643 3.014 2.945
     */
    eos.Press_PPEOS[4] = 34.126; eos.Gamma_PPEOS[4] = 2.643; eos.Gamma_PPEOS[5] = 3.014; eos.Gamma_PPEOS[6] = 2.945;
  }
    
  else if( EOSname == "APR3" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.392 3.166 3.573 3.281
     */
    eos.Press_PPEOS[4] = 34.392; eos.Gamma_PPEOS[4] = 3.166; eos.Gamma_PPEOS[5] = 3.573; eos.Gamma_PPEOS[6] = 3.281;
  }
    
  else if( EOSname == "APR4" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.269 2.830 3.445 3.348
     */
    eos.Press_PPEOS[4] = 34.269; eos.Gamma_PPEOS[4] = 2.830; eos.Gamma_PPEOS[5] = 3.445; eos.Gamma_PPEOS[6] = 3.348;
  }
    
  else if( EOSname == "FPS" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.283 2.985 2.863 2.600
     */
    eos.Press_PPEOS[4] = 34.283; eos.Gamma_PPEOS[4] = 2.985; eos.Gamma_PPEOS[5] = 2.863; eos.Gamma_PPEOS[6] = 2.600;
  }
    
  else if( EOSname == "WFF1" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.031 2.519 3.791 3.660
     */
    eos.Press_PPEOS[4] = 34.031; eos.Gamma_PPEOS[4] = 2.519; eos.Gamma_PPEOS[5] = 3.791; eos.Gamma_PPEOS[6] = 3.660;
  }
    
  else if( EOSname == "WFF2" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.233 2.888 3.475 3.517
     */
    eos.Press_PPEOS[4] = 34.233; eos.Gamma_PPEOS[4] = 2.888; eos.Gamma_PPEOS[5] = 3.475; eos.Gamma_PPEOS[6] = 3.517;
  }
    
  else if( EOSname == "WFF3" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.283 3.329 2.952 2.589
     */
    eos.Press_PPEOS[4] = 34.283; eos.Gamma_PPEOS[4] = 3.329; eos.Gamma_PPEOS[5] = 2.952; eos.Gamma_PPEOS[6] = 2.589;
  }
    
  else if( EOSname == "BBB2" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.331 3.418 2.835 2.832
     */
    eos.Press_PPEOS[4] = 34.331; eos.Gamma_PPEOS[4] = 3.418; eos.Gamma_PPEOS[5] = 2.835; eos.Gamma_PPEOS[6] = 2.832;
  }
    
  else if( EOSname == "BPAL12" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.358 2.209 2.201 2.176
     */
    eos.Press_PPEOS[4] = 34.358; eos.Gamma_PPEOS[4] = 2.209; eos.Gamma_PPEOS[5] = 2.201; eos.Gamma_PPEOS[6] = 2.176;
  }
    
  else if( EOSname == "ENG" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.437 3.514 3.130 3.168
     */
    eos.Press_PPEOS[4] = 34.437; eos.Gamma_PPEOS[4] = 3.514; eos.Gamma_PPEOS[5] = 3.130; eos.Gamma_PPEOS[6] = 3.168;
  }
    
  else if( EOSname == "MPA1" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.495 3.446 3.572 2.887
     */
    eos.Press_PPEOS[4] = 34.495; eos.Gamma_PPEOS[4] = 3.446; eos.Gamma_PPEOS[5] = 3.572; eos.Gamma_PPEOS[6] = 2.887;
  }
    
  else if( EOSname == "MS1" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.858 3.224 3.033 1.325
     */
    eos.Press_PPEOS[4] = 34.858; eos.Gamma_PPEOS[4] = 3.224; eos.Gamma_PPEOS[5] = 3.033; eos.Gamma_PPEOS[6] = 1.325;
  }
    
  else if( EOSname == "MS2" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.605 2.447 2.184 1.855
     */
    eos.Press_PPEOS[4] = 34.605; eos.Gamma_PPEOS[4] = 2.447; eos.Gamma_PPEOS[5] = 2.184; eos.Gamma_PPEOS[6] = 1.855;
  }
    
  else if( EOSname == "MS1b" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.855 3.456 3.011 1.425
     */
    eos.Press_PPEOS[4] = 34.855; eos.Gamma_PPEOS[4] = 3.456; eos.Gamma_PPEOS[5] = 3.011; eos.Gamma_PPEOS[6] = 1.425;
  }
    
  else if( EOSname == "PS" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.671 2.216 1.640 2.365
     */
    eos.Press_PPEOS[4] = 34.671; eos.Gamma_PPEOS[4] = 2.216; eos.Gamma_PPEOS[5] = 1.640; eos.Gamma_PPEOS[6] = 2.365;
  }
    
  else if( EOSname == "GS1" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.504 2.350 1.267 2.421
     */
    eos.Press_PPEOS[4] = 34.504; eos.Gamma_PPEOS[4] = 2.350; eos.Gamma_PPEOS[5] = 1.267; eos.Gamma_PPEOS[6] = 2.421;
  }
    
  else if( EOSname == "GS2" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.642 2.519 1.571 2.314
     */
    eos.Press_PPEOS[4] = 34.642; eos.Gamma_PPEOS[4] = 2.519; eos.Gamma_PPEOS[5] = 1.571; eos.Gamma_PPEOS[6] = 2.314;
  }

  else if( EOSname == "BGN1H1" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.623 3.258 1.472 2.464
     */
    eos.Press_PPEOS[4] = 34.623; eos.Gamma_PPEOS[4] = 3.258; eos.Gamma_PPEOS[5] = 1.472; eos.Gamma_PPEOS[6] = 2.464;
  }
    
  else if( EOSname == "GNH3" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.648 2.664 2.194 2.304
     */
    eos.Press_PPEOS[4] = 34.648; eos.Gamma_PPEOS[4] = 2.664; eos.Gamma_PPEOS[5] = 2.194; eos.Gamma_PPEOS[6] = 2.304;
  }
    
  else if( EOSname == "H1" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.564 2.595 1.845 1.897
     */
    eos.Press_PPEOS[4] = 34.564; eos.Gamma_PPEOS[4] = 2.595; eos.Gamma_PPEOS[5] = 1.845; eos.Gamma_PPEOS[6] = 1.897;
  }
    
  else if( EOSname == "H2" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.617 2.775 1.855 1.858
     */
    eos.Press_PPEOS[4] = 34.617; eos.Gamma_PPEOS[4] = 2.775; eos.Gamma_PPEOS[5] = 1.855; eos.Gamma_PPEOS[6] = 1.858;
  }
    
  else if( EOSname == "H3" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.646 2.787 1.951 1.901
     */
    eos.Press_PPEOS[4] = 34.646; eos.Gamma_PPEOS[4] = 2.787; eos.Gamma_PPEOS[5] = 1.951; eos.Gamma_PPEOS[6] = 1.901;
  }
    
  else if( EOSname == "H4" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.669 2.909 2.246 2.144
     */
    eos.Press_PPEOS[4] = 34.669; eos.Gamma_PPEOS[4] = 2.909; eos.Gamma_PPEOS[5] = 2.246; eos.Gamma_PPEOS[6] = 2.144;
  }
    
  else if( EOSname == "H5" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.609 2.793 1.974 1.915
     */
    eos.Press_PPEOS[4] = 34.609; eos.Gamma_PPEOS[4] = 2.793; eos.Gamma_PPEOS[5] = 1.974; eos.Gamma_PPEOS[6] = 1.915;
  }
    
  else if( EOSname == "H6" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.593 2.637 2.121 2.064
     */
    eos.Press_PPEOS[4] = 34.593; eos.Gamma_PPEOS[4] = 2.637; eos.Gamma_PPEOS[5] = 2.121; eos.Gamma_PPEOS[6] = 2.064;
  }
    
  else if( EOSname == "H7" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.559 2.621 2.048 2.006
     */
    eos.Press_PPEOS[4] = 34.559; eos.Gamma_PPEOS[4] = 2.621; eos.Gamma_PPEOS[5] = 2.048; eos.Gamma_PPEOS[6] = 2.006;
  }
    
  else if( EOSname == "PCL2" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.507 2.554 1.880 1.977
     */
    eos.Press_PPEOS[4] = 34.507; eos.Gamma_PPEOS[4] = 2.554; eos.Gamma_PPEOS[5] = 1.880; eos.Gamma_PPEOS[6] = 1.977;
  }
  
  else if( EOSname == "ALF1" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.055 2.013 3.389 2.033
     */
    eos.Press_PPEOS[4] = 34.055; eos.Gamma_PPEOS[4] = 2.013; eos.Gamma_PPEOS[5] = 3.389; eos.Gamma_PPEOS[6] = 2.033;
  }
  
  else if( EOSname == "ALF2" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.616 4.070 2.411 1.890
     */
    eos.Press_PPEOS[4] = 34.616; eos.Gamma_PPEOS[4] = 4.070; eos.Gamma_PPEOS[5] = 2.411; eos.Gamma_PPEOS[6] = 1.890;
  }
  
  else if( EOSname == "ALF3" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.283 2.883 2.653 1.952
     */
    eos.Press_PPEOS[4] = 34.283; eos.Gamma_PPEOS[4] = 2.883; eos.Gamma_PPEOS[5] = 2.653; eos.Gamma_PPEOS[6] = 1.952;
  }
  
  else if( EOSname == "ALF4" )
  {
    /* Copy & paste from Table III of Read et al.
     * 34.314 3.009 3.438 1.803
     */
    eos.Press_PPEOS[4] = 34.314; eos.Gamma_PPEOS[4] = 3.009; eos.Gamma_PPEOS[5] = 3.438; eos.Gamma_PPEOS[6] = 1.803;
  }
  else
  {
    return 2;
  }

  /* The input of log(p1) from table III of Read et al. is in CGS units */
  REAL log10_pressure_in_cgs  = eos.Press_PPEOS[4];

  /* Convert this to units where [P] = [rho] by dividing P by c^(2) */
  /* Note: The speed of light is a global variable, defined in TOV_headers.h */
  REAL log10_pressure_over_c2 = log10_pressure_in_cgs - 2.0*log10(SPEEDOFLIGHT);

  /* Update the eos value of P_(4) */
  eos.Press_PPEOS[4] = pow(10.0,log10_pressure_over_c2);

  /* Return 0, which indicates the Piecewise Polytropic EOS case */
  return 0;

}
