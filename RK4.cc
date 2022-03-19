/* .-----------------------------------------.
 * |  Copyright (c) 2019, Leonardo Werneck   |
 * | Licensed under the BSD 2-Clause License |
 * .-----------------------------------------.
 */

/* Program     : TOV Solver
 * File        : RK4.C
 * Author      : Leo Werneck (werneck@if.usp.br)
 * Date        : October 29, 2019
 *
 * Description : This file implements a fourth-order Runge-Kutta
 *               integrator that is used to solve the TOV equations.
 *
 * Dependencies: tov_headers.h & TOV_RHSs.C
 *
 * Reference(s): https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
 */

/* Include dependencies */
#include "tov_headers.h"

/* Function   : RK4_method()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Performs a single RK4 integration step
 *
 * Input(s)   : rr                - Current value of the radius, the integration variable
 *            : dr                - Current value of the step size
 *            : gfs_at_rr         - Array containing the gridfunctions at r = rr
 *            : gfs_at_rr_plus_dr - Array to store the gridfunctions at r = rr + dr
 *
 * Outputs(s) : gfs_at_rr_plus_dr - Fully populated array
 */
void RK4_method( eos_struct eos, const REAL rr, const REAL dr,  const REAL *restrict gfs_at_rr, REAL *restrict gfs_at_rr_plus_dr ) {

  /******************
   * The RK4 Method *
   ******************
   *
   * The RK4 method consists of solving a differential equation,
   * or systems thereof, of the form
   * .----------------.
   * | y'(r) = f(r,y) |
   * .----------------.
   * from a given initial condition y(0). The algorithm then
   * computes:
   * .------------------------------------------------.
   * | k1 = dr*f(r,y)              <- 1st RK4 step    |
   * .------------------------------------------------.
   * | k2 = dr*f(r+dr/2,y+k1/2)    <- 2nd RK4 step    |
   * .------------------------------------------------.
   * | k3 = dr*f(r+dr/2,y+k2/2)    <- 3rd RK4 step    |
   * .------------------------------------------------.
   * | k4 = dr*f(r+dr,y+k3)        <- 4th RK4 step    |
   * .------------------------------------------------.
   * | y(r+dr) = y(r) + (1/6)*(k1 + 2*k2 + 2*k3 + k4) |
   * .------------------------------------------------.
   *
   * We start by declarating variables to store the RHSs
   * of the differential equations and {k1,k2,k3,k4}.
   */
  REAL RHS[NGFS], k1[NGFS], k2[NGFS], k3[NGFS], k4[NGFS];

  /****************
   * 1st RK4 STEP *
   ****************
   * Compute the RHS of dP/dr, dnu/dr, and dm/dr,
   * and store them in RHS.
   * Remember: gfs_at_rr = { P(r), nu(r), m(r) }.
   */
  TOV_RHSs( eos, rr, gfs_at_rr, RHS);
  /* Then compute k1 */
  for(int GF=0; GF<NGFS; GF++) { k1[GF] = dr * RHS[GF]; }

  /****************
   * 2nd RK4 STEP *
   ****************
   * Set up gfs_at_rr_plus_dr as the auxiliary variable:
   * gfs_at_rr_plus_dr = { P(r) + k1_P/2 , nu(r) + k1_nu/2, m(r) + k1_m/2 }.
   */
  for(int GF=0; GF<NGFS; GF++) { gfs_at_rr_plus_dr[GF] = gfs_at_rr[GF] + k1[GF]/2.0; }
  /*
   * Compute the RHS of dP/dr, dnu/dr, and dm/dr
   * Remember: gfs_at_rr_plus_dr = gfs_at_rr + k1/2.
   */
  TOV_RHSs( eos, rr+dr/2.0, gfs_at_rr_plus_dr, RHS );
  /* Then compute k2 */
  for(int GF=0; GF<NGFS; GF++) { k2[GF]  = dr * RHS[GF]; }

  /****************
   * 3rd RK4 STEP *
   ****************
   * Set up gfs_at_rr_plus_dr as the auxiliary variable:
   * gfs_at_rr_plus_dr = { P(r) + k2_P/2 , nu(r) + k2_nu/2, m(r) + k2_m/2 }.
   */
  for(int GF=0; GF< NGFS; GF++) { gfs_at_rr_plus_dr[GF] = gfs_at_rr[GF] + k2[GF]/2.0; }
  /*
   * Compute the RHS of dP/dr, dnu/dr, and dm/dr
   * Remember: gfs_at_rr_plus_dr = gfs_at_rr + k2/2.
   */
  TOV_RHSs( eos, rr+dr/2.0, gfs_at_rr_plus_dr, RHS);
  /* Then compute k3 */
  for(int GF=0; GF<NGFS; GF++) { k3[GF]  = dr * RHS[GF]; }

  /****************
   * 4th RK4 STEP *
   ****************
   * Set up gfs_at_rr_plus_dr as the auxiliary variable:
   * gfs_at_rr_plus_dr = { P(r) + k3_P , nu(r) + k3_nu, m(r) + k3_m }.
   */
  for(int GF=0; GF<NGFS; GF++) { gfs_at_rr_plus_dr[GF] = gfs_at_rr[GF] + k3[GF]; }
  /*
   * Compute the RHS of dP/dr, dnu/dr, and dm/dr
   * Remember: gfs_at_rr_plus_dr = gfs_at_rr + k3.
   */
  TOV_RHSs( eos, rr+dr, gfs_at_rr_plus_dr, RHS);
  /* Then compute k4 */
  for(int GF=0; GF<NGFS; GF++) { k4[GF]  = dr * RHS[GF]; }

  /* Finally, Update gfs_at_rr_plus_dr with the values of
   * the gridfunctions at the radius rr + dr
   */
  for(int GF=0; GF<NGFS; GF++) { gfs_at_rr_plus_dr[GF] = gfs_at_rr[GF] + (1.0/6.0) * ( k1[GF] + 2.0*k2[GF] + 2.0*k3[GF] + k4[GF] ); }

  if(DEBUG) printf("(DEBUGGING INFO - RK4_method()) P(r) = %.15e | nu(r) = %.15e | m(r) = %.15e \n",gfs_at_rr[PRESSURE],gfs_at_rr[NU],gfs_at_rr[MASS]);
  if(DEBUG) printf("(DEBUGGING INFO - RK4_method()) k1_P = %.15e | k1_nu = %.15e | k1_m = %.15e \n",k1[PRESSURE],k1[NU],k1[MASS]);
  if(DEBUG) printf("(DEBUGGING INFO - RK4_method()) k2_P = %.15e | k2_nu = %.15e | k2_m = %.15e \n",k2[PRESSURE],k2[NU],k2[MASS]);
  if(DEBUG) printf("(DEBUGGING INFO - RK4_method()) k3_P = %.15e | k3_nu = %.15e | k3_m = %.15e \n",k3[PRESSURE],k3[NU],k3[MASS]);
  if(DEBUG) printf("(DEBUGGING INFO - RK4_method()) k4_P = %.15e | k4_nu = %.15e | k4_m = %.15e \n",k4[PRESSURE],k4[NU],k4[MASS]);
  if(DEBUG) printf("(DEBUGGING INFO - RK4_method()) P(r+dr) = %.15e | nu(r+dr) = %.15e | m(r+dr) = %.15e \n",gfs_at_rr_plus_dr[PRESSURE],gfs_at_rr_plus_dr[NU],gfs_at_rr_plus_dr[MASS]);

}

/* Function   : recompute_step_size()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Performs a single RK4 integration step
 *
 * Input(s)   : gfs_at_rr         - Array containing the gridfunctions at r = rr
 *            : gfs_at_rr_plus_dr - Array to store the gridfunctions at r = rr + dr
 *            : dr                - Current value of the step size
 *
 * Outputs(s) : dr                - Modified, if necessary, for adaptive integration
 */
void recompute_step_size( const eos_struct eos, const REAL rr, const REAL *gfs_at_rr, REAL &dr ) {

  /* We will now write the algorithm to change the step size of
   * our integration, when necessary, making our method adaptive.
   * We will choose to to so based on the ratio of our current
   * integration variables and their derivatives. We will, also,
   * set up a maximum step size, to avoid our solution to be too
   * poorly sampled.
   *
   * We start by computing the derivatives ( dP/dr, dm/dr, dnu/dr )
   * for our current value of r.
   */
  REAL RHS[NGFS];
  TOV_RHSs(eos, rr, gfs_at_rr, RHS);

  /* Compute our new step size */
  REAL new_step_size = 0.1*MIN( fabs(gfs_at_rr[PRESSURE]/RHS[PRESSURE]), fabs(gfs_at_rr[MASS]/RHS[MASS]) );

  /* Make sure our new step size does not exceed the maximum allowed step size */
  new_step_size = MIN( new_step_size, MAXIMUMSTEPSIZE );

  /* Update our step size */
  dr = new_step_size;

}
