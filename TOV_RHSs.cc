/* .-----------------------------------------.
 * |  Copyright (c) 2019, Leonardo Werneck   |
 * | Licensed under the BSD 2-Clause License |
 * .-----------------------------------------.
 */

/* Program     : TOV Solver
 * File        : TOV_RHSs.C
 * Author      : Leo Werneck (werneck@if.usp.br)
 * Date        : October 29, 2019
 *
 * Description : This file implements the RHSs of the Tolman–Oppenheimer–Volkoff
 *               equations
 *
 * Dependencies: math.h, tov_headers.h & Polytropic_EOS__lowlevel_functions.C
 *
 * Reference(s): https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation
 */

/* Include dependencies */
#include "tov_headers.h"

/* Function   : TOV_RHSs()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Evaluates the RHSs of the TOV equations
 *
 * Input(s)   : rr        - Current value of the radius, the integration variable
 *            : gfs_at_rr - Array containing the gridfunctions at r = rr
 *            : RHS       - Array that will hold the RHSs of the TOV equations at r = rr
 *
 * Outputs(s) : RHS       - Fully populated array with the RHSs of the TOV equations at r = rr
 */
void TOV_RHSs( const eos_struct eos, const REAL rr, const REAL *restrict gfs_at_rr, REAL *restrict RHS ) {

  /**************************************************
   * The Tolman–Oppenheimer–Volkoff (TOV) Equations *
   **************************************************
   *
   * The TOV equations are:
   * .-----------------------------------------------------.
   * | dP/dr  = -(mu+P)(m + 4 pi r^3 P)/( r^2 (1 - 2m/r) ) |
   * .-----------------------------------------------------.
   * | dnu/dr = -(2 / (mu+P) ) dP/dr                       |
   * .-----------------------------------------------------.
   * | dm/dr  = 4 pi r^2 mu                                |
   * .-----------------------------------------------------.
   * We start by definiing useful auxiliary variables:
   */
  REAL P  = gfs_at_rr[PRESSURE];
  REAL m  = gfs_at_rr[MASS];
  REAL mu = compute_mu_from_P(eos,gfs_at_rr[PRESSURE]);
  REAL r2 = rr*rr;
  REAL r3 = rr*r2;

  if(DEBUG) printf("(DEBUGGING INFO - TOV_RHSs()) Radius = %.15e | P = %.15e | M = %.15e | mu = %.15e\n",rr,P,m,mu);

  if( rr > 0.0 && P > 0.0 ) {
    /* Compute dP/dr */
    RHS[PRESSURE] = -( mu + P )*( m + 4 * M_PI * r3 * P ) / ( r2*( 1.0 - 2*m/rr ) );

    /* Compute dnu/dr */
    RHS[NU] = - ( 2.0 / (mu + P) ) * RHS[PRESSURE];

    /* Compute dm/dr */
    RHS[MASS] = 4.0 * M_PI * r2 * mu;

    /* Debugging output */
    if(DEBUG) printf("(DEBUGGING INFO - TOV_RHSs()) dP/dr = %.15e\n",RHS[PRESSURE]);
    if(DEBUG) printf("(DEBUGGING INFO - TOV_RHSs()) dnu/dr = %.15e\n",RHS[NU]);
    if(DEBUG) printf("(DEBUGGING INFO - TOV_RHSs()) dm/dr = %.15e\n",RHS[MASS]);
  }
  else {
    RHS[PRESSURE] = 0.0;
    RHS[NU]       = 0.0;
    RHS[MASS]     = 0.0;
  }

}
