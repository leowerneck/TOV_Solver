/* .-----------------------------------------.
 * |  Copyright (c) 2019, Leonardo Werneck   |
 * | Licensed under the BSD 2-Clause License |
 * .-----------------------------------------.
 */

/* Program     : TOV Solver
 * File        : Polytropic_EOS__lowlevel_functions.C
 * Author      : Leo Werneck (werneck@if.usp.br)
 * Date        : October 29, 2019
 *
 * Description : This file implements many useful 'lowlevel' functions that
 *               handle single and piecewise polytropic EOSs.
 *
 * Dependencies: math.h & tov_headers.h
 *
 * Reference(s): Read et al., PRD 79, 124032 (2009) | (https://arxiv.org/pdf/0812.2163.pdf)
 *
 */
#include "tov_headers.h"

/* Function   : polytropic_index_from_rhob()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Determines the appropriate polytropic index from baryonic density
 *
 * Input(s)   : eos              - Struct containing EOS information
 *            : rho_b_in         - Value of rhob for which we want to determine the polytropic index
 *
 * Outputs(s) : polytropic_index - Index corresponding to the appropriate EOS to be used with rho_in
 */
int polytropic_index_from_rhob( eos_struct eos, REAL rho_b_in ) {

  if(eos.neos == 1) { return 0; }

  int polytropic_index = 0;
  for(int j=0; j<eos.neos-1; j++) { polytropic_index += ( rho_b_in > eos.rho_b_PPEOS[j] ); }

  return polytropic_index;

}

/* Function   : polytropic_index_from_pressure()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Determines the appropriate polytropic index from pressure
 *
 * Input(s)   : eos              - Struct containing EOS information
 *            : P_in             - Value of pressure for which we want to determine the polytropic index
 *
 * Outputs(s) : polytropic_index - Index corresponding to the appropriate EOS to be used with P_in
 */
int polytropic_index_from_pressure( eos_struct eos, REAL P_in ) {

  if(eos.neos == 1) { return 0; }

  int polytropic_index = 0;
  for(int j=0; j<eos.neos-1; j++) { polytropic_index += ( P_in > eos.Press_PPEOS[j] ); }

  return polytropic_index;

}

/* Function   : compute_P_from_rhob()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Computes P(rho) for a polytropic EOS
 *
 * Input(s)   : eos      - Struct containing EOS information
 *            : rho_b_in - Value of rho_b for which we want to compute the pressure
 *
 * Outputs(s) : Pressure - P(rho_b_in)
 */
REAL compute_P_from_rhob( eos_struct eos, REAL rho_b_in ) {

  if( rho_b_in <= 0.0 ) { return 0.0; }

  /* Compute the polytropic index from rho_b_in */
  int j = polytropic_index_from_rhob(eos,rho_b_in);

  /* Compute the polytropic presssure for the appropriate EOS:
   * .--------------------------------------.
   * | P(rho_b) = K_{j} * rho_b^(Gamma_{j}) |
   * .--------------------------------------.
   */
  return ( eos.Kpoly_PPEOS[j] * pow( rho_b_in, eos.Gamma_PPEOS[j] ) );

}

/* Function   : compute_rhob_from_P()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Computes rho(P) for a polytropic EOS
 *
 * Input(s)   : eos   - Struct containing EOS information
 *            : P_in  - Value of rho_b for which we want to compute the pressure
 *
 * Outputs(s) : rho_b - rho_b(P_in)
 */
REAL compute_rhob_from_P( eos_struct eos, REAL P_in ) {

  if( P_in <= 0.0 ) { return 0.0; }

  /* Compute the polytropic index from P_in */
  int j = polytropic_index_from_pressure(eos,P_in);

  /* Compute the polytropic presssure for the appropriate EOS:
   * .---------------------------------.
   * | rho_b = (P/K_{j})^(1/Gamma_{j}) |
   * .---------------------------------.
   */
  return ( pow( P_in/eos.Kpoly_PPEOS[j], 1.0/eos.Gamma_PPEOS[j] ) );

}

/* Function   : compute_mu_from_P()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Computes rho(P) for a polytropic EOS
 *
 * Input(s)   : eos  - Struct containing EOS information
 *            : P_in - Value of rho_b for which we want to compute the pressure
 *
 * Outputs(s) : mu   - mu(P_in)
 */
REAL compute_mu_from_P( eos_struct eos, REAL P_in ) {

  if( P_in <= 0.0 ) { return 0.0; }

  /* Compute rho_b from P */
  REAL rho_b = compute_rhob_from_P( eos, P_in );

  /* Compute the polytropic index from P_in */
  int j = polytropic_index_from_pressure(eos,P_in);

  /* Compute eps for a polytropic EOS:
   * .-------------------------------------.
   * | eps = epsIC + P / ( rho*(Gamma-1) ) |
   * .-------------------------------------.
   */
  REAL eps = eos.epsIC_PPEOS[j] + P_in / ( rho_b * (eos.Gamma_PPEOS[j] - 1.0) );

  /* Compute the total energy density
   * .----------------------.
   * | mu = (1 + eps)*rho_b |
   * .----------------------.
   */
  return ( (1.0 + eps)*rho_b );

}
