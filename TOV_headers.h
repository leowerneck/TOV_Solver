/* .-----------------------------------------.
 * |  Copyright (c) 2019, Leonardo Werneck   |
 * | Licensed under the BSD 2-Clause License |
 * .-----------------------------------------.
 */

/* Program     : TOV Solver
 * File        : tov_headers.h
 * Author      : Leo Werneck (werneck@if.usp.br)
 * Date        : October 29, 2019
 *
 * Description : This file contain the necessary function and variable
 *               definitions for the TOV Solver.
 *
 * Dependencies: None
 *
 * Reference(s): Read et al., PRD 79, 124032 (2009) | (https://arxiv.org/pdf/0812.2163.pdf)
 */
#ifndef __TOV_HEADERS_h__
#define __TOV_HEADERS_h__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>

using namespace std;

/* .---------.
 * | Defines |
 * .---------.
 */

#ifdef __cplusplus
#define restrict __restrict__
#endif

/* Debugging option */
#define DEBUG 0

/* Definition of REAL */
#define REAL double

/* Physical constants */
#define SPEEDOFLIGHT 2.99792e+10
#define GNEWTON 6.67408e-08
#define SUNMASS 1.989e+33

/* Gridfunction definitions */
#define NGFS 3
#define PRESSURE 0
#define NU       1
#define MASS     2

/* EOS related definitions */
#define MAXEOSNAMESIZE 6
#define MAXEOSPARAMETERS 7

/* Integration related definitions */
#define MAXIMUMSTEPSIZE 1e-2
#define MAXINTEGRATIONS 300
#define MIN(a,b) a < b ? a : b

/* .------------.
 * | EOS Struct |
 * .------------.
 */
struct eos_struct
{
  int neos;
  REAL rho_b_PPEOS[MAXEOSPARAMETERS];
  REAL Gamma_PPEOS[MAXEOSPARAMETERS];
  REAL Kpoly_PPEOS[MAXEOSPARAMETERS];
  REAL Press_PPEOS[MAXEOSPARAMETERS];
  REAL epsIC_PPEOS[MAXEOSPARAMETERS];
};

/*
 * .---------------------.
 * | Function prototypes |
 * .---------------------.
 */
void RK4_method( eos_struct eos, const REAL rr, const REAL dr,  const REAL *restrict gfs_at_rr, REAL *restrict gfs_at_rr_plus_dr );
void TOV_RHSs( const eos_struct eos, const REAL rr, const REAL *restrict gfs_at_rr, REAL *restrict RHS );
void print_TOV_logo(void);
void print_EOS_table( const eos_struct eos );
REAL convert_mass_to_solar_masses( REAL Mass_G_units );
REAL convert_radius_to_kilometers( REAL Radius_G_units );
void output_rhob_mass_radius_to_file( string EOSname, const int iteration, const REAL rho_b_central, const REAL Mass_in_solar_masses, const REAL Radius_in_kilometers );
void output_current_data_to_file(string outfile_name,
                                 eos_struct eos,
                                 int iteration,
                                 REAL Radius_G_units,
                                 REAL step_size_G_units,
                                 REAL *restrict gfs_at_rr);
int polytropic_index_from_rhob( eos_struct eos, REAL rho_b_in );
int polytropic_index_from_pressure( eos_struct eos, REAL P_in );
REAL compute_P_from_rhob( eos_struct eos, REAL rho_b_in );
REAL compute_rhob_from_P( eos_struct eos, REAL P_in );
REAL compute_mu_from_P( eos_struct eos, REAL P_in );
void initialize_EOS_struct( string EOSname, eos_struct &eos );
void convert_EOS_struct_to_geometrized_units( eos_struct &eos );
void populate_Kpoly_PPEOS( eos_struct &eos );
void populate_Press_PPEOS( eos_struct &eos );
void populate_epsIC_PPEOS( eos_struct &eos );
int get_EOS_parameters_from_EOSname( string EOSname, eos_struct &eos );
void recompute_step_size( const eos_struct eos, const REAL rr, const REAL *gfs_at_rr, REAL &dr );


#endif // __TOV_HEADERS_h__
