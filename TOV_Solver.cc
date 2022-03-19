/* .-----------------------------------------.
 * |  Copyright (c) 2019, Leonardo Werneck   |
 * | Licensed under the BSD 2-Clause License |
 * .-----------------------------------------.
 */

/* Program     : TOV Solver
 * File        : TOV_Solver.C
 * Author      : Leo Werneck (werneck@if.usp.br)
 * Date        : October 29, 2019
 *
 * Description : This file implements the RHSs of the Tolman–Oppenheimer–Volkoff
 *               equations
 *
 * Dependencies: stdio.h, stdlib.h, math.h, tov_headers.h, RK4.C, TOV_RHSs.C, Polytropic_EOS__struct_initialization.C, & Polytropic_EOS__lowlevel_functions.C
 *
 * Reference(s): Read et al., PRD 79, 124032 (2009) | (https://arxiv.org/pdf/0812.2163.pdf)
 *               https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation
 *
 */

/* TOV_Solver files */
#include "tov_headers.h"

/* Function   : main()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: TOV Solver - main driver function
 *
 * Input(s)   : argc           - Number of command line input parameters (expected: EOSname)
 *            : argv           - Array containing the command line input parameters
 *
 * Outputs(s) : return value:
 *                           0 - Program finished with no errors
 *                           1 - Input error, wrong usage of the program
 *                           2 - Input error, invalid EOS name
 *
 *            : output file:
 *                        name - EOSname-TOV_Solver_Output.dat
 *                    Column 1 - Iteration number
 *                    Column 2 - Step size   (Geometrized units)
 *                    Column 3 - Radius      (Geometrized units)
 *                    Column 4 - Mass        (Geometrized units)
 *                    Column 5 - rho_baryon  (Geometrized units)
 *                    Column 6 - Pressure    (Geometrized units)
 *                    Column 7 - Mass        (Solar masses)
 *                    Column 8 - Radius      (Kilometers)
 */
int main (int argc, char *argv[]) {

  print_TOV_logo();

  /* .--------------------.
   * | Reading user input |
   * .--------------------.
   *
   * Verify correct usage of the program */
  if(argc != 3) {
    printf("Error! Incorrect usage. Run the program as: ./TOV_Solver EOSname CentralDensity\n");
    exit(1);
  }

  /* Read in the EOSname */
  string EOSname = argv[1];

  printf("Input EOS: %s\n",EOSname.c_str());
  printf("Initializing EOS parameters...\n");

   /* .----------------.
   * | EOS Parameters |
   * .----------------.
   *
   * Declare and initialize the EOS struct */
  eos_struct eos;
  initialize_EOS_struct( EOSname, eos );

  print_EOS_table( eos );

  /* .----------------------.
   * | Variable declaration |
   * .----------------------.
   *
   * Gridfunctions */
  REAL gfs_at_rr[NGFS], gfs_at_rr_plus_dr[NGFS];

  /* Radial variable */
  REAL rr = 0.0;

  /* RK4 initial step size */
  REAL dr = 1.0e-5;

  /* RK4 integration limiter */
  int max_integration_steps = (int)(1.0/dr+0.5);

  /* .---------------------------.
   * | Output initial parameters |
   * .---------------------------.
   *
   * Declare outputfile name */
  string outfile_name = EOSname+"TOV_Solver_Output.dat";

  /* Declare output checkpoint */
  int output_checkpoint = 1e6;

  /* Output initial data to file */
  output_current_data_to_file(outfile_name,eos, 0,rr,dr,gfs_at_rr);

  /* Print integration related information to the user */
  printf("(TOV_Solver INFO) Iteration: %5d/%5d | Radius = %9.6lf | Mass = %9.6lf | Compactness = %9.6lf\n",
  	 0,max_integration_steps,rr,gfs_at_rr_plus_dr[MASS],gfs_at_rr_plus_dr[MASS] == 0 ? 0.0 : gfs_at_rr_plus_dr[MASS]/rr);

  /* .----------------------.
   * | Mass vs. Radius data |
   * .----------------------.
   * We now start the algorithm which loops over central densities, solving
   * the TOV equations for a given EOS. For each value of the central density,
   * the algorithm below finds the maximum value of the mass and its respective
   * radius. The algorithm continues until Mass(rho_c + drho_c) < Mass(rho_c),
   * where rho_c is the central density and drho_c how much we increment it in
   * between iterations.
   *
   * Start by setting the initial central density */
  REAL rho_b_central = atof(argv[2]);

  /* Then set the increment in the central density
   * Experimental: setting to *fixed*, 1% of the initial
   * central density.
   */
  REAL drho_b_central = 0.05 * rho_b_central;

  /* Set up variabless to store the outputs */
  REAL Mass_in_solar_masses = 0.0;
  REAL Radius_in_kilometers = 0.0;
  REAL Mass, Radius;

  for(int k=0; k<500; k++) {

    /* Print integration related information to the user */
    printf("(TOV_Solver INFO) EOSname: %s | Interation: %3d | rho_b_central = %6.4lf\n",EOSname.c_str(),k+1,rho_b_central);

    /* .-------------.
     * | Integration |
     * .-------------.
     *
     * Set up initial conditions */

    /* Pressure associated with central density: this is
     * the initial condition for our RK4 integration method
     */
    REAL Press_central = compute_P_from_rhob( eos, rho_b_central );

    /* Mass associated with central density: because the
     * central density corresponds to the density at r=0,
     * this means that m(r=0) = 0.
     */
    REAL mass_central = 0.0;

    /* To complete the set of initial conditions, we also
     * set nu(r=0) = 0.
     */
    REAL nu_central = 0.0;

    /* Set the initial conditions */
    gfs_at_rr[PRESSURE] = Press_central;
    gfs_at_rr[MASS]     = mass_central;
    gfs_at_rr[NU]       = nu_central;

    /* Perform the first RK4 step manually */
    RK4_method( eos, rr, dr, gfs_at_rr, gfs_at_rr_plus_dr );
    rr += dr;

    /* Begin the integration loop  */
    for(int n=2; n < max_integration_steps; n++) {

      /* Update the gridfunctions */
      for(int GF=PRESSURE; GF < NGFS; GF++) {
	gfs_at_rr[GF] = gfs_at_rr_plus_dr[GF];
      }

      /* Call upon the adaptive step size routine */
      recompute_step_size(eos,rr,gfs_at_rr, dr);

      /* Perform an RK4 integration step */
      RK4_method( eos, rr, dr, gfs_at_rr, gfs_at_rr_plus_dr );
      rr += dr;

      /* Output to file, if we have reached the checkpoint */
      if( n%output_checkpoint == 0 ) {

	/* Print integration related information to the user */
	printf("(TOV_Solver INFO) Iteration: %5d/%5d | Radius = %9.6lf | Mass = %9.6lf | Compactness = %9.6lf\n",
	       n,max_integration_steps,rr,gfs_at_rr_plus_dr[MASS],gfs_at_rr_plus_dr[MASS]/rr);

	/* Output to the current file */
	output_current_data_to_file(outfile_name,eos, n,rr,dr,gfs_at_rr_plus_dr);

      }

      if( gfs_at_rr_plus_dr[MASS] <= gfs_at_rr[MASS] ) {
	       output_current_data_to_file(outfile_name,eos, n,rr,dr,gfs_at_rr_plus_dr);
         Mass   = convert_mass_to_solar_masses( gfs_at_rr[MASS] );
         Radius = convert_radius_to_kilometers( rr );
         printf("(TOV_Solver INFO) EOS: %s | Results: Mass = %.5e Solar masses | Radius = %.5e km\n",EOSname.c_str(),Mass,Radius);
         break;
      }
      else{
        Mass   = 0.0;
        Radius = 0.0;
      }
    }

    if( Mass > Mass_in_solar_masses ) {
      Mass_in_solar_masses = Mass;
      Radius_in_kilometers = Radius;
      //output_rhob_mass_radius_to_file( EOSname, k, rho_b_central, gfs_at_rr[MASS], rr );
      output_rhob_mass_radius_to_file( EOSname, k, rho_b_central, Mass_in_solar_masses, Radius_in_kilometers );

      rho_b_central += drho_b_central;
      rr = 0.0;
      dr = 1.0e-5;
      continue;
    }
    else{
      printf("(TOV_Solver INFO) Maximum mass found for EOS: %s | Mass_max = %6.4lf Solar masses | Radius(Mass_max) = %7.4lf km\n",
	     EOSname.c_str(),Mass_in_solar_masses,Radius_in_kilometers);
      //output_rhob_mass_radius_to_file( EOSname, k, rho_b_central, gfs_at_rr[MASS], rr );
      output_rhob_mass_radius_to_file( EOSname, k, rho_b_central, Mass_in_solar_masses, Radius_in_kilometers );
      break;
    }

  }

  return 0;

}
