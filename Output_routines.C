#ifndef __OUTPUT_ROUTINES__
#define __OUTPUT_ROUTINES__
/* Program     : TOV Solver
 * File        : Output_routines.C
 * Author      : Leo Werneck (werneck@if.usp.br)
 * Date        : October 29, 2019
 *
 * Description : This file implements functions to be used
 *               when outputting our results to file
 *
 * Dependencies: stdio.h, stdlib.h, string.h, TOV_headers.h, & Polytropic_EOS__lowlevel_functions.C
 *
 * Reference(s): None
 *
 */

/* .--------------.
 * | Dependencies |
 * .--------------.
 */

/* Standard C libraries */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* TOV_Solver files */
#include "TOV_headers.h"
#include "Polytropic_EOS__lowlevel_functions.C"

/* .---------------------.
 * | Function prototypes |
 * .---------------------.
 */
inline REAL convert_mass_to_solar_masses( REAL Mass_G_units );
inline REAL convert_radius_to_kilometers( REAL Radius_G_units );


/* Function   : output_program_logo()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Outputs the TOV logo
 *
 * Input(s)   : None
 *
 * Outputs(s) : Prints the program logo to screen
 *
 */
inline void print_TOV_logo() {

  printf(".------------.  .----------. .-.            .-.   \n");
  printf("|            |  |          |  \\ \\          / /  \n");
  printf(".----.  .----.  |   .--.   |   \\ \\        / /   \n");
  printf("     |  |       |   |  |   |    \\ \\      / /    \n");
  printf("     |  |       |   .--.   |     \\ \\    / /     \n");
  printf("     |  |       |          |      \\ \\  / /      \n");
  printf("     .--.       .----------.       \\_\\/_/       \n");
  printf("                                                  \n");
  printf(" TOV Solver by Leo Werneck (werneck@if.usp.br)  \n\n");

}

/* Function   : output_current_data_to_file()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Outputs current data to file
 *
 * Input(s)   : outfile_name      - Output file name
 *            : eos               - Struct containing EOS parameters
 *            : iteration         - Current iteration
 *            : Radius_G_units    - Current value of the radius (Geometrized units)
 *            : step_size_G_units - Current value of the step size (Geometrized units)
 *            : gfs_at_rr         - Gridfunctions at the current radius (Geometrized units)
 *
 * Outputs(s) : output file:
 *                           name - EOSname-TOV_Solver_Output.dat
 *                       Column 1 - Iteration number
 *                       Column 2 - Step size   (Geometrized units)
 *                       Column 3 - Radius      (Geometrized units)
 *                       Column 4 - Mass        (Geometrized units)
 *                       Column 5 - rho_baryon  (Geometrized units)
 *                       Column 6 - Pressure    (Geometrized units)
 *                       Column 7 - Mass        (Solar masses)
 *                       Column 8 - Radius      (Kilometers)
 */
inline void output_current_data_to_file(string outfile_name,eos_struct eos,
              				int iteration, REAL Radius_G_units, REAL step_size_G_units, REAL *gfs_at_rr) {

  /* Pressure in geometrized units */
  REAL Pressure = gfs_at_rr[PRESSURE];

  /* Baryonic density in geometrized units */
  REAL rho_baryon = compute_rhob_from_P(eos,Pressure);

  /* Declare the mass and radius in geometrized units */
  REAL Mass_G_units   = gfs_at_rr[MASS];

  /* Declare the mass in solar masses andd the radius in kilometers */
  REAL Mass_solar_masses = convert_mass_to_solar_masses( Mass_G_units );
  REAL Radius_km         = convert_radius_to_kilometers( Radius_G_units );

  /* If we are printing initial data, create the
   * output file, otherwise append to file
   */
  FILE *outfileTOV;
  if( iteration == 0 ) {
    outfileTOV = fopen(outfile_name.c_str(),"w");
  }
  else {
    outfileTOV = fopen(outfile_name.c_str(),"a");
  }

  /* .----------------.
   * | Output to file |
   * .----------------.
   */
  fprintf(outfileTOV,"%d %e %e %e %e %e %e %e\n",
          iteration,step_size_G_units,Radius_G_units,Mass_G_units,rho_baryon,Pressure,Mass_solar_masses,Radius_km
	 );


  fclose(outfileTOV);
}

/* Function   : print_EOS_table()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 29, 2019
 *
 * Description: Prints the EOS table input by the user
 *
 * Input(s)   : eos - Struct containing EOS parameters
 *
 * Outputs(s) : Prints EOS table to terminal
 *
 */
inline void print_EOS_table( const eos_struct eos ) {

  printf(".--------------------.-----------------.\n");
  printf("|   EOS Parameter    | Parameter Value |\n");
  printf(".--------------------.-----------------.\n");
  printf("|     eos.neos       |        %d        |\n",eos.neos);
  printf(".--------------------.-----------------.\n");
  for(int j=0; j<eos.neos-1; j++) {
    printf("| eos.rho_b_PPEOS[%d] | %.9e |\n",j,eos.rho_b_PPEOS[j]);
  }
  if(eos.neos != 1) {
    printf(".--------------------.-----------------.\n");
  }
  for(int j=0; j<eos.neos; j++) {
    printf("| eos.Kpoly_PPEOS[%d] | %.9e |\n",j,eos.Kpoly_PPEOS[j]);
  }
  printf(".--------------------.-----------------.\n");
  for(int j=0; j<eos.neos; j++) {
    printf("| eos.Gamma_PPEOS[%d] | %.9e |\n",j,eos.Gamma_PPEOS[j]);
  }
  printf(".--------------------.-----------------.\n");

}

/* Function   : convert_mass_to_solar_masses()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 30, 2019
 *
 * Description: Converts mass from code units to solar massses
 *
 * Input(s)   : Mass_G_units      - Mass in code units
 *
 * Outputs(s) : Mass_solar_masses - Mass in solar masses
 *
 */
inline REAL convert_mass_to_solar_masses( REAL Mass_G_units ) {

  /* .--------------------.
   * | Physical constants |
   * .--------------------.
   *
   * Speed of light in cgs (cm/s) units */
  REAL c = SPEEDOFLIGHT;

  /* Gravitational constant in cgs (cm^3/g/s^2) */
  REAL G_Newton = GNEWTON;

  /* Mass of the Sun in cgs (grams) */
  REAL Sun_mass = SUNMASS;

  /* We have chosen units such that
   *
   * 1 <-> 1e15 g/cm^3 .
   *
   * Now, because of the units of G and c, namely
   *
   * [G] = cm^3/g/s^2 , [c] = cm/s ,
   *
   * it is easy to see that the combination
   * .----------------.
   * | [G/c^2] = cm/g | ,
   * .----------------.
   * offers a converter back to cgs units.
   */
  REAL one_to_one_over_cm2 = 1.0e+15 * G_Newton / c / c;
  REAL one_to_one_over_cm  = sqrt(one_to_one_over_cm2);
  REAL one_to_cm           = 1.0/one_to_one_over_cm;
  REAL one_to_g            = one_to_cm * c * c / G_Newton;

  /* To obtain the current mass in solar masses,
   * we need to first convert it to grams, then
   * simply rewrite it in Solar masses units.
   * .-------------------------------------------.
   * | M_grams        = (one_to_g) * M_G_units   |
   * | M_Solar_masses = M_grams / M_Sun_in_grams |
   * .-------------------------------------------.
   */
  return(one_to_g * Mass_G_units / Sun_mass);

}

/* Function   : convert_radius_to_kilometers()
 * Author     : Leo Werneck (werneck@if.usp.br)
 * Date       : October 30, 2019
 *
 * Description: Converts radius from code units to kilometers
 *
 * Input(s)   : Radius_G_units - Radius in code units
 *
 * Outputs(s) : Radius_km      - Radius in kilometers
 */
inline REAL convert_radius_to_kilometers( REAL Radius_G_units ) {

  /* .--------------------.
   * | Physical constants |
   * .--------------------.
   *
   * Speed of light in cgs (cm/s) units */
  REAL c = SPEEDOFLIGHT;

  /* Gravitational constant in cgs (cm^3/g/s^2) */
  REAL G_Newton = GNEWTON;

  /* We have chosen units such that
   *
   * 1 <-> 1e15 g/cm^3 ,
   *
   * and with G = 1 = c, we have
   *
   * 1 <-> 1e15 1/cm^2 .
   *
   * Now, because of the units of G and c, namely
   *
   * [G] = cm^3/g/s^2 , [c] = cm/s ,
   *
   * it is easy to see that the combination
   * .----------------.
   * | [G/c^2] = cm/g | ,
   * .----------------.
   * offers a converter back to cgs units.
   */
  REAL one_to_one_over_cm2 = 1.0e+15 * G_Newton / c / c;
  REAL one_to_one_over_cm  = sqrt(one_to_one_over_cm2);
  REAL one_to_cm           = 1.0/one_to_one_over_cm;
  REAL one_to_km           = one_to_cm * 1.0e-5;

  /* Finally, in geometrized units we still have
   * [Radius] = Length = cm (because we started
   * with cgs units). Thus
   * .------------------------------------------.
   * | Radius_km = (one_to_km) * Radius_G_units |
   * .------------------------------------------.
   */
  return(Radius_G_units * one_to_km);

}

inline void output_rhob_mass_radius_to_file( string EOSname, const int iteration, const REAL rho_b_central, const REAL Mass_in_solar_masses, const REAL Radius_in_kilometers ) {

  string filename = EOSname+"_mass_vs_radius.dat";
  FILE *outfilemassvsradius;

  if( iteration == 0 ) {
    outfilemassvsradius = fopen(filename.c_str(),"w");
  }
  else {
    outfilemassvsradius = fopen(filename.c_str(),"a");
  }

  fprintf(outfilemassvsradius,"%d %.15e %.15e %.15e\n",iteration,rho_b_central,Mass_in_solar_masses,Radius_in_kilometers);

  fclose(outfilemassvsradius);

}

#endif
