#ifndef EOS_H_
#define EOS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>

// .--------------------.
// | Physical constants |
// .--------------------.
#define SPEEDOFLIGHT (2.997924580000000e+10)

// .---------.
// | Structs |
// .---------.
typedef struct eos_poly_t {
  unsigned n_pieces;
  double *restrict poly_rho;
  double *restrict poly_press;
  double *restrict poly_k;
  double *restrict poly_gamma;
} eos_poly_t;

typedef struct eos_table_t {
  int nr, nt, ny;
  double energy_shift, inv_dlr, inv_dlt, inv_dye, inv_dlp_of_lr;
  double *restrict lr;
  double *restrict lt;
  double *restrict ye;
  double *restrict lp;
  double *restrict le;
  double *restrict munu;
  double *restrict lp_of_lr;
  double *restrict le_of_lr;
  double *restrict ye_of_lr;
} eos_table_t;

// .---------------------.
// | Function prototypes |
// .---------------------.
// From eos_poly_init.c

// From eos_table_reader.c
eos_table_t *read_eos_table( const char *restrict filename );
void free_eos_table( eos_table_t *restrict eos_table );

// From eos_table_interp1d.c
int eos_table_P_of_rho( const eos_table_t *restrict eos_table, const double rho, double *restrict P );
int eos_table_eps_of_rho( const eos_table_t *restrict eos_table, const double rho, double *restrict eps );
int eos_table_rho_of_P( const eos_table_t *restrict eos_table, const double P, double *restrict rho );

// From eos_table_interp3d.c
int eos_table_munu_of_rho_Ye_T(
      const eos_table_t *restrict eos_table,
      const double rho,
      const double ye,
      const double T,
      double *restrict munu );

// From eos_free.c
void eos_poly_free( eos_poly_t *restrict eos_poly );
void eos_table_free( eos_table_t *restrict eos_table );

#endif // EOS_H_
