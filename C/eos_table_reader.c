#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>

typedef struct eos_table_t {
  int nr, nt, ny;
  double energy_shift;
  double *lr, *lt, *ye, *lp, *le, *munu;
} eos_table_t;

static inline
void read_table_ints( const hid_t fid, const char *name, int *data ) {

  hid_t dataset_id   = H5Dopen2(fid, name, H5P_DEFAULT);
  hid_t dataspace_id = H5Dget_space(dataset_id);
  if( H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0 )
    fprintf(stderr, "Failed to read dataset %s from EOS table\n", name);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);
}

static inline
void read_table_doubles( const hid_t fid, const char *name, double *data ) {

  hid_t dataset_id   = H5Dopen2(fid, name, H5P_DEFAULT);
  hid_t dataspace_id = H5Dget_space(dataset_id);
  if( H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0 )
    fprintf(stderr, "Failed to read dataset %s from EOS table\n", name);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);
}

void read_table( const char *filename, eos_table_t *eos_table ) {

  int nr, nt, ny;
  double energy_shift;
  double *lr, *lt, *ye, *lp, *le, *munu;
  hid_t file_id;

  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  read_table_ints(file_id, "pointsrho" , &nr);
  read_table_ints(file_id, "pointstemp", &nt);
  read_table_ints(file_id, "pointsye"  , &ny);

  lr   = malloc(sizeof(double)*nr);
  lt   = malloc(sizeof(double)*nt);
  ye   = malloc(sizeof(double)*ny);
  lp   = malloc(sizeof(double)*nr*nt*ny);
  le   = malloc(sizeof(double)*nr*nt*ny);
  munu = malloc(sizeof(double)*nr*nt*ny);

  read_table_doubles(file_id, "energy_shift", &energy_shift);
  read_table_doubles(file_id, "logrho"      , lr);
  read_table_doubles(file_id, "logtemp"     , lt);
  read_table_doubles(file_id, "ye"          , ye);
  read_table_doubles(file_id, "logpress"    , lp);
  read_table_doubles(file_id, "logenergy"   , le);
  read_table_doubles(file_id, "munu"        , munu);

  double log_of_10 = log(10.0);
  for(int i=0;i<nr;i++) lr[i] *= log_of_10;
  for(int j=0;j<nt;j++) lt[j] *= log_of_10;
  for(int n=0;n<nr*nt*ny;n++) {
    lp[n] *= log_of_10;
    le[n] *= log_of_10;
  }

  H5Fclose(file_id);

  eos_table->energy_shift = energy_shift;
  eos_table->nr   = nr;
  eos_table->nt   = nt;
  eos_table->ny   = ny;
  eos_table->lr   = lr;
  eos_table->lt   = lt;
  eos_table->ye   = ye;
  eos_table->lp   = lp;
  eos_table->le   = le;
  eos_table->munu = munu;
}

int main( int argc, char **argv ) {
  eos_table_t eos_table;
  read_table(argv[1], &eos_table);

  printf("(%d, %d, %d)\n", eos_table.nr, eos_table.nt, eos_table.ny);

  free(eos_table.lr);
  free(eos_table.lt);
  free(eos_table.ye);
  free(eos_table.lp);
  free(eos_table.le);
  free(eos_table.munu);
  return 0;
}
