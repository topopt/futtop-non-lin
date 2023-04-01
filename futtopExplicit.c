#include <math.h>
#include <stdalign.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "io.h"
#include "libexplicit.h"

#ifdef _OPENMP
#include <omp.h>
#endif

double getTime() {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}

void handle_fut_err(struct futhark_context *fut_ctx, int err) {
  if (err != 0) {
    char *err_string = futhark_context_get_error(fut_ctx);
    printf("error in fut call!! error message:\n%s", err_string);
    free(err_string);
  }
}

void explicitPrintTimesteps(const uint_fast32_t nelx, const uint_fast32_t nely,
                            const uint_fast32_t nelz) {

  int err;
  const uint_fast64_t nelem = nelx * nely * nelz;

  const uint_fast32_t nx = nelx + 1;
  const uint_fast32_t ny = nely + 1;
  const uint_fast32_t nz = nelz + 1;

  const uint_fast64_t ndof = 3 * nx * ny * nz;

  float *F = malloc(sizeof(float) * ndof);
  float *U = malloc(sizeof(float) * ndof);

  for (uint_fast64_t i = 0; i < ndof; i++) {
    F[i] = 0.0;
    U[i] = 0.0;
  }

  const double forceval = -0.1;
  for (int j = 0; j < ny; j++) {
    const int i = nx - 1;
    const int k = 0;
    const uint_fast32_t nidx = i * ny * nz + j * nz + k;
    F[3 * nidx + 2] = forceval;
  }

  // set corners
  {
    const int j1 = 0;
    const int j2 = ny - 1;
    const int i = nx - 1;
    const int k = 0;
    const uint_fast32_t nidx1 = i * ny * nz + j1 * nz + k;
    const uint_fast32_t nidx2 = i * ny * nz + j2 * nz + k;
    F[3 * nidx1 + 2] = forceval / 2;
    F[3 * nidx2 + 2] = forceval / 2;
  }

  float *x = malloc(sizeof(float) * nelem);
  for (uint_fast64_t i = 0; i < nelem; i++) {
    x[i] = 1.0;
  }

  unsigned int loop = 0;
  float change = 1;

#ifdef _OPENMP
  printf(" OpenMP enabled with %d threads\n", omp_get_max_threads());
#endif
  const double totalTime = getTime();

  struct futhark_context_config *fut_cfg = futhark_context_config_new();
  struct futhark_context *fut_ctx = futhark_context_new(fut_cfg);

  // setup forces and displacements on device
  struct futhark_f32_4d *force_gpu =
      futhark_new_f32_4d(fut_ctx, F, nx, ny, nz, 3);
  struct futhark_f32_4d *uold_gpu =
      futhark_new_f32_4d(fut_ctx, U, nx, ny, nz, 3);
  struct futhark_f32_4d *u_gpu = futhark_new_f32_4d(fut_ctx, U, nx, ny, nz, 3);

  // setup densities on device
  struct futhark_f32_3d *x_gpu =
      futhark_new_f32_3d(fut_ctx, x, nelx, nely, nelz);
  struct futhark_f32_3d *stress =
      futhark_new_f32_3d(fut_ctx, x, nelx, nely, nelz);

  struct futhark_f32_4d *Mass;
  err = futhark_entry_assembleMass(fut_ctx, &Mass, x_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_context_sync(fut_ctx);
  handle_fut_err(fut_ctx, err);

  /* %% START ITERATION */
  const int num_timesteps = 20000;
  int filenumber = 0;
  int file_freq = 100;
  for (int i = 0; i < num_timesteps; i++) {

    const double loopTime = getTime();

    // start writeout of last iterations data
    if (i % file_freq == 0) {
      err = futhark_values_f32_4d(fut_ctx, u_gpu, U);
      handle_fut_err(fut_ctx, err);
      err = futhark_context_sync(fut_ctx);
      handle_fut_err(fut_ctx, err);
    }

    struct futhark_f32_4d *unew_gpu;
    float loadscale;
    err = futhark_entry_TimeStepNL(fut_ctx, &loadscale, &unew_gpu, 2000, i,
                                   x_gpu, force_gpu, Mass, u_gpu, uold_gpu);
    handle_fut_err(fut_ctx, err);

    if (i % file_freq == 0) {
      futhark_entry_computeVMstress(fut_ctx, &stress, x_gpu, unew_gpu);
      handle_fut_err(fut_ctx, err);
      err = futhark_context_sync(fut_ctx);
      handle_fut_err(fut_ctx, err);

      err = futhark_values_f32_3d(fut_ctx, stress, x);
      handle_fut_err(fut_ctx, err);
      err = futhark_context_sync(fut_ctx);
      handle_fut_err(fut_ctx, err);

      char filename[50];
      sprintf(filename, "timeseries_%d.vtr", filenumber);
      filenumber++;
      writeDensityAndDisplacement(nelx, nely, nelz, x, U, filename);
    }

    err = futhark_context_sync(fut_ctx);
    handle_fut_err(fut_ctx, err);

    if (i % file_freq == 0) {
      float stopval;
      err = futhark_entry_stoppingCriteria(fut_ctx, &stopval, x_gpu, force_gpu,
                                           Mass, unew_gpu, u_gpu);
      handle_fut_err(fut_ctx, err);
      err = futhark_context_sync(fut_ctx);
      handle_fut_err(fut_ctx, err);
      printf("Step %i Stopval:%6.3e ls: %6.3e\n", i, stopval, loadscale);
    }

    // clear uneeded data
    err = futhark_free_f32_4d(fut_ctx, uold_gpu);
    handle_fut_err(fut_ctx, err);
    err = futhark_context_sync(fut_ctx);
    handle_fut_err(fut_ctx, err);
    uold_gpu = u_gpu;
    u_gpu = unew_gpu;

    // printf("It.:%4i Obj.:%6.3e Vol.:%6.3f ch.:%4.2f relres: %4.2e iters:
    // %4i\n",
    //        loop, compliance, vol, change, cgres, cgiter);
  }

  printf("End time: %6.3f \n", getTime() - totalTime);

  err = futhark_values_f32_4d(fut_ctx, u_gpu, U);
  handle_fut_err(fut_ctx, err);
  err = futhark_context_sync(fut_ctx);
  handle_fut_err(fut_ctx, err);

  err = futhark_free_f32_3d(fut_ctx, x_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_free_f32_3d(fut_ctx, stress);
  handle_fut_err(fut_ctx, err);
  err = futhark_free_f32_4d(fut_ctx, force_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_free_f32_4d(fut_ctx, uold_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_free_f32_4d(fut_ctx, u_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_free_f32_4d(fut_ctx, Mass);
  handle_fut_err(fut_ctx, err);
  err = futhark_context_sync(fut_ctx);
  handle_fut_err(fut_ctx, err);

  futhark_context_free(fut_ctx);
  futhark_context_config_free(fut_cfg);

  char filename[50];
  sprintf(filename, "timeseries_%d.vtr", filenumber);
  writeDensityAndDisplacement(nelx, nely, nelz, x, U, filename);

  free(x);
  free(F);
  free(U);
}

int main(int argc, char *argv[]) {
  int nelx_coarse = 4;
  int nely_coarse = 2;
  int nelz_coarse = 2;
  float volfrac = 0.2;
  float rmin = 1.5;
  int iters = 20;
  int opt;

  while ((opt = getopt(argc, argv, "x:y:z:r:v:i:")) != -1) {
    switch (opt) {
    case 'x':
      nelx_coarse = atoi(optarg);
      break;
    case 'y':
      nely_coarse = atoi(optarg);
      break;
    case 'z':
      nelz_coarse = atoi(optarg);
      break;
    case 'r':
      rmin = atof(optarg);
      break;
    case 'v':
      volfrac = atof(optarg);
      break;
    case 'i':
      iters = atoi(optarg);
      break;
    default:
      fprintf(stderr, "Usage: %s [-xyzrv]\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  int nelx = nelx_coarse * 8;
  int nely = nely_coarse * 8;
  int nelz = nelz_coarse * 8;

  printf("Running topopt with:\n number of coarse els x (-x): %i (fine els = "
         "%i)\n number of coarse els y (-y): %i (fine els = %i)\n number of "
         "coarse els z (-z): %i (fine els = %i)\n total number of elements: "
         "%i\n volume fraction (-v): %f\n filter radius in elements (-r): %f\n "
         "number of design iterations (-i): %i\n\n",
         nelx_coarse, nelx, nely_coarse, nely, nelz_coarse, nelz,
         nelx * nely * nelz, volfrac, rmin, iters);
  explicitPrintTimesteps(nelx, nely, nelz);
}
