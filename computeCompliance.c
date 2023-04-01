#include <getopt.h>
#include <stdlib.h>

#include "io.h"
#include "libcompliance.h"

void setupLoad(const int nelx, const int nely, const int nelz,
               const double forceMagnitude, double *F) {
  const uint_fast32_t nx = nelx + 1;
  const uint_fast32_t ny = nely + 1;
  const uint_fast32_t nz = nelz + 1;

  if (true) { // cantilever
    const int k = 0;
    const double radius = ((double)ny) / 5.0; // snap
    const double radius2 = radius * radius;
    const double center_x = (double)nelx;
    const double center_y = ((double)nely - 1.0) / 2.0;

    int num_elements = 0;
    for (int i = 0; i < nelx; i++) {
      for (int j = 0; j < nely; j++) {
        const double dx = (double)i - center_x;
        const double dy = (double)j - center_y;
        const double dist2 = dx * dx + dy * dy;
        if (dist2 < radius2) {
          num_elements++;
        }
      }
    }

    double nodalForce = forceMagnitude / (4.0 * (double)num_elements);
    int isum = 0;
    for (int i = 0; i < nelx; i++) {
      for (int j = 0; j < nely; j++) {

        const double dx = (double)i - center_x;
        const double dy = (double)j - center_y;
        const double dist2 = dx * dx + dy * dy;

        if (dist2 < radius2) {
          const int nidx1 = (i + 1) * ny * nz + (j + 1) * nz + k;
          const int nidx2 = (i + 1) * ny * nz + j * nz + k;
          const int nidx3 = (i)*ny * nz + (j + 1) * nz + k;
          const int nidx4 = (i)*ny * nz + j * nz + k;
          F[3 * nidx1 + 2] += nodalForce;
          F[3 * nidx2 + 2] += nodalForce;
          F[3 * nidx3 + 2] += nodalForce;
          F[3 * nidx4 + 2] += nodalForce;
        }
      }
    }
  } else if (true) { // MBB Beam and snapthrough with center load

    const int k = nelz;
    const double radius = ((double)ny) / 5.0;
    const double radius2 = radius * radius;
    const double center_x = 0;
    const double center_y = ((double)nely - 1.0) / 2.0;

    int num_elements = 0;
    for (int i = 0; i < nelx; i++) {
      for (int j = 0; j < nely; j++) {
        const double dx = (double)i - center_x;
        const double dy = (double)j - center_y;
        const double dist2 = dx * dx + dy * dy;
        if (dist2 < radius2) {
          num_elements++;
        }
      }
    }

    double nodalForce = forceMagnitude / (4.0 * (double)num_elements);
    int isum = 0;
    for (int i = 0; i < nelx; i++) {
      for (int j = 0; j < nely; j++) {

        const double dx = (double)i - center_x;
        const double dy = (double)j - center_y;
        const double dist2 = dx * dx + dy * dy;

        if (dist2 < radius2) {
          const int nidx1 = (i + 1) * ny * nz + (j + 1) * nz + k;
          const int nidx2 = (i + 1) * ny * nz + j * nz + k;
          const int nidx3 = (i)*ny * nz + (j + 1) * nz + k;
          const int nidx4 = (i)*ny * nz + j * nz + k;
          F[3 * nidx1 + 2] += nodalForce;
          F[3 * nidx2 + 2] += nodalForce;
          F[3 * nidx3 + 2] += nodalForce;
          F[3 * nidx4 + 2] += nodalForce;
        }
      }
    }
  }
}

void handle_fut_err(struct futhark_context *fut_ctx, int err) {
  if (err != 0) {
    char *err_string = futhark_context_get_error(fut_ctx);
    printf("error in fut call!! error message:\n%s", err_string);
    free(err_string);
  }
}

int main(int argc, char *argv[]) {

  char *filename = NULL;
  char *outFilename = NULL;
  int opt;

  double forceMagnitude = -0.002;
  int nonlinear = 0;
  int zeroInitial = 0;

  while ((opt = getopt(argc, argv, "f:m:n:z:o:")) != -1) {
    switch (opt) {
    case 'f':
      filename = optarg;
      break;
    case 'o':
      outFilename = optarg;
      break;
    case 'm':
      forceMagnitude = atof(optarg);
      break;
    case 'n':
      nonlinear = atoi(optarg);
      break;
    case 'z':
      zeroInitial = atoi(optarg);
      break;
    default:
      fprintf(stderr, "Usage: %s [-xyzrv]\n", argv[0]);
      exit(1);
    }
  }

  if (filename == NULL) {
    fprintf(stderr, "Please specify a file with -f <filename>\n");
    exit(1);
  }

  if (outFilename == NULL) {
    fprintf(stderr, "Please specify a output file with -o <filename>\n");
    exit(1);
  }

  int nelx, nely, nelz;
  float *density;
  double *fileDisplacement;
  readDensity(&nelx, &nely, &nelz, &density, &fileDisplacement, filename);

  const int nx = nelx + 1;
  const int ny = nely + 1;
  const int nz = nelz + 1;

  const int ndof = 3 * nx * ny * nz;
  double *F = malloc(sizeof(double) * ndof);
  for (uint_fast64_t i = 0; i < ndof; i++) {
    F[i] = 0.0;
  }

  if (zeroInitial) { // set 0 initial condition
    printf("Setting 0 initial guess\n");
    for (uint_fast64_t i = 0; i < ndof; i++) {
      fileDisplacement[i] = 0.0;
    }
  }
  setupLoad(nelx, nely, nelz, forceMagnitude, F);

  struct futhark_context_config *fut_cfg = futhark_context_config_new();
  struct futhark_context *fut_ctx = futhark_context_new(fut_cfg);
  int err;

  struct futhark_f32_3d *xPhys_gpu =
      futhark_new_f32_3d(fut_ctx, density, nelx, nely, nelz);
  struct futhark_f64_4d *force_gpu =
      futhark_new_f64_4d(fut_ctx, F, nx, ny, nz, 3);
  struct futhark_f64_4d *u_init_gpu =
      futhark_new_f64_4d(fut_ctx, fileDisplacement, nx, ny, nz, 3);
  err = futhark_context_sync(fut_ctx);
  handle_fut_err(fut_ctx, err);

  int cgiter, linearIter;
  double cgres, compliance;

  if (nonlinear) { // compute nonlinear U
    printf("Computing nonlinear compliance\n");
    fflush(stdout);

    struct futhark_f64_4d *u_step_gpu;

    int nstep = 10;
    for (int i = 1; i < nstep + 1; i++) {
      const double loadIncrement = ((double)i) / ((double)nstep);

      err = futhark_entry_solveNonlinearLoadStep(
          fut_ctx, &u_step_gpu, &cgres, &cgiter, &linearIter, xPhys_gpu,
          u_init_gpu, force_gpu, loadIncrement);
      handle_fut_err(fut_ctx, err);
      err = futhark_context_sync(fut_ctx);
      handle_fut_err(fut_ctx, err);

      err = futhark_free_f64_4d(fut_ctx, u_init_gpu);
      handle_fut_err(fut_ctx, err);
      err = futhark_context_sync(fut_ctx);
      handle_fut_err(fut_ctx, err);

      u_init_gpu = u_step_gpu;

      printf("Solved loadstep %i with %e load fraction, to a relative "
             "tolerance of %e, in %i linear "
             "iterations "
             "(%i nonlinear)\n",
             i, loadIncrement, cgres, linearIter, cgiter);
      fflush(stdout);
    }

    err = futhark_entry_computeCompliance(fut_ctx, &compliance, u_init_gpu,
                                          force_gpu);
    handle_fut_err(fut_ctx, err);

  } else { // compute linear U
    printf("Computing linear compliance\n");
    fflush(stdout);

    struct futhark_f64_4d *u_result;

    err = futhark_entry_solveLinear(fut_ctx, &u_result, &cgres, &cgiter,
                                    &linearIter, xPhys_gpu, u_init_gpu,
                                    force_gpu);
    handle_fut_err(fut_ctx, err);
    err = futhark_context_sync(fut_ctx);
    handle_fut_err(fut_ctx, err);

    err = futhark_free_f64_4d(fut_ctx, u_init_gpu);
    handle_fut_err(fut_ctx, err);
    err = futhark_context_sync(fut_ctx);
    handle_fut_err(fut_ctx, err);

    u_init_gpu = u_result;

    err = futhark_entry_computeCompliance(fut_ctx, &compliance, u_init_gpu,
                                          force_gpu);
    handle_fut_err(fut_ctx, err);

    printf(
        "Solved system to a relative tolerance of %e, in %i linear iterations "
        "(%i nonlinear)\n",
        cgres, linearIter, cgiter);
  }
  err = futhark_context_sync(fut_ctx);
  handle_fut_err(fut_ctx, err);
  printf("Compliance is %e\n", compliance);

  // write displacement to disk
  err = futhark_values_f64_4d(fut_ctx, u_init_gpu, F);
  handle_fut_err(fut_ctx, err);
  err = futhark_context_sync(fut_ctx);
  handle_fut_err(fut_ctx, err);
  writeDensityAndDisplacementDouble(nelx, nely, nelz, density, F, outFilename);

  err = futhark_free_f64_4d(fut_ctx, u_init_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_free_f64_4d(fut_ctx, force_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_free_f32_3d(fut_ctx, xPhys_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_context_sync(fut_ctx);
  handle_fut_err(fut_ctx, err);

  futhark_context_free(fut_ctx);
  futhark_context_config_free(fut_cfg);

  free(fileDisplacement);
  free(density);
  free(F);

  return 0;
}