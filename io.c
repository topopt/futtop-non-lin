#include "io.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// writes a file with a snapshot of the density field (x,xPhys), can be opened
// with paraview temperature: very cold, usually called once only
void writeDensity(const int nelx, const int nely, const int nelz,
                  const float *densityArray, const char *filename) {
  int nx = nelx + 1;
  int ny = nely + 1;
  int nz = nelz + 1;

  int numberOfNodes = nx * ny * nz;
  int numberOfElements = nelx * nely * nelz;

  FILE *fid = fopen(filename, "w");

  // write header
  fprintf(fid, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" "
               "byte_order=\"LittleEndian\">\n");
  fprintf(fid, "<RectilinearGrid WholeExtent=\"1 %i 1 %i 1 %i\">\n ", nx, ny,
          nz);
  fprintf(fid, "<Piece Extent=\"1 %i 1 %i 1 %i\">\n", nx, ny, nz);

  // points
  fprintf(fid, "<Coordinates>\n");
  fprintf(fid,
          "<DataArray type=\"Float32\" NumberOfComponents=\" %i\" "
          "format=\"ascii\">\n",
          1);
  for (int i = 0; i < nx; i++)
    fprintf(fid, "%e\n", (float)i);

  fprintf(fid, "</DataArray>\n");
  fprintf(fid,
          "<DataArray type=\"Float32\" NumberOfComponents=\"%i\" "
          "format=\"ascii\">\n",
          1);
  for (int j = 0; j < ny; j++)
    fprintf(fid, "%e\n", (float)j);
  fprintf(fid, "</DataArray>\n");

  fprintf(fid,
          "<DataArray type=\"Float32\" NumberOfComponents=\"%i\" "
          "format=\"ascii\">\n",
          1);
  for (int k = 0; k < nz; k++)
    fprintf(fid, "%e\n", (float)k);
  fprintf(fid, "</DataArray>\n");
  fprintf(fid, "</Coordinates>\n");

  fprintf(fid, "<CellData>\n");
  fprintf(fid, "<DataArray type=\"Float32\" NumberOfComponents=\"1\" "
               "Name=\"density\" format=\"ascii\">\n");
  for (int k = 0; k < nelz; k++)
    for (int j = 0; j < nely; j++)
      for (int i = 0; i < nelx; i++) {
        const int index = i * nely * nelz + j * nelz + k;
        fprintf(fid, "%e\n", densityArray[index]);
      }

  fprintf(fid, "</DataArray>\n");
  fprintf(fid, "</CellData>\n");

  fprintf(fid, "</Piece>\n");
  fprintf(fid, "</RectilinearGrid>\n");
  fprintf(fid, "</VTKFile>\n");

  fclose(fid);
}

// writes a file with a snapshot of the density field (x,xPhys) and the
// displacement field, can be opened with paraview temperature: cold,
// usually called once only
void writeDensityAndDisplacement(const int nelx, const int nely, const int nelz,
                                 const float *densityArray,
                                 const float *displacementArray,
                                 const char *filename) {
  int nx = nelx + 1;
  int ny = nely + 1;
  int nz = nelz + 1;

  int numberOfNodes = nx * ny * nz;
  int numberOfElements = nelx * nely * nelz;

  FILE *fid = fopen(filename, "w");

  // write header
  fprintf(fid, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" "
               "byte_order=\"LittleEndian\">\n");
  fprintf(fid, "<RectilinearGrid WholeExtent=\"1 %i 1 %i 1 %i\">\n ", nx, ny,
          nz);
  fprintf(fid, "<Piece Extent=\"1 %i 1 %i 1 %i\">\n", nx, ny, nz);

  // points
  fprintf(fid, "<Coordinates>\n");
  fprintf(fid,
          "<DataArray type=\"Float32\" NumberOfComponents=\" %i\" "
          "format=\"ascii\">\n",
          1);
  for (int i = 0; i < nx; i++)
    fprintf(fid, "%e\n", (float)i);

  fprintf(fid, "</DataArray>\n");
  fprintf(fid,
          "<DataArray type=\"Float32\" NumberOfComponents=\"%i\" "
          "format=\"ascii\">\n",
          1);
  for (int j = 0; j < ny; j++)
    fprintf(fid, "%e\n", (float)j);
  fprintf(fid, "</DataArray>\n");

  fprintf(fid,
          "<DataArray type=\"Float32\" NumberOfComponents=\"%i\" "
          "format=\"ascii\">\n",
          1);
  for (int k = 0; k < nz; k++)
    fprintf(fid, "%e\n", (float)k);
  fprintf(fid, "</DataArray>\n");
  fprintf(fid, "</Coordinates>\n");

  fprintf(fid, "<CellData>\n");
  fprintf(fid, "<DataArray type=\"Float32\" NumberOfComponents=\"1\" "
               "Name=\"density\" format=\"ascii\">\n");
  for (int k = 0; k < nelz; k++)
    for (int j = 0; j < nely; j++)
      for (int i = 0; i < nelx; i++) {
        const int index = i * nely * nelz + j * nelz + k;
        fprintf(fid, "%e\n", densityArray[index]);
      }
  fprintf(fid, "</DataArray>\n");
  fprintf(fid, "</CellData>\n");

  fprintf(fid, "<PointData>\n");
  fprintf(fid, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
               "Name=\"displacement\" format=\"ascii\">\n");
  for (int k = 0; k < nz; k++)
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++) {
        const int index = i * ny * nz * 3 + j * nz * 3 + k * 3;
        fprintf(fid, "%e %e %e\n", displacementArray[index],
                displacementArray[index + 1], displacementArray[index + 2]);
      }
  fprintf(fid, "</DataArray>\n");
  fprintf(fid, "</PointData>\n");

  fprintf(fid, "</Piece>\n");
  fprintf(fid, "</RectilinearGrid>\n");
  fprintf(fid, "</VTKFile>\n");

  fclose(fid);
}

// writes a file with a snapshot of the density field (x,xPhys) and the
// displacement field, can be opened with paraview temperature: cold,
// usually called once only
void writeDensityAndDisplacementDouble(const int nelx, const int nely,
                                       const int nelz,
                                       const float *densityArray,
                                       const double *displacementArray,
                                       const char *filename) {
  int nx = nelx + 1;
  int ny = nely + 1;
  int nz = nelz + 1;

  int numberOfNodes = nx * ny * nz;
  int numberOfElements = nelx * nely * nelz;

  float node_distance = 1.0 / ((float)nely);

  FILE *fid = fopen(filename, "w");

  // write header
  fprintf(fid, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" "
               "byte_order=\"LittleEndian\">\n");
  fprintf(fid, "<RectilinearGrid WholeExtent=\"1 %i 1 %i 1 %i\">\n ", nx, ny,
          nz);
  fprintf(fid, "<Piece Extent=\"1 %i 1 %i 1 %i\">\n", nx, ny, nz);

  // points
  fprintf(fid, "<Coordinates>\n");
  fprintf(fid,
          "<DataArray type=\"Float32\" NumberOfComponents=\" %i\" "
          "format=\"ascii\">\n",
          1);
  for (int i = 0; i < nx; i++)
    fprintf(fid, "%e\n", node_distance * i);

  fprintf(fid, "</DataArray>\n");
  fprintf(fid,
          "<DataArray type=\"Float32\" NumberOfComponents=\"%i\" "
          "format=\"ascii\">\n",
          1);
  for (int j = 0; j < ny; j++)
    fprintf(fid, "%e\n", node_distance * j);
  fprintf(fid, "</DataArray>\n");

  fprintf(fid,
          "<DataArray type=\"Float32\" NumberOfComponents=\"%i\" "
          "format=\"ascii\">\n",
          1);
  for (int k = 0; k < nz; k++)
    fprintf(fid, "%e\n", node_distance * k);
  fprintf(fid, "</DataArray>\n");
  fprintf(fid, "</Coordinates>\n");

  fprintf(fid, "<CellData>\n");
  fprintf(fid, "<DataArray type=\"Float32\" NumberOfComponents=\"1\" "
               "Name=\"density\" format=\"ascii\">\n");
  for (int k = 0; k < nelz; k++)
    for (int j = 0; j < nely; j++)
      for (int i = 0; i < nelx; i++) {
        const int index = i * nely * nelz + j * nelz + k;
        fprintf(fid, "%e\n", densityArray[index]);
      }
  fprintf(fid, "</DataArray>\n");
  fprintf(fid, "</CellData>\n");

  fprintf(fid, "<PointData>\n");
  fprintf(fid, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
               "Name=\"displacement\" format=\"ascii\">\n");
  for (int k = 0; k < nz; k++)
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++) {
        const int index = i * ny * nz * 3 + j * nz * 3 + k * 3;
        fprintf(fid, "%e %e %e\n", displacementArray[index],
                displacementArray[index + 1], displacementArray[index + 2]);
      }
  fprintf(fid, "</DataArray>\n");
  fprintf(fid, "</PointData>\n");

  fprintf(fid, "</Piece>\n");
  fprintf(fid, "</RectilinearGrid>\n");
  fprintf(fid, "</VTKFile>\n");

  fclose(fid);
}

void readDensity(int *nelx_out, int *nely_out, int *nelz_out,
                 float **densityArray, double **displacementArray,
                 const char *filename) {
  int nx, ny, nz, nelx, nely, nelz;

  FILE *fid = fopen(filename, "r");

  char lineBuffer[256];

  // run file through to find size
  while (fgets(lineBuffer, sizeof(lineBuffer), fid)) {
    if (strstr(lineBuffer, "WholeExtent") != NULL)
      break;
  }

  sscanf(lineBuffer, "<RectilinearGrid WholeExtent=\"1 %i 1 %i 1 %i\">\n ", &nx,
         &ny, &nz);

  nelx = nx - 1;
  nely = ny - 1;
  nelz = nz - 1;
  const int nelem = nelx * nely * nelz;
  const int ndof = 3 * nx * ny * nz;

  *nelx_out = nelx;
  *nely_out = nely;
  *nelz_out = nelz;
  *densityArray = malloc(sizeof(float) * nelem);
  *displacementArray = malloc(sizeof(double) * ndof);

  printf("Reading density of size [%i,%i,%i] (= %i) \n", nelx, nely, nelz,
         nelem);

  // run file through to find size
  while (fgets(lineBuffer, sizeof(lineBuffer), fid)) {
    if (strstr(lineBuffer, "CellData") != NULL) {
      fgets(lineBuffer, sizeof(lineBuffer), fid);
      break;
    }
  }

  for (int k = 0; k < nelz; k++)
    for (int j = 0; j < nely; j++)
      for (int i = 0; i < nelx; i++) {
        const int index = i * (nely) * (nelz) + j * (nelz) + k;
        fgets(lineBuffer, sizeof(lineBuffer), fid);
        sscanf(lineBuffer, "%e\n", (*densityArray) + index);
      }

  // run file through to find size
  while (fgets(lineBuffer, sizeof(lineBuffer), fid)) {
    if (strstr(lineBuffer, "PointData") != NULL) {
      fgets(lineBuffer, sizeof(lineBuffer), fid);
      break;
    }
  }

  for (int k = 0; k < nelz; k++)
    for (int j = 0; j < nely; j++)
      for (int i = 0; i < nelx; i++) {
        const int index = i * ny * nz * 3 + j * nz * 3 + k * 3;
        fgets(lineBuffer, sizeof(lineBuffer), fid);
        sscanf(lineBuffer, "%e %e %e\n", (*displacementArray) + index,
               (*displacementArray) + index + 1,
               (*displacementArray) + index + 2);
      }

  fclose(fid);
}