#define _CRT_SECURE_NO_WARNINGS

// Final Project - BinaryClassification
// Shalev Lazarof - 300209202
#include <cstdio>
#include <cstdlib>
#include <mpi.h>

struct FileProperties 
{
	int numberOfPoints;
	int coords;
	double alpha0;
	double alphaMax;
	int limit;
	double qc;
}
typedef FileProperties;

struct Points
{
	int numberOfPoints;
	int coords;
	int *classify;
	double *X;

}
typedef Points;

Points readInputFile(FileProperties* fileProps, const char* path);
void writeOutputFile(double alpha, double alphaMax, double q, double* weightsVector, double w0, int Coords, const char* path);