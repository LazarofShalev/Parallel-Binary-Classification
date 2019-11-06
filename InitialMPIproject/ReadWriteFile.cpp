#define _CRT_SECURE_NO_WARNINGS

// Final Project - BinaryClassification
// Shalev Lazarof - 300209202

#include <mpi.h>
#include "StructsAndFuncsDecl.h"

Points readInputFile(FileProperties* fileProps, const char* path) 
{
	FILE* f = fopen(path, "r");
	Points points;

	// check the opened Correctly
	if (f == NULL)
	{
		printf("File is missing");
		exit(0);
	}

	// scan the first line in the input file
	fscanf(f, "%d%d%lf%lf%d%lf\n", &fileProps->numberOfPoints, &fileProps->coords, &fileProps->alpha0, &fileProps->alphaMax, &fileProps->limit, &fileProps->qc);

	points.X = (double*)calloc(fileProps->numberOfPoints * fileProps->coords, sizeof(double));
	points.classify = (int*)calloc(fileProps->numberOfPoints, sizeof(int));
	points.numberOfPoints = fileProps->numberOfPoints;
	points.coords = fileProps->coords;

	// for loop on the file input points number
	for (int i = 0; i < fileProps->numberOfPoints; i++)
	{
		for (int j = 0; j < fileProps->coords; j++) 
		{
			fscanf(f, "%lf ", &points.X[(i * fileProps->coords) + j]);
		}
		fscanf(f, "%d", &points.classify[i]);
	}

	// close the input file
	fclose(f);
	return points;
}

void writeOutputFile(double alpha, double alphaMax, double q, double* weightsVector, double w0, int coords, const char* path) 
{
	FILE* f = fopen(path, "w");

	// check the opened Correctly
	if (f == NULL) 
	{
		printf("File is missing");
		exit(0);
	}

	// check alphaMax <= alpha condition
	// write to output file
	if (alphaMax <= alpha) 
	{
		printf("\ncouldn't find Alpha\n");
		fprintf(f, "couldn't find Alpha\n");
	}

	// write Alpha minimum to output file
	printf("Alpha minimum = %lf, q = %lf\n\n", alpha, q);
	fprintf(f, "Alpha minimum = %lf, q = %lf\n\n", alpha, q);

	// write w0 to output file
	printf("w0 = %lf\n", w0);
	fprintf(f, "w0 = %lf\n", w0);

	// write w1,w2,.. to output file
	for (int i = 0; i < coords; i++)
	{
		fprintf(f, "w%d = %lf\n", (i + 1), weightsVector[i]);
		printf("w%d = %lf\n", (i + 1), weightsVector[i]);
	}

	// close the output file
	fclose(f);
}