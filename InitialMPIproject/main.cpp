#define _CRT_SECURE_NO_WARNINGS

// Final Project - BinaryClassification
// Shalev Lazarof - 300209202

#include <cmath>
#include <omp.h>
#include "mpi.h"

#include "StructsAndFuncsDecl.h"

#define INPUT "C:\\Users\\Lazar\\Downloads\\data\\data1.txt"
#define OUTPUT "C:\\Users\\Lazar\\Downloads\\data\\output.txt"

// calculate Weight Vector
void calculateWeightVector(int coords, double alpha, int type, double *weightsVector, double *w0, double *points) 
{
#pragma omp parallel for
	for (int i = 0; i < coords; ++i) 
	{
		weightsVector[i] = weightsVector[i] + ((alpha * type) * points[i]);
	}
	*w0 = *w0 + (alpha * type);
}

// calculate q value
double CalculateQ(double* sums, Points* points) 
{
	int cnt = 0;
#pragma omp parallel for reduction(+:cnt)
	for (int i = 0; i < points->numberOfPoints; i++) 
	{
		int type;

		// update type value according to the sum value
		if (0 < sums[i]) {
			type = 1;
		}
		else {
			type = -1;
		}

		// update cnt
		if (points->classify[i] != type) {
			cnt++;
		}
	}
	return (double)cnt / points->numberOfPoints;
}

// calculate Sum
double calculateSum(int coords, double w0, double *weightsVector, double *points)
{
	double sum = 0;
	for (int i = 0; i < coords; i++)
	{
		sum = sum + (points[i] * weightsVector[i]);
	}
	sum = sum + w0;
	return sum;
}

void procsGather(int numberOfProcs, int rank, double *qa, double *qaPairs, double *w0, double *weightsVector, int *MinAlphaProc, FileProperties *fileProps) 
{
	int coords = fileProps->coords;

	*MinAlphaProc = -1;

	// All procs gather
	MPI_Allgather(qa, 2, MPI_DOUBLE, qaPairs, 2, MPI_DOUBLE, MPI_COMM_WORLD);

	for (int i = 0; i < numberOfProcs; i++)
	{
		if ((fileProps->alphaMax - fileProps->alpha0 < qaPairs[(i * 2) + 1]) || (qaPairs[i * 2] < fileProps->qc))
		{
			*MinAlphaProc = i;
			break;
		}
	}

	if (*MinAlphaProc != -1 && *MinAlphaProc != 0) 
	{
		// slaves send
		if (rank == *MinAlphaProc) 
		{
			MPI_Send(w0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(weightsVector, coords, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		// master receive 
		else if (rank == 0) 
		{
			MPI_Recv(w0, 1, MPI_DOUBLE, *MinAlphaProc, 0, MPI_COMM_WORLD, 0);
			MPI_Recv(weightsVector, coords, MPI_DOUBLE, *MinAlphaProc, 0, MPI_COMM_WORLD, 0);
		}
	}
}

bool classify(int numberOfPoints, double alpha, double *weightsVector, double *w0, double *sums, Points *points) 
{
	int type;
	int coords = points->coords;

#pragma omp parallel for
	for (int i = 0; i < numberOfPoints; ++i) 
	{
		sums[i] = calculateSum(coords, *w0, weightsVector, points->X + (points->coords * i));
	}

	for (int i = 0; i < numberOfPoints; i++)
	{
		// update type value according to the sum value
		if (0 < sums[i]) {
			type = 1;
		}
		else  {
			type = -1;
		}

		// each iteration check if need to update weight vector
		if (type != points->classify[i]) 
		{
			calculateWeightVector(coords, alpha, points->classify[i], weightsVector, w0, points->X + (points->coords * i));
			return true;
		}
	}
	return false;
}

void BinaryClassification(int rank, int numberOfProcs, FileProperties *fileProps, Points *points)
{
	bool flag;

	double q;
	double w0;
	double alpha;
	double minimumAlpha = -1;

	int MinAlphaProc;
	int coords;

	double* weightsVector;
	double* qaPairs;
	double qa[2];
	double* sums;

	alpha = fileProps->alpha0 * (1 + rank);
	coords = fileProps->coords;
	weightsVector = (double*)calloc(coords, sizeof(double));
	qaPairs = (double*)calloc(numberOfProcs * 2, sizeof(double));
	sums = (double*)calloc(fileProps->numberOfPoints, sizeof(double));

	MinAlphaProc = -1;
	w0 = 0;

	// master prints qc and alpha max to console
	if (rank == 0)
	{
		printf("qc = %lf alphaMax = %lf\n", fileProps->qc, fileProps->alphaMax);
	}

	do 
	{
		flag = true;
		for (int i = 0; i < (fileProps->limit) && (flag == true); i++)
		{
			flag = classify(fileProps->numberOfPoints, alpha, weightsVector, &w0, sums, points);
		}

		q = CalculateQ(sums, points);
		qa[0] = q;
		qa[1] = alpha;
		alpha = alpha + (fileProps->alpha0 * numberOfProcs);

		// all procs gather
		procsGather(numberOfProcs, rank, qa, qaPairs, &w0, weightsVector, &MinAlphaProc, fileProps);

	} while ((qa[1] < fileProps->alphaMax) && (fileProps->qc < q));
	alpha = alpha - (numberOfProcs * fileProps->alpha0);

	// master update and write output file
	if (rank == 0) 
	{
		// update minimumAlpha value based on MinAlphaProc
		if (MinAlphaProc == -1) {
			minimumAlpha = alpha;
		}
		else {
			minimumAlpha = qaPairs[MinAlphaProc * 2 + 1];
		}

		// update q value based on MinAlphaProc value
		if (MinAlphaProc == -1) {
			q = q;
		}
		else {
			q = qaPairs[MinAlphaProc * 2];
		}

		// master write output file
		writeOutputFile(minimumAlpha, fileProps->alphaMax, q, weightsVector, w0, coords, OUTPUT);
	}
}

int main(int argc, char* argv[]) 
{
	int rank;
	int numberOfProcs;

	double time1;
	double time2;

	FileProperties fileProps;
	Points points;

	// 1,1,1,1,1,1 because has only one of each
	int blocklen[] = { 1, 1, 1, 1, 1, 1, 1 };

	// define a new MPI_TYPE name
	MPI_Datatype PropsMPIType;

	// the data type of the specifiec struct
	MPI_Datatype type[] = { MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE };

	// the locations
	MPI_Aint disp[6];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcs);

	// all the locations
	disp[0] = offsetof(FileProperties, numberOfPoints);
	disp[1] = offsetof(FileProperties, coords);
	disp[2] = offsetof(FileProperties, alpha0);
	disp[3] = offsetof(FileProperties, alphaMax);
	disp[4] = offsetof(FileProperties, limit);
	disp[5] = offsetof(FileProperties, qc);

	// 6 is size of veriables, blocklen is the ammount of veriables in the struct, 
	// disp is the locations of the veriables, type is all types, &PropsMPIType the MPI_NEWTYPE address
	MPI_Type_create_struct(6, blocklen, disp, type, &PropsMPIType);

	// the new MPI_TYPE
	MPI_Type_commit(&PropsMPIType);

	// master read input file
	if (rank == 0) 
	{
		points = readInputFile(&fileProps, INPUT);
	}

	// broadcast
	MPI_Bcast(&fileProps, 1, PropsMPIType, 0, MPI_COMM_WORLD);

	// all slaves init points 
	if (rank != 0) 
	{
		points.numberOfPoints = fileProps.numberOfPoints;
		points.coords = fileProps.coords;
		points.classify = (int*)calloc(fileProps.numberOfPoints, sizeof(int));
		points.X = (double*)calloc(fileProps.numberOfPoints * fileProps.coords, sizeof(double));	
	}

	// broadcast
	MPI_Bcast(points.X, fileProps.numberOfPoints * fileProps.coords, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(points.classify, fileProps.numberOfPoints, MPI_INT, 0, MPI_COMM_WORLD);

	// measure times
	time1 = MPI_Wtime();
	BinaryClassification(rank, numberOfProcs, &fileProps, &points);
	time2 = MPI_Wtime();

	// master print time difference
	if (rank == 0)
	{
		printf("\nTime: %lf\n", time2 - time1);
		fflush(stdout);
	}

	MPI_Finalize();
	return 0;
}




