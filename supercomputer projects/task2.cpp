#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

#define A 2.0
#define B 2.0

double *XNodes, *YNodes;	

#define hx(i)  (XNodes[i+1]-XNodes[i])
#define hy(j)  (YNodes[j+1]-YNodes[j])
#define Max(A,B) ((A)>(B)?(A):(B))

inline double FCalc(double x,double y){
    return (2 * (x*x + y*y) * (1 - 2 * x*x * y*y) * exp(1 - x*x * y*y));
}

inline double Solution(double x, double y){
    return exp(1 - x*x * y*y);
}

inline double Borders(double x, double y){
    return exp(1 - x*x * y*y);
}

void FMatrCalc(double *FMatr, int Nx, int Ny){
    int i, j;
    for (i=0; i<Nx; i++)
        for (j=0; j<Ny; j++)
            FMatr[i*Ny+j] = FCalc(XNodes[i], YNodes[j]);
}

double Delta(double * P,int i, int j, int Nx, int Ny){
    double a;

    a = ( ((P[i*Ny+j]-P[(i-1)*Ny+j])/hx(i-1)-(P[(i+1)*Ny+j]-P[i*Ny+j])/hx(i)) / (0.5*(hx(i)+hx(i-1))) + 
			((P[i*Ny+j]-P[i*Ny+(j-1)])/hy(j-1)-(P[i*Ny+j+1]-P[i*Ny+j])/hy(j)) / (0.5*(hy(j)+hy(j-1))) );

    return -a;
}

void ComputeDim (int size, int& dim_x, int& dim_y){
	if (size >= 512) {
        dim_x = 16;
        dim_y = 32;
	} else if (size >= 256) {
       dim_x = 16;
       dim_y = 16;
	} else if (size >= 128) {
          dim_x = 8;
          dim_y = 16;
 	 } else if (size >= 64) {
			dim_x = 8;
			dim_y = 8;
		} 
		else if (size >= 32) {
				dim_x = 4;
				dim_y = 8;
			} else if (size >= 16) {
					dim_x = 4;
					dim_y = 4;
				} else if (size >= 8){
						dim_x = 2;
						dim_y = 4;
					} else if (size >= 4) {
							dim_x = 2;
							dim_y = 2;
						} else if (size >= 2) {
								dim_x = 1;
								dim_y = 2;
							} else if (size >= 1) {
									dim_x = 1;
									dim_y = 1;
								} // else { throw Exception("Wrong process amount given for computations");
}

void SplitFunction(int N, int p0, int p1, int& n0, int& n1){
//must split only interior nodes => their number is N-1 (A1=x_0 < x_1 <..< x_N=A2: x_0 and x_N - boundary)
	n0 = (N-1) / p0; 
	n1 = (N-1) / p1;
	return;
}

int GridEvenGen(int XAll, int YAll, int XBeg, int XCount, int YBeg, int YCount, double* Xnodes, double* Ynodes){
	// generetion of part of even grig (XCount)*(YCount)
	for(int i=0; i<XCount; i++){
		Xnodes[i] = A*(i+XBeg)/(1.0*XAll);
	}
	for(int j=0; j<YCount; j++)
		Ynodes[j] = B(j+YBeg)/(1.0*YAll);
	return 0;
}

void SendColumn(double *matr, int nx, int ny, int ncol, int pto, int tag, MPI_Comm comm){
//function for sending column from current processor to pto
//MPI_Send(&matr[ncol], 1, Column, pto, tag, comm);
//	MPI_Request request;
	int kx;
	for (kx=0; kx<nx; kx++)
		MPI_Send(&matr[kx*ny+ncol], 1, MPI_DOUBLE, pto, tag, comm);
		//MPI_Isend(&matr[kx*ny+ncol], 1, MPI_DOUBLE, pto, tag, comm, &request);
	return;
}

void RecvColumn(double *matr, int nx, int ny, int ncol, int pfrom, int tag, MPI_Comm comm){
//function for sending column from current processor to pto
//MPI_Recv(&matr[ncol], 1, Column, pfrom, tag, comm, status);
	MPI_Status status;
	int kx;
	for (kx=0; kx<nx; kx++)
		MPI_Recv(&matr[kx*ny+ncol], 1, MPI_DOUBLE, pfrom, tag, comm, &status);
	return;
}

void SendRow(double *matr, int nx, int ny, int nrow, int pto, int tag, MPI_Comm comm){
//function for sending column from current processor to pto
//MPI_Send(&matr[nrow], 1, Row, pto, tag, comm);
//	MPI_Request request;
	int ky;
	for (ky=0; ky<ny; ky++)
		MPI_Send(&matr[nrow*ny+ky], 1, MPI_DOUBLE, pto, tag, comm);
		//MPI_Isend(&matr[nrow*ny+ky], 1, MPI_DOUBLE, pto, tag, comm, &request);
	return;
}

void RecvRow(double *matr, int nx, int ny, int nrow, int pfrom, int tag, MPI_Comm comm){
//function for sending column from current processor to pto
//MPI_Recv(&matr[nrow], 1, Column, pfrom, tag, comm, status);
	MPI_Status status;
	int ky;
	for (ky=0; ky<ny; ky++)
		MPI_Recv(&matr[nrow*ny+ky], 1, MPI_DOUBLE, pfrom, tag, comm, &status);
	return;
}


int main(int argc, char **argv)
{
	// MPI initialization 
	int rank, size;
	int i, j;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;
	double time_e, time_b;
	time_b = MPI_Wtime();

	// Creation of topology
	MPI_Comm comm_2D;
	int const ndim = 2;       // Topology dimention
	int dims[ndim];           
	int periods[ndim]={0,0};
	int coords[ndim];         // Current process coordinates in topology
	int p0, p1;		// dimention of topology 
	// How to split nods among processors
	int N = 1000;   // N - number of nods
	int n0, n1;		// number of nodes on the current processor along 2 dimentions  
	int k0, k1;		// modulo from dividing N over n0 and n1
	int NX, NY;		// number of inner points in the current process

	ComputeDim(size, p0, p1);
	SplitFunction (N, p0, p1, n0, n1);
	k0 = (N-1) - p0*n0; k1 = (N-1) - p1*n1;
	NX = n0; NY = n1;
	dims[0] = p0, dims[1] = p1;

	MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, 0, &comm_2D);  // Topology creation
	MPI_Cart_coords(comm_2D, rank, ndim, coords);                       // Identify coordinates for current
	int XBeg, YBeg;
	if (coords[0] < k0) { NX++; XBeg = 1+coords[0]*NX;}
	else XBeg = 1 + k0*(NX+1) + (coords[0]-k0)*NX;
    if (coords[1] < k1) { NY++; YBeg = 1+coords[1]*NY;}
	else YBeg = 1 + k0*(NY+1) + (coords[1]-k1)*NY;

	int NX_p = NX+2; int NY_p = NY+2;
	XNodes = (double *)malloc((NX_p)*sizeof(double));
	YNodes = (double *)malloc((NY_p)*sizeof(double));
	double *CurVec = (double *)malloc((NX_p)*(NY_p)*sizeof(double)); // Zero approximation
	double *RVec = (double *)malloc((NX_p)*(NY_p)*sizeof(double)); // Next step approximation
	double *FMatr = (double *)malloc((NX_p)*(NY_p)*sizeof(double));

	// fprintf(fp, "Process %d: [0,%f]x[0,%f], number of points: N[0,A] = %d, N[0,B] = %d, Time=%f;\n", rank, A, B, NX, NY, time_b);

	//-------------------------------------------------------------------------------------

	// Derivative datatypes
	MPI_Datatype Column, Row;
	MPI_Type_vector(NX_p, 1, NY_p, MPI_DOUBLE, &Column);
	MPI_Type_commit(&Column);
	MPI_Type_vector(1, NY_p, NY_p, MPI_DOUBLE, &Row);
	MPI_Type_commit(&Row);

	// Neighbours
	int left, right, lower, upper;
	MPI_Cart_shift(comm_2D, 1, 1, &left, &right); //shift along dim[1] - along colomns
	MPI_Cart_shift(comm_2D, 0, 1, &upper, &lower); //shift along dim[0] - along rows

	// fprintf(fp, "Rank = %d, Coord[0] = %d, Coord[1] = %d, Left = %d, Right = %d, Lower = %d, Upper = %d\n", 
											//rank, coords[0], coords[1], left, right, lower, upper);
	// Initialization of Arrays
	// fprintf(fp, "Rank = %d, XBeg = %d, YBeg = %d, NX = %d, NY = %d\n", 
											//rank, XBeg, YBeg, NX, NY);
	GridEvenGen(N, N, XBeg-1, NX_p, YBeg-1, NY_p, XNodes, YNodes);

// --- Initial conditions ---
	for (i=0; i < NX_p; i++)
        for (j=0; j < NY_p; j++) {
            CurVec[i*NY_p+j] = 3.0;
            RVec[i*NY_p+j] = 0.0;
        }

    if (lower < 0)
		for (j=0; j<NY_p; j++)
			CurVec[(NX_p-1)*NY_p+j] = Borders(A, YNodes[j]);

	if (upper < 0)
		for (j=0; j<NY_p; j++)
			CurVec[j] = Borders(0.0, YNodes[j]);

	if (left < 0)
		for (i=0; i<NX_p; i++)
			CurVec[i*NY_p] = Borders(XNodes[i], 0.0);
	
	if (right < 0)
		for (i=0; i<NX_p; i++)
			CurVec[i*NY_p+(NY_p-1)] = Borders(XNodes[i], B);
// --- End of Initial conditions ---


// --- First step ---

	// printf("First step begins %d\n", rank);
	FMatrCalc(FMatr, NX_p, NY_p);

	// r(k) = Dx(k)-F(x_i, y_j)
	for(i=1; i<NX_p-1; i++)
		for (j=1; j<NY_p-1; j++){
			RVec[i*NY_p+j] = Delta(CurVec, i, j, NX_p, NY_p)-FMatr[i*NY_p+j];
			// printf("%d\n", FMatr[NX_p*i+j]);
		}
	
	// (r(k),r(k)) 
	double sp = 0.0;
	double sp_p = 0.0;
	#pragma omp parallel for reduction (+ : sp_p)
	for(i=1; i<NX_p-1; i++)
		for(j=1; j<NY_p-1; j++)
			sp_p += SQR(RVec[i*NY_p+j])*(0.5*(hx(i)+hx(i-1)))*(0.5*(hy(j)+hy(j-1)));
	double tau_p = sp_p;
	double tau;
	MPI_Allreduce(&tau_p, &tau, 1, MPI_DOUBLE, MPI_SUM, comm_2D);

	// sp = (Dr(k),r(k)) 
	sp = 0.0;
	sp_p = 0.0;
	#pragma omp parallel for reduction (+ : sp_p)
	for(i=1; i<NX_p-1; i++)
		for(j=1; j<NY_p-1; j++)
			sp_p += Delta(RVec,i,j,NX_p,NY_p)*RVec[i*NY_p+j]*(0.5*(hx(i)+hx(i-1)))*(0.5*(hy(j)+hy(j-1)));
    MPI_Allreduce(&sp_p, &sp, 1, MPI_DOUBLE, MPI_SUM, comm_2D);
	tau = tau/sp;

	// The p(k+1) = p(k) - tau*r(k) + count error 
	double err = 0.0;
	double err_p = 0.0;
	double NewValue;
	for (i=1; i<NX_p-1; i++)
		for(j=1; j<NY_p-1; j++){
			NewValue = CurVec[i*NY_p+j]-tau*RVec[i*NY_p+j];
			err_p = Max(err_p, fabs(Solution(XNodes[i], YNodes[j])-NewValue));
			CurVec[i*NY_p+j] = NewValue;
		}
    MPI_Allreduce(&err_p, &err, 1, MPI_DOUBLE, MPI_MAX, comm_2D);
	
// --- End First step ---
// --- Iterations ---
	//printf("Iterations begin\n");

	double *GVec = RVec;    // g(k) = r(k-1).
	RVec = (double *)malloc(NX_p*NY_p*sizeof(double));
	memset(RVec, 0, NX_p*NY_p*sizeof(double));

    int tag = 10;
	int iter = 0;

	while((err > 0.0001)){

		// exchange of border values
		//MPI_Send(&CurVec[1], 1, Column, left, tag, comm_2D);
		if (left >= 0)  SendColumn(CurVec, NX_p, NY_p, 1, left, tag, comm_2D);
		//MPI_Recv(&CurVec[NY_p-1], 1, Column, right, tag, comm_2D, &status);
		if (right >= 0) RecvColumn(CurVec, NX_p, NY_p, NY_p-1, right, tag, comm_2D);
		
		//MPI_Send(&CurVec[NY], 1, Column, right, tag, comm_2D);
		if (right >= 0) SendColumn(CurVec, NX_p, NY_p, NY, right, tag, comm_2D);
		//MPI_Recv(&CurVec[0], 1, Column, left, tag, comm_2D, &status);
		if (left >= 0)  RecvColumn(CurVec, NX_p, NY_p, 0, left, tag, comm_2D);

		//MPI_Send(&CurVec[NX*NY_p], 1, Row, lower, tag, comm_2D);
		if (lower >= 0) SendRow(CurVec, NX_p, NY_p, NX, lower, tag, comm_2D);
		//MPI_Recv(&CurVec[0], 1, Row, upper, tag, comm_2D, &status);
		if (upper >= 0) RecvRow(CurVec, NX_p, NY_p, 0, upper, tag, comm_2D);

		//MPI_Send(&CurVec[1*NY_p], 1, Row, upper, tag, comm_2D);
		if (upper >= 0) SendRow(CurVec, NX_p, NY_p, 1, upper, tag, comm_2D);
		//MPI_Recv(&CurVec[(NX_p-1)*NY_p], 1, Row, lower, tag, comm_2D, &status);
		if (lower >= 0) RecvRow(CurVec, NX_p, NY_p, NX_p-1, lower, tag, comm_2D);

		iter++;

		// r(k) = Dp(k)-F(x_i, y_j)
		for(i=1; i<NX_p-1; i++)
			for(j=1; j<NY_p-1; j++)
				RVec[i*NY_p+j] = Delta(CurVec,i,j,NX_p,NY_p)-FMatr[i*NY_p+j];

		// (Dr(k),g(k-1))
		double alpha = 0.0;
		double alpha_p = 0.0;
		#pragma omp parallel for reduction (+ : alpha_p)
		for(i=1; i<NX_p-1; i++)
			for(j=1; j<NY_p-1; j++)
				alpha_p += Delta(RVec,i,j,NX_p,NY_p)*GVec[i*NY_p+j]*(0.5*(hx(i)+hx(i-1)))*(0.5*(hy(j)+hy(j-1)));
		MPI_Allreduce(&alpha_p, &alpha, 1, MPI_DOUBLE, MPI_SUM, comm_2D);

		// sp = (Ag(k-1),g(k-1))
		sp = 0.0;
		sp_p = 0.0;
		#pragma omp parallel for reduction (+ : sp_p)
		for(i=1; i<NX_p-1; i++)
			for(j=1; j<NY_p-1; j++)
				sp_p += Delta(GVec,i,j,NX_p,NY_p)*GVec[i*NY_p+j]*(0.5*(hx(i)+hx(i-1)))*(0.5*(hy(j)+hy(j-1)));
		MPI_Allreduce(&sp_p, &sp, 1, MPI_DOUBLE, MPI_SUM, comm_2D);
		alpha = alpha/sp;

		// g(k)= r(k)-alpha*g(k-1)
		for(i=1; i<NX_p-1; i++)
			for(j=1; j<NY_p-1; j++)
				GVec[i*NY_p+j] = RVec[i*NY_p+j]-alpha*GVec[i*NY_p+j];

		// tau_p = (r(k),g(k))
		tau = 0.0;
		tau_p = 0.0;
		#pragma omp parallel for reduction (+ : tau_p)
		for(i=1; i<NX_p-1; i++)
			for(j=1; j<NY_p-1; j++)
				tau_p += RVec[i*NY_p+j]*GVec[i*NY_p+j]*(0.5*(hx(i)+hx(i-1)))*(0.5*(hy(j)+hy(j-1)));
        MPI_Allreduce(&tau_p, &tau, 1, MPI_DOUBLE, MPI_SUM, comm_2D);

		// sp = (Dg(k),g(k)) 
		sp = 0.0;
		sp_p = 0.0;
		#pragma omp parallel for reduction (+ : sp_p)
		for(i=1; i<NX_p-1; i++)
			for(j=1; j<NY_p-1; j++)
				sp_p += Delta(GVec,i,j,NX_p,NY_p)*GVec[i*NY_p+j]*(0.5*(hx(i)+hx(i-1)))*(0.5*(hy(j)+hy(j-1)));
		MPI_Allreduce(&sp_p, &sp, 1, MPI_DOUBLE, MPI_SUM, comm_2D);
		tau = tau/sp;

		// p(k+1) = p(k) - tau*g(k) 
		err = 0.0;
		err_p = 0.0;
		for(i=1; i<NX_p-1; i++)
			for(j=1; j<NY_p-1; j++){
				NewValue = CurVec[i*NY_p+j]-tau*GVec[i*NY_p+j];
				err_p = Max(err_p, fabs(CurVec[i*NY_p + j] - NewValue)); 
				CurVec[i*NY_p+j] = NewValue;
			}
        MPI_Allreduce(&err_p, &err, 1, MPI_DOUBLE, MPI_MAX, comm_2D);
		//if ((rank == 0) && (iter%500==1)) // fprintf(fp,"%d :: , it = %d, error is %f\n", rank, iter, err);
	}
// --- End Iterations ---
	// if (rank == 0) // fprintf(fp,"%d :: , it = %d, error is %f\n", rank, iter, err);

	err = 0.0;
    err_p = 0.0;
	for(i=1; i<NX_p-1; i++)
		for(j=1; j<NY_p-1; j++)
			err_p = Max(err_p, fabs(Solution(XNodes[i], YNodes[j]) - CurVec[i*NY_p + j]));
    MPI_Allreduce(&err_p, &err, 1, MPI_DOUBLE, MPI_MAX, comm_2D);

	//if (rank == 0) // fprintf(fp,"%d :: The iterations ended. The residual error is %f, iterations = %d\n", rank, err, iter);
	if (rank == 0) printf("%d :: The iterations ended. The residual error is %f, iterations = %d\n", rank, err, iter);

// print solution to file
/*	FILE *fp_ij; FILE *fp_p;
	fp_ij = fopen("pixc.txt", "a");
	fp_p = fopen("pixp.txt", "a");
	// fprintf(fp_ij,"# index x_i y_j p_ij\n");
	// fprintf(fp_p,"# index x_i y_j pres\n");
		for(i=0; i<NX_p; i++){
			for(j=0; j<NY_p; j++){
				// fprintf(fp_ij,"%f %f %f\n", XNodes[i], YNodes[j], CurVec[i*NY_p+j]);
				// fprintf(fp_p,"%f %f %f\n", XNodes[i], YNodes[j], Solution(XNodes[i],YNodes[j]));
			}
			// fprintf(fp_ij, "\n");
			// fprintf(fp_p, "\n");
		}
	fclose(fp_ij);
	fclose(fp_p);
*/	
	free(CurVec);  
	free(FMatr);
	free(GVec); 
	free(RVec);
	free(XNodes); free(YNodes);

	time_e = MPI_Wtime() - time_b;

	//if (rank == 0) // fprintf(fp,"%d :: The time = %f\n", rank, time_e);
	if (rank == 0) printf("%d :: Time = %f\n", rank, time_e);
	
	//fclose(fp);
	MPI_Finalize();
	return 0;
}
