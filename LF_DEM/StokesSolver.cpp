#include "StokesSolver.h"
#ifdef TRILINOS
#include <BelosCGIteration.hpp>
#endif
using namespace std;
#define DELETE(x) if(x){delete [] x; x = NULL;}

/******************************************************
 *                                                     *
 *                   Public Methods                    *
 *                                                     *
 ******************************************************/

StokesSolver::~StokesSolver(){
    if (!dblocks) {
		delete [] dblocks;
	}
	if (!odblocks) {
		delete [] odblocks;
	}

	if(!chol_solution) {
		cholmod_free_dense(&chol_solution, &chol_c);
	}
	if(!chol_rhs) {
		cholmod_free_dense(&chol_rhs, &chol_c);
	}
	if (brownian) {
		cholmod_free_dense(&chol_brownian_rhs, &chol_c);
	}
	if(!chol_res_matrix) {
		cholmod_free_sparse(&chol_res_matrix, &chol_c);
	}
	if(chol_init) {
		cholmod_finish(&chol_c);
	}
	
#ifdef TRILINOS
	for (int i=0; i<res_matrix_linear_size; i++) {
		delete [] columns[i];
	}
	delete [] columns;
	delete [] columns_nb;
	for (int i=0; i<res_matrix_linear_size; i++){
		delete [] values[i];
	}
	delete [] values;
#endif
}

void
StokesSolver::init(int n, bool is_brownian){
	np = n;
    np3 = 3*np;
	brownian = is_brownian;
	// initializing values that can be changed later
	_direct = true;
	_iterative = false;
	chol_init = false;
	FTcoupling = false;
}

void
StokesSolver::initialize(){

	// CHOLMOD parameters
	stype = -1; // 1 is symmetric, stored upper triangular (UT), -1 is LT
	sorted = 0;		/* TRUE if columns sorted, FALSE otherwise*/
	packed = 1;		/* TRUE if matrix packed, FALSE otherwise */
	xtype = CHOLMOD_REAL;

	
#ifdef TRILINOS
	// TRILINOS init and parameters
	// initialize solver
	RCP<ParameterList> solverParams = parameterList();
	// parameters to be tuned (and understood!)
	int blocksize = 10;
	int maxiters = 400;
	double tol = 1.e-6;
	solverParams->set("Block Size", blocksize);              // Blocksize to be used by iterative solver
	solverParams->set("Maximum Iterations", maxiters);       // Maximum number of iterations allowed
	solverParams->set("Convergence Tolerance", tol);         // Relative convergence tolerance requested
	solverParams->set("Verbosity", Belos::Errors + Belos::Warnings);
	tril_solver = tril_factory.create ("CG", solverParams);
	// initialize empty linear problem
    tril_stokes_equation = rcp(new Belos::LinearProblem <SCAL, VEC, MAT> ());
#endif


	// resistance matrix characteristics (see header for matrix description)
	if(FTcoupling){
		res_matrix_linear_size = 6*np;
		dblocks_nb = 2*np;
		sdblocks_nb = 2*np;
		dblocks_element_nb = 12*np;
	}
	else{
		res_matrix_linear_size = 3*np;
		dblocks_nb = np;
		sdblocks_nb = 0;
		dblocks_element_nb = 6*np;
	}



	allocateRessources();
	chol_L_to_be_freed = false;
}

/************* Matrix filling methods **********************/

// Diagonal Terms, FT/UW version
void
StokesSolver::addToDiag(int ii, double FUvalue, double TWvalue){
	if (direct()) {
		int ii18 = 18*ii;
		dblocks[ii18   ] += FUvalue;
		dblocks[ii18+6 ] += FUvalue;
		dblocks[ii18+10] += FUvalue;
		dblocks[ii18+12] += TWvalue;
		dblocks[ii18+15] += TWvalue;
		dblocks[ii18+17] += TWvalue;
	}
#ifdef TRILINOS
	if (iterative()) {

		cerr << " Error : StokesSolver::addToDiag(const vec3d &nvec, int ii, double FUvalue, double TWvalue) not implemented for TRILINOS yet ! " << endl;
		exit(1);

		int ii3 = 3*ii;
		for (int j = 0; j < 3; j ++) {
			values[ii3+j][j] += FUvalue;
		}
	}
#endif
}


// Diagonal Blocks Terms, FT/UW version
void
StokesSolver::addToDiagBlock(const vec3d &nvec, int ii, double XA, double YB, double YC){
	
    double XA_n0 = XA*nvec.x;
    double XA_n1 = XA*nvec.y;
    double XA_n2 = XA*nvec.z;
    double XA_n1n0 = XA_n0*nvec.y;
    double XA_n2n1 = XA_n1*nvec.z;
    double XA_n0n2 = XA_n2*nvec.x;
	if (direct()) {
		int ii18 = 18*ii;
		dblocks[ii18   ] += XA_n0*nvec.x;   // 00 element of the dblock
		dblocks[ii18+1 ] += XA_n1n0;        // 10
		dblocks[ii18+2 ] += XA_n0n2;        // 20
		dblocks[ii18+6 ] += XA_n1*nvec.y;   // 11
		dblocks[ii18+7 ] += XA_n2n1;        // 21
		dblocks[ii18+10] += XA_n2*nvec.z;   // 22
	}
#ifdef TRILINOS
	if (iterative()) {
		int iidof = dof*ii;
		cerr << " Error : StokesSolver::addToDiagBlock(const vec3d &nvec, int ii, double XA, double YB, double YC) not implemented for TRILINOS yet ! " << endl;
		exit(1);
		values[iidof  ][0] += XA_n0*nvec.x; // 00
		values[iidof  ][1] += XA_n1n0; // 01
		values[iidof  ][2] += XA_n0n2; // 02
		values[iidof+1][0] += XA_n1n0; // 10
		values[iidof+1][1] += XA_n1*nvec.y; // 11
		values[iidof+1][2] += XA_n2n1; // 12
		values[iidof+2][0] += XA_n0n2; // 20
		values[iidof+2][1] += XA_n2n1; // 21
		values[iidof+2][2] += XA_n2*nvec.z; // 20
	}
#endif
}


// Off-Diagonal Blocks Terms, FT/UW version
void
StokesSolver::setOffDiagBlock(const vec3d &nvec, int ii, int jj, double XA,double YB, double YC){
	if (direct()) {
		setColumn(nvec, ii, jj, XA, YB, YC);
	}
#ifdef TRILINOS
	if (iterative()) {
		setRow(nvec, ii, jj, XA, YB, YC);
	}
#endif
}


/*************** Cholmod Matrix Filling *************
 Cholmod matrices we are using are defined in column major order (index j is column index)
 
 Cholmod matrices are defined as follows:
 - all values are stored in array x ( size nzmax )
 - locations of values are encoded in array p ( size np ):
 values corresponding to column j are x[ p[j] ]  to x[ p[j+1] - 1 ]
 - corresponding rows are stored in array i ( size nzmax ):
 rows corresponding to column j are i[ p[j] ]  to i[ p[j+1] - 1 ]
 
 Hence:
 with p[j] < a < p[j+1]-1
        . . . . j . . . . . .
     .|         .            |
     .|         .            |
     .|         .            |
 i[a] | . . . .x[a]          |
     .|                      |
     .|                      |
 
 
 *****************************************************/
void
StokesSolver::completeResistanceMatrix_cholmod(){
	// this function is commented, but you are strongly advised to read 
	// the description of storage in the header file first :)

    // declare the last 2 values of odbrows_table
    odbrows_table[np-1] = (unsigned int)odbrows.size();
    odbrows_table[np] = (unsigned int)odbrows.size();
    allocateResistanceMatrix();
    // fill
	for (int j=0; j<np; j++) {
		
		// associated with particle j are 6 columns in the matrix:
		// { 6j, ... , 6j+5 }
		int j6 = 6*j;
		int j18 = 18*j;
		int j21 = 21*j;
		
		// the number of non-zero elements before column 6j is:
		// - 21*j from j diagonal blocks
		// - 21*odbrows_table[j] from odbrows_table[j] off-diagonal blocks
		//
		// the number of non-zero elements before column 6j+1 is:
		// - number of non-zero before column 6j + number of non-zero in column 6*j
		// (in 6j: 6 elements in diagonal block, plus 6*(odbrows_table[j+1]-odbrows_table[j])
		//
		// for 6j+2 --> 6j+5: same idea
			
		((int*)chol_res_matrix->p)[j6  ] = j21   + 21*odbrows_table[j];
		((int*)chol_res_matrix->p)[j6+1] = ((int*)chol_res_matrix->p)[j6] + 6 + 6*(odbrows_table[j+1]-odbrows_table[j]);
		((int*)chol_res_matrix->p)[j6+2] = ((int*)chol_res_matrix->p)[j6+1] + 5 + 6*(odbrows_table[j+2]-odbrows_table[j+1]);
		((int*)chol_res_matrix->p)[j6+3] = ((int*)chol_res_matrix->p)[j6+2] + 4 + 6*(odbrows_table[j+3]-odbrows_table[j+2]);
		((int*)chol_res_matrix->p)[j6+4] = ((int*)chol_res_matrix->p)[j6+3] + 3 + 6*(odbrows_table[j+4]-odbrows_table[j+3]);
		((int*)chol_res_matrix->p)[j6+5] = ((int*)chol_res_matrix->p)[j6+4] + 2 + 6*(odbrows_table[j+5]-odbrows_table[j+4]);	

		int pj6   = ((int*)chol_res_matrix->p)[j6];
		int pj6_1 = ((int*)chol_res_matrix->p)[j6+1];
		int pj6_2 = ((int*)chol_res_matrix->p)[j6+2];
		int pj6_3 = ((int*)chol_res_matrix->p)[j6+3];
		int pj6_4 = ((int*)chol_res_matrix->p)[j6+4];
		int pj6_5 = ((int*)chol_res_matrix->p)[j6+5];
			
		// diagonal block row indices (21)
		((int*)chol_res_matrix->i)[pj6  ] = j6;   // column j6
		((int*)chol_res_matrix->i)[pj6+1] = j6+1;
		((int*)chol_res_matrix->i)[pj6+2] = j6+2;
		((int*)chol_res_matrix->i)[pj6+3] = j6+3;
		((int*)chol_res_matrix->i)[pj6+4] = j6+4;
		((int*)chol_res_matrix->i)[pj6+5] = j6+5;
		((int*)chol_res_matrix->i)[pj6_1  ] = j6+1;    // column j6+1
		((int*)chol_res_matrix->i)[pj6_1+1] = j6+2;
		((int*)chol_res_matrix->i)[pj6_1+2] = j6+3;
		((int*)chol_res_matrix->i)[pj6_1+3] = j6+4;
		((int*)chol_res_matrix->i)[pj6_1+4] = j6+5;
		((int*)chol_res_matrix->i)[pj6_2  ] = j6+2;    // column j6+2
		((int*)chol_res_matrix->i)[pj6_2+1] = j6+3;
		((int*)chol_res_matrix->i)[pj6_2+2] = j6+4;
		((int*)chol_res_matrix->i)[pj6_2+3] = j6+5;
		((int*)chol_res_matrix->i)[pj6_3  ] = j6+3;    // column j6+3
		((int*)chol_res_matrix->i)[pj6_3+1] = j6+4;
		((int*)chol_res_matrix->i)[pj6_3+2] = j6+5;
		((int*)chol_res_matrix->i)[pj6_4  ] = j6+4;    // column j6+4
		((int*)chol_res_matrix->i)[pj6_4+1] = j6+5;
		((int*)chol_res_matrix->i)[pj6_5  ] = j6+5;    // column j6+5

		// diagonal blocks row values (21)
		((double*)chol_res_matrix->x)[pj6  ] = dblocks[j18];   // column j6
		((double*)chol_res_matrix->x)[pj6+1] = dblocks[j18+1];
		((double*)chol_res_matrix->x)[pj6+2] = dblocks[j18+2];
		((double*)chol_res_matrix->x)[pj6+3] = dblocks[j18+3];
		((double*)chol_res_matrix->x)[pj6+4] = dblocks[j18+4];
		((double*)chol_res_matrix->x)[pj6+5] = dblocks[j18+5];
		((double*)chol_res_matrix->x)[pj6_1  ] = dblocks[j18+6];   // column j6+1
		((double*)chol_res_matrix->x)[pj6_1+1] = dblocks[j18+7];
		((double*)chol_res_matrix->x)[pj6_1+2] = -dblocks[j18+3];   // anti-symmetry
		((double*)chol_res_matrix->x)[pj6_1+3] = dblocks[j18+8];
		((double*)chol_res_matrix->x)[pj6_1+4] = dblocks[j18+9];
		((double*)chol_res_matrix->x)[pj6_2  ] = dblocks[j18+10];   // column j6+2
		((double*)chol_res_matrix->x)[pj6_2+1] = -dblocks[j18+5];   // anti-symmetry
		((double*)chol_res_matrix->x)[pj6_2+2] = -dblocks[j18+9];   // anti-symmetry
		((double*)chol_res_matrix->x)[pj6_2+3] = dblocks[j18+11];
		((double*)chol_res_matrix->x)[pj6_3  ] = dblocks[j18+12];   // column j6+3
		((double*)chol_res_matrix->x)[pj6_3+1] = dblocks[j18+13];
		((double*)chol_res_matrix->x)[pj6_3+2] = dblocks[j18+14];
		((double*)chol_res_matrix->x)[pj6_4  ] = dblocks[j18+15];   // column j6+4
		((double*)chol_res_matrix->x)[pj6_4+1] = dblocks[j18+16];
		((double*)chol_res_matrix->x)[pj6_5  ] = dblocks[j18+17];   // column j6+5
			
		/*****  2  : off-diagonal blocks row indices and values ***********/
		// 36 non-zero elements per block

		for(int k = odbrows_table[j]; k < odbrows_table[j+1]; k++){
			int u = 6*(k-odbrows_table[j]); 

			// we are filling the "k-odbFrows_table[j]"th off-diag block of the column.
			// For column j6, for exemple, the indices of the non-zero values are:
			// pj6 for all non-zero elements before column j6,
			// + 6 for the diagonal block of column j6
			// + u (=6*(k-odbFrows_table[j])) for the off-diag blocks of j6
			// + index inside the current block
			for(int s=0; s<6;s++){
				((int*)chol_res_matrix->i)[pj6+6 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_1+5 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_2+4 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_3+3 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_4+2 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_5+1 + u +s] = odbrows[k]+s;
			}

			int k6 = 6*k;
			int k4 = 4*k;
			int k3 = 3*k;
			int k2 = 2*k;

			((double*)chol_res_matrix->x)[pj6+6 + u   ]   = odblocks[0][k6  ];    // column j6
			((double*)chol_res_matrix->x)[pj6+6 + u +1]   = odblocks[0][k6+1];
			((double*)chol_res_matrix->x)[pj6+6 + u +2]   = odblocks[0][k6+2];
			((double*)chol_res_matrix->x)[pj6+6 + u +3]   = odblocks[0][k6+3];
			((double*)chol_res_matrix->x)[pj6+6 + u +4]   = odblocks[0][k6+4];
			((double*)chol_res_matrix->x)[pj6+6 + u +5]   = odblocks[0][k6+5];

			((double*)chol_res_matrix->x)[pj6_1+5 + u   ] = odblocks[0][k6+1];  // symmetry  // column j6+1
			((double*)chol_res_matrix->x)[pj6_1+5 + u +1] = odblocks[1][k4  ];
			((double*)chol_res_matrix->x)[pj6_1+5 + u +2] = odblocks[1][k4+1];
			((double*)chol_res_matrix->x)[pj6_1+5 + u +3] = -odblocks[0][k6+4];  // antisymmetry
			((double*)chol_res_matrix->x)[pj6_1+5 + u +4] = odblocks[1][k4+2];
			((double*)chol_res_matrix->x)[pj6_1+5 + u +5] = odblocks[1][k4+3];

			((double*)chol_res_matrix->x)[pj6_2+4 + u   ] = odblocks[0][k6+2];    // symmetry  // column j6+2
			((double*)chol_res_matrix->x)[pj6_2+4 + u +1] = odblocks[1][k4+1];    // symmetry
			((double*)chol_res_matrix->x)[pj6_2+4 + u +2] = odblocks[2][k2  ];
			((double*)chol_res_matrix->x)[pj6_2+4 + u +3] = -odblocks[0][k6+5];  // antisymmetry
			((double*)chol_res_matrix->x)[pj6_2+4 + u +4] = -odblocks[1][k4+3];  // antisymmetry
			((double*)chol_res_matrix->x)[pj6_2+4 + u +5] = odblocks[2][k2+1];
			
			((double*)chol_res_matrix->x)[pj6_3+3 + u   ] = odblocks[0][k6+3];    // symmetry // column j6+3
			((double*)chol_res_matrix->x)[pj6_3+3 + u +1] = -odblocks[0][k6+4];    // symmetry+antisymmetry
			((double*)chol_res_matrix->x)[pj6_3+3 + u +2] = -odblocks[0][k6+5];    // symmetry+antisymmetry
			((double*)chol_res_matrix->x)[pj6_3+3 + u +3] = odblocks[3][k3  ];
			((double*)chol_res_matrix->x)[pj6_3+3 + u +4] = odblocks[3][k3+1];
			((double*)chol_res_matrix->x)[pj6_3+3 + u +5] = odblocks[3][k3+2];
			
			((double*)chol_res_matrix->x)[pj6_4+2 + u   ] = odblocks[0][k6+4];  // symmetry  // column j6+4
			((double*)chol_res_matrix->x)[pj6_4+2 + u +1] = odblocks[1][k4+2]; // symmetry
			((double*)chol_res_matrix->x)[pj6_4+2 + u +2] = -odblocks[1][k4+3]; // symmetry+antisymmetry
			((double*)chol_res_matrix->x)[pj6_4+2 + u +3] = odblocks[3][k3+1]; // symmetry
			((double*)chol_res_matrix->x)[pj6_4+2 + u +4] = odblocks[4][k2  ];
			((double*)chol_res_matrix->x)[pj6_4+2 + u +5] = odblocks[4][k2+1];
			
			((double*)chol_res_matrix->x)[pj6_5+1 + u   ] = odblocks[0][k6+5];  // symmetry  // column j6+5
			((double*)chol_res_matrix->x)[pj6_5+1 + u +1] = odblocks[1][k4+3]; // symmetry
			((double*)chol_res_matrix->x)[pj6_5+1 + u +2] = odblocks[2][k2+1]; // symmetry
			((double*)chol_res_matrix->x)[pj6_5+1 + u +3] = odblocks[3][k3+2]; // symmetry
			((double*)chol_res_matrix->x)[pj6_5+1 + u +4] = odblocks[4][k2+1]; // symmetry
			((double*)chol_res_matrix->x)[pj6_5+1 + u +5] = odblocks[5][k   ];
			
		}
	}

    ((int*)chol_res_matrix->p)[np6] = ((int*)chol_res_matrix->p)[np6-1]+1;
	
	factorizeResistanceMatrix();
}


/*************** Epetra_CrsMatrix Filling *************
 Epetra_CrsMatrix we are using are defined in row major order.
 
 Epetra_CrsMatrix must be stored completely.
 
 Epetra_CrsMatrix elements are not accessed directly for filling.
 Instead we use user friendly methods, that take one row at a time.
 
 *****************************************************/

void
StokesSolver::completeResistanceMatrix_trilinos(){
#ifdef TRILINOS
    for (int i = 0; i < res_matrix_linear_size; i++) {
		tril_res_matrix->InsertGlobalValues(i, columns_nb[i] , values[i], columns[i]);
    }
    // FillComplete matrix before building the preconditioner
    tril_res_matrix->FillComplete();
	tril_stokes_equation->setOperator(rcp(tril_res_matrix, false));
	//setDiagBlockPreconditioner();
	//	setIncCholPreconditioner();
	//	setSpInvPreconditioner();
#endif
}

void
StokesSolver::completeResistanceMatrix(){
	if (direct()) {
		completeResistanceMatrix_cholmod();
	}
#ifdef TRILINOS
	if (iterative()) {
		completeResistanceMatrix_trilinos();
	}
#endif
}

void
StokesSolver::resetResistanceMatrix(string solver_type){
	setSolverType(solver_type);
	
	if (direct()) {
		for (int k=0; k<dblocks_element_nb; k++){
			dblocks[k] = 0;
		}
		odbrows.clear();
		odblocks[0].clear();
		odblocks[1].clear();
		odblocks[2].clear();
		odbrows_table[0] = 0;
	}
#ifdef TRILINOS
	if (iterative()) {
		tril_res_matrix = new Epetra_CrsMatrix(Copy, *Map, 20*dof+dof );
		tril_l_precond = new Epetra_CrsMatrix(Copy, *Map, 3);
		tril_res_matrix->PutScalar(0.);
		for (int i=0; i<res_matrix_linear_size; i++){
			for (int j=0; j<columns_max_nb; j++){
				columns[i][j] = -1;
				values[i][j] = 0;
			}
		}
		// declare the diagonal blocks
		for (int i=0; i<np; i++){
			int idof = dof*i;
			for (int j=0; j<dof; j++){
				columns[idof  ][j] = idof+j;
				columns[idof+1][j] = idof+j;
				columns[idof+2][j] = idof+j;
			}
			columns_nb[idof  ] = dof;
			columns_nb[idof+1] = dof;
			columns_nb[idof+2] = dof;
		}
	}
#endif
}

void
StokesSolver::resetRHS(){
	if (direct()) {
		for (int i=0; i<res_matrix_linear_size; i++){
			((double*)chol_rhs->x)[i] = 0;
		}
	}
#ifdef TRILINOS
	if (iterative()) {
		tril_rhs->PutScalar(0);
	}
#endif
}

void
StokesSolver::addToRHS(int i, double val){
	if (direct()) {
		((double*)chol_rhs->x)[i] += val;
	}
#ifdef TRILINOS
	if (iterative()) {
		tril_rhs->SumIntoGlobalValue(i, 0, val);
	}
#endif
}

void
StokesSolver::addToRHS(double *rhs){
	if (direct()) {
		for (int i=0; i<res_matrix_linear_size; i++) {
			((double*)chol_rhs->x)[i] += rhs[i];
		}
	}
	
#ifdef TRILINOS
	if (iterative()) {
		for (int i=0; i<res_matrix_linear_size; i++) {
			tril_rhs->SumIntoGlobalValue(i, 0, rhs[i]);
		}
	}
#endif
}

void
StokesSolver::setRHS(double* rhs){
	if (direct()) {
		for (int i=0; i<res_matrix_linear_size; i++) {
			((double*)chol_rhs->x)[i] = rhs[i];
		}
	}
	if (iterative()) {
		cerr << " Error : StokesSolver::setRHS(double* rhs) not implemented for TRILINOS yet ! " << endl;
		exit(1);
	}
}

void
StokesSolver::getRHS(double* rhs){
	if (direct()) {
		for (int i=0; i<res_matrix_linear_size; i++) {
			rhs[i] = ((double*)chol_rhs->x)[i];
		}
	}
#ifdef TRILINOS
	if (iterative()) {
		tril_rhs->ExtractCopy(rhs);
	}
#endif
}

void
StokesSolver::solve_CholTrans(double* velocity){
	if (direct()) {
		chol_PTsolution = cholmod_solve (CHOLMOD_Lt, chol_L, chol_rhs, &chol_c) ;
		chol_solution = cholmod_solve (CHOLMOD_Pt, chol_L, chol_PTsolution, &chol_c) ;
		for (int i=0; i<res_matrix_linear_size; i++) {
			velocity[i] = ((double*)chol_solution->x)[i];
		}
		cholmod_free_dense(&chol_solution, &chol_c);
		cholmod_free_dense(&chol_PTsolution, &chol_c);
	}
#ifdef TRILINOS
	if (iterative()) {
		cerr << " StokesSolver::solve_CholTrans(double* velocity) not implemented for iterative solver." << endl;
		exit(1);
	}
#endif
}

void
StokesSolver::solve(double* velocity){
	if (direct()) {
		chol_solution = cholmod_solve (CHOLMOD_A, chol_L, chol_rhs, &chol_c) ;
		for (int i=0; i<res_matrix_linear_size; i++) {
			velocity[i] = ((double*)chol_solution->x)[i];
		}
		cholmod_free_dense(&chol_solution, &chol_c);
	}
#ifdef TRILINOS
	if (iterative()) {
		tril_stokes_equation->setLHS(tril_solution);
		tril_stokes_equation->setRHS(tril_rhs);
		bool set_success = tril_stokes_equation->setProblem();
		if (!set_success) {
			cerr << "ERROR:  Belos::LinearProblem failed to set up correctly" << endl;
			exit(1);
		}
		tril_solver->setProblem (tril_stokes_equation);
		Belos::ReturnType ret = tril_solver->solve();
		if (ret != Belos::Converged) {
			cerr << " Warning: Belos::Solver did not converge" << endl;
		}
		tril_solver->getNumIters();
		//		int iter_steps = tril_solver->getNumIters();
		//		cout << " iterations " << iter_steps << endl;
		tril_solution->ExtractCopy(&velocity);
	}
#endif
}

void
StokesSolver::convertDirectToIterative(){
#ifdef TRILINOS
	// don't free the Cholesky factor, but rememver to do it when solvingIsDone
	cholmod_free_sparse(&chol_res_matrix, &chol_c);
	cholmod_free_dense(&chol_solution, &chol_c);
	chol_L_to_be_freed = true;
	// convert RHS
	for (int i=0; i<res_matrix_linear_size;i++) {
		tril_rhs->ReplaceGlobalValue(i, 0, ((double*)chol_rhs->x)[i]);
	}
	setSolverType("iterative");
#else
	cerr << " Error: StokesSolver::convertDirectToIterative() : no iterative solver. Compile withe Trilinos. " << endl;
	exit(1);
#endif
}

void
StokesSolver::solvingIsDone(){
	if (direct()) {
		cholmod_free_factor(&chol_L, &chol_c);
		cholmod_free_sparse(&chol_res_matrix, &chol_c);
		//	cholmod_free_dense(&chol_rhs, &chol_c);
	}
#ifdef TRILINOS
	if (iterative()) {
		delete tril_res_matrix;
		delete tril_l_precond;
		if (chol_L_to_be_freed) {
			cholmod_free_factor(&chol_L, &chol_c);
		}
	}
#endif
}

/******************************************************
 *                                                     *
 *                  Private Methods                    *
 *                                                     *
 ******************************************************/

void
StokesSolver::allocateRessources(){
#ifdef TRILINOS
    int maxnum_interactionpair_per_particle = 20;
    columns_max_nb = dof*maxnum_interactionpair_per_particle;
    int numlhs = 1;
    int numrhs = 1;
    Map = rcp(new Epetra_Map(res_matrix_linear_size, 0, Comm));
    tril_solution = rcp(new Epetra_Vector(*Map, numlhs));
    tril_rhs = rcp(new Epetra_Vector(*Map, numrhs));
    //	tril_res_matrix = rcp( new MAT(res_matrix_linear_size) );
    columns = new int* [res_matrix_linear_size];
    for (int i=0; i<res_matrix_linear_size; i++) {
		columns[i] = new int [columns_max_nb];
		for (int j=0; j<columns_max_nb; j++) {
			columns[i][j] = -1;
		}
    }
    values = new double* [res_matrix_linear_size];
    for (int i=0; i<res_matrix_linear_size; i++) {
		values[i] = new double [columns_max_nb];
		for (int j=0; j<columns_max_nb; j++) {
			values[i][j] = 0.;
		}
    }
    columns_nb = new int [res_matrix_linear_size];
    for (int i=0; i<res_matrix_linear_size; i++) {
		columns_nb[i] = 0;
	}
#endif
    cholmod_start (&chol_c);
	chol_init = true;
    dblocks = new double [6*np];
    odblocks = new vector <double> [3];
    odbrows_table = new int [np+1];
    chol_rhs = cholmod_allocate_dense(np3, 1, np3, xtype, &chol_c);
	for (int i=0; i<np3; i++) {
		((double*)chol_rhs->x)[i] = 0;
	}
    chol_L = NULL;
}

// only needed for Cholmod
void
StokesSolver::allocateResistanceMatrix(){
	// allocate
	int nzmax; // non-zero values
	nzmax = 6*np; // diagonal blocks
	for (int s=0; s<3; s++) {
		nzmax += odblocks[s].size();  // off-diagonal
	}
	chol_res_matrix = cholmod_allocate_sparse(np3, np3, nzmax, sorted, packed, stype,xtype, &chol_c);
}

void
StokesSolver::doneBlocks(int i){
	if (direct()) {
		odbrows_table[i+1] = (unsigned int)odbrows.size();
	}
}


// odblocks fillings, for FT/UW version
void
StokesSolver::setColumn(const vec3d &nvec, int ii, int jj, double XA, double YB, double YC){
    double XA_n0 = XA*nvec.x;
    double XA_n1 = XA*nvec.y;
    double XA_n2 = XA*nvec.z;
    double XA_n1n0 = XA_n0*nvec.y;
    double XA_n2n1 = XA_n1*nvec.z;
    double XA_n0n2 = XA_n2*nvec.x;


    double mYB_n0 = -YB*nvec.x;
    double mYB_n1 = -YB*nvec.y;
    double mYB_n2 = -YB*nvec.z;

    double YC_n0 = YC*nvec.x;
    double YC_n1 = YC*nvec.y;
    double YC_n2 = YC*nvec.z;
    double mYC_n1n0 = -YC_n0*nvec.y;
    double mYC_n2n1 = -YC_n1*nvec.z;
    double mYC_n0n2 = -YC_n2*nvec.x;
    double YC_1_m_n0n0 = YC*(1-nvec.x*nvec.x);
    double YC_1_m_n1n1 = YC*(1-nvec.y*nvec.y);
    double YC_1_m_n2n2 = YC*(1-nvec.y*nvec.y);


	odbrows.push_back(6*jj);

	odblocks[0].push_back(XA_n0*nvec.x); // 00
	odblocks[0].push_back(XA_n1n0); // 10
	odblocks[0].push_back(XA_n0n2); // 20
	odblocks[1].push_back(XA_n1n0); // 01
	odblocks[1].push_back(XA_n1*nvec.y); //11
	odblocks[1].push_back(XA_n2n1); // 21
	odblocks[2].push_back(XA_n0n2); // 02
	odblocks[2].push_back(XA_n2n1); // 12
	odblocks[2].push_back(XA_n2*nvec.z); //22

}


void
StokesSolver::setRow(const vec3d &nvec, int ii, int jj, double XA, double YB, double YC){
	cerr << " Error : StokesSolver::addToDiag(const vec3d &nvec, int ii, double FUvalue, double TWvalue) not implemented for TRILINOS yet ! " << endl;
	exit(1);

    int ii3 = 3*ii;
    int ii3_1 = ii3+1;
    int ii3_2 = ii3+2;
    int jj3 = 3*jj;
    int jj3_1 = jj3+1;
    int jj3_2 = jj3+2;
    
    double XA_n0 = XA*nvec.x;
    double XA_n1 = XA*nvec.y;
    double XA_n2 = XA*nvec.z;
    double XA_n1n0 = XA_n0*nvec.y;
    double XA_n2n1 = XA_n1*nvec.z;
    double XA_n0n2 = XA_n2*nvec.x;
	
    // declare ii and jj new columns, and update column nb
    int last_col_nb_ii = columns_nb[ii3];
    int last_col_nb_jj = columns_nb[jj3];
	
    columns[ii3  ][last_col_nb_ii  ] = jj3  ;
    columns[ii3  ][last_col_nb_ii+1] = jj3_1;
    columns[ii3  ][last_col_nb_ii+2] = jj3_2;
    columns[ii3_1][last_col_nb_ii  ] = jj3  ;
    columns[ii3_1][last_col_nb_ii+1] = jj3_1;
    columns[ii3_1][last_col_nb_ii+2] = jj3_2;
    columns[ii3_2][last_col_nb_ii  ] = jj3  ;
    columns[ii3_2][last_col_nb_ii+1] = jj3_1;
    columns[ii3_2][last_col_nb_ii+2] = jj3_2;
    
    columns[jj3  ][last_col_nb_jj  ] = ii3  ;
    columns[jj3  ][last_col_nb_jj+1] = ii3_1;
    columns[jj3  ][last_col_nb_jj+2] = ii3_2;
    columns[jj3_1][last_col_nb_jj  ] = ii3  ;
    columns[jj3_1][last_col_nb_jj+1] = ii3_1;
    columns[jj3_1][last_col_nb_jj+2] = ii3_2;
    columns[jj3_2][last_col_nb_jj  ] = ii3  ;
    columns[jj3_2][last_col_nb_jj+1] = ii3_1;
    columns[jj3_2][last_col_nb_jj+2] = ii3_2;
    
    columns_nb[ii3]   += 3;
    columns_nb[ii3_1] += 3;
    columns_nb[ii3_2] += 3;
    columns_nb[jj3]   += 3;
    columns_nb[jj3_1] += 3;
    columns_nb[jj3_2] += 3;
    
    // set values
    values[ii3  ][last_col_nb_ii  ] = XA_n0*nvec.x; // 00
    values[ii3  ][last_col_nb_ii+1] = XA_n1n0;      // 01
    values[ii3  ][last_col_nb_ii+2] = XA_n0n2;      // 02
    values[ii3_1][last_col_nb_ii  ] = XA_n1n0;      // 10
    values[ii3_1][last_col_nb_ii+1] = XA_n1*nvec.y; // 11
    values[ii3_1][last_col_nb_ii+2] = XA_n2n1;      // 12
    values[ii3_2][last_col_nb_ii  ] = XA_n0n2;      // 20
    values[ii3_2][last_col_nb_ii+1] = XA_n2n1;      // 21
    values[ii3_2][last_col_nb_ii+2] = XA_n2*nvec.z; // 22
    
    values[jj3  ][last_col_nb_jj  ] = XA_n0*nvec.x; // 00
    values[jj3  ][last_col_nb_jj+1] = XA_n1n0;      // 01
    values[jj3  ][last_col_nb_jj+2] = XA_n0n2;      // 02
    values[jj3_1][last_col_nb_jj  ] = XA_n1n0;      // 10
    values[jj3_1][last_col_nb_jj+1] = XA_n1*nvec.y; // 11
    values[jj3_1][last_col_nb_jj+2] = XA_n2n1;      // 12
    values[jj3_2][last_col_nb_jj  ] = XA_n0n2;      // 20
    values[jj3_2][last_col_nb_jj+1] = XA_n2n1;      // 21
    values[jj3_2][last_col_nb_jj+2] = XA_n2*nvec.z; // 22
}

void
StokesSolver::factorizeResistanceMatrix(){

	chol_L = cholmod_analyze(chol_res_matrix, &chol_c);
	cholmod_factorize(chol_res_matrix, chol_L, &chol_c);

    if (chol_c.status) {
		// Cholesky decomposition has failed: usually because matrix is incorrectly found to be positive-definite
		// It is very often enough to force another preconditioner to solve the problem.
		cerr << " factorization failed. forcing simplicial algorithm... " << endl;
		chol_c.supernodal = CHOLMOD_SIMPLICIAL;
		chol_L = cholmod_analyze (chol_res_matrix, &chol_c);
		cholmod_factorize (chol_res_matrix, chol_L, &chol_c);
		cerr << " factorization status " << chol_c.status << " final_ll ( 0 is LDL, 1 is LL ) " <<  chol_c.final_ll <<endl;
		chol_c.supernodal = CHOLMOD_SUPERNODAL;
    }
}

#ifdef TRILINOS
/*************  Preconditioners *********************/

/*
 buildDiagBlockPreconditioner() :
 A block diagonal (left-)preconditioner, almost similar to the
 one described in Amit Kumar's PhD Thesis: it is zero everywhere
 except along the diagonal where diagonal 3x3 block are the
 ones of R_FU:
 
        ..........
       |          .                             |
       | R_FU(i,j).     0                       |
       |          .                             |
       | .....................                  |
       |          .          .                  |
       |     0    . R_FU(i,j).                  |
       |          .          .                  |
       |          ............                  |
 P =   |                       .                |
       |                         .              |
       |                           .            |
       |                              ..........|
       |                             .          |
       |                             . R_FU(i,j)|
       |                             .          |
	                                  ...........
 
 This method stores P^{-1} in tril_l_precond.
 */
void
StokesSolver::setDiagBlockPreconditioner(){

    double a00, a01, a02, a11, a12, a22;
    double det, idet;
    double *precond_row = new double [3];
    int *indices = new int [3];
    for (int i = 0; i < np; i++) {
		int i3 = 3*i;
		
		indices[0] = i3;
		indices[1] = i3+1;
		indices[2] = i3+2;
		
	    // +2.5*r in the diagonal ? --> Amit Kumar, PhD Thesis
		a00 = values[i3][0];
		a01 = values[i3][1];
		a02 = values[i3][2];
		a11 = values[i3+1][1];
		a12 = values[i3+1][2];
		a22 = values[i3+2][2];
		
		det = a00*(a22*a11-a12*a12)+a01*(-a01*a22+2*a12*a02)-a02*a02*a11;
		idet = 1/det;
		
		// row i3
		precond_row[0] = idet*(a11*a22-a12*a12);
		precond_row[1] = idet*(a02*a12-a01*a22);
		precond_row[2] = idet*(a01*a12-a02*a11);
		
		tril_l_precond->InsertGlobalValues(i3, 3, precond_row, indices);
		
		// row i3+1
		precond_row[0] = precond_row[1]; // symmetric matrix!
		precond_row[1] = idet*(a00*a22-a02*a02);
		precond_row[2] = idet*(a02*a01-a00*a12);
		
		tril_l_precond->InsertGlobalValues(i3+1, 3, precond_row, indices);
		
		// row i3+2
		precond_row[1] = precond_row[2]; // symmetric matrix!
		precond_row[0] = idet*(a01*a12-a02*a11);
		precond_row[2] = idet*(a00*a11-a01*a01);
		
		tril_l_precond->InsertGlobalValues(i3+2, 3, precond_row, indices);
    }
	
    tril_l_precond->FillComplete();
	
	// give it to the LinearProblem
    tril_stokes_equation->setLeftPrec(rcp(tril_l_precond, false));
	
    delete [] precond_row;
    delete [] indices;
}

/*
 setIncCholPreconditioner() :
 A incomplete Cholesky factorization (left-)preconditioner.
 It uses Ifpack routines.
 
 Right now it seems that the routine sends back a diagonal preconditionner,
 with values being the inverse of the ones on R_FU diagonal.
 I (Romain) don't understand this behavior so far.
 */
void
StokesSolver::setIncCholPreconditioner(){
	//  parameters to be tuned
    int fill_level = 0;
	//    double drop_tolerance = 1.;
	
	RCP <Ifpack_IC> tril_ICT_precond = rcp(new Ifpack_IC(tril_res_matrix));
	
	ParameterList precondParams;
	//	precondParams.set("fact: drop tolerance", drop_tolerance);
	precondParams.set("fact: ict level-of-fill", fill_level);
	
	tril_ICT_precond->SetParameters(precondParams);
	tril_ICT_precond->Initialize();
	tril_ICT_precond->Compute();
	
	
	/*****	 TESTING *****
	 cout << " non zeros : " << tril_ICT_precond->NumGlobalNonzeros() << " " << tril_ICT_precond->IsInitialized() << " " << tril_ICT_precond->IsComputed() << endl;
	 int nb;
	 double *values = new double [res_matrix_linear_size];
	 int *indices = new int [res_matrix_linear_size];
	 
	 Epetra_CrsMatrix precU = tril_ICT_precond->U();
	 Epetra_Vector precD = tril_ICT_precond->D();
	 precD.ExtractCopy(values);
	 
	 cout << " Diagonal " << endl;
	 for(int i=0; i<res_matrix_linear_size; i++)
	 cout << i << " " << values[i] << endl;
	 
	 cout << " Upper " << endl;
	 for(int i=0; i<res_matrix_linear_size; i++){
	 precU.ExtractGlobalRowCopy(i, res_matrix_linear_size, nb, values, indices);
	 for(int j=0; j<nb; j++)
	 cout << i << " " << indices[j] << " " << values[j] << endl;
	 }
	 
	 cout << " Original Matrix Diagonal " << endl;
	 for(int i=0; i<res_matrix_linear_size; i++){
	 
	 tril_res_matrix->ExtractGlobalRowCopy(i, res_matrix_linear_size, nb, values, indices);
	 for(int j=0; j<nb; j++){
	 if(indices[j] == i )
	 cout << i << " " << indices[j] << " " << 1./values[j] << endl;
	 }
	 }
	 
	 delete [] values;
	 delete [] indices;
	 ***** END TESTING ********/
	
	// template conversion, to make Ifpack preconditioner compatible with belos
	RCP<Belos::EpetraPrecOp> belos_ICT_precond = rcp (new Belos::EpetraPrecOp(tril_ICT_precond));
	
	tril_stokes_equation->setLeftPrec(belos_ICT_precond);
}
#endif

#ifdef TRILINOS
void
StokesSolver::matrixChol2Tril(const cholmod_sparse *C, Epetra_CrsMatrix* &T){
	vector <double> row_values;
	vector <int> row_indices;
	for(int i=0; i< res_matrix_linear_size;i++){
		int nz = C->
		C->x
	}
}
#endif

#ifdef CHOLMOD_EXTRA
/*
 setSpInvPreconditioner() :
 A sparse inverse (left-)preconditioner.
 cholmod-extra routine by Jaakko Luttinen
 
 */
void
StokesSolver::setSpInvPreconditioner(){
	cholmod_sparse *sparse_inv = cholmod_spinv(chol_L, &chol_c);
	cholmod_free_sparse(&sparse_inv, &chol_c);
}
#endif

void
StokesSolver::setSolverType(string solver_type){

	if (solver_type == "direct") {
		_direct = true;
		_iterative = false;
	} else {
		if (solver_type == "iterative") {
			cerr << " Error : StokesSolver::setSolverType(string solver_type) : extended Linear algebra not implemented with TRILINOS yet. "<< endl;
			exit(1);
#ifdef TRILINOS
			_direct = false;
			_iterative = true;
#else
			cerr << " Error : StokesSolver::setSolverType(string solver_type) : 'iterative' solver asked, but needs to be compiled with Trilinos library."<< endl;
			exit(1);
#endif
		} else {
			cerr << " Error : StokesSolver::setSolverType(string solver_type) : Unknown solver type '" << solver_type << "'"<< endl;
			exit(1);
		}
	}
}

// testing
void
StokesSolver::printResistanceMatrix(){
	if (direct()) {
		//		cout << endl<< " chol res " << endl;
		for (int i = 0; i < res_matrix_linear_size; i++) {
			for (int k =((int*)chol_res_matrix->p)[i] ; k < ((int*)chol_res_matrix->p)[i+1]; k++) {
				cout << i << " " << ((int*)chol_res_matrix->i)[k] << " " << ((double*)chol_res_matrix->x)[k] << endl;
			}
		}
	}
	exit(1);
#ifdef TRILINOS
	if (iterative()) {
		int int_nb = 100;
		double *val = new double [int_nb];
		int *ind = new int [int_nb];
		int nz;
		cout << endl<< " tril res " << endl;
		for (int i = 0; i < res_matrix_linear_size; i++) {
			tril_res_matrix->ExtractGlobalRowCopy(i, int_nb, nz, val, ind);
			//	   cout << i << " " << nz << endl;
			for (int j = 0; j < nz; j++) {
				cout << i << " " << ind[j] << " " << val[j] << endl;
			}
		}
		// cout << "precond " << endl;
		// for(int i = 0; i < res_matrix_linear_size; i++){
		//   tril_l_precond->ExtractGlobalRowCopy(i, int_nb, nz, val, ind);
		//   cout << " line " << i << " " << nz << endl;
		//    for(int j = 0; j < nz; j++){
		//      cout << i << " " << ind[j] << " " << val[j] << endl;
		//    }
		// }
		delete [] val;
		delete [] ind;
	}
#endif
}
