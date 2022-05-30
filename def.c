#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <lapacke.h>
#include <lapacke_utils.h>
#include "include/def.h"

/* SQUARE MATRIX
typedef struct {
	float ** data;		// (2D-array) matrix data
	int nrow;			// number of rows
	int ncol;			// number of columns
} tsmatrix;
*/

// INITIALIZE SQUARE MATRIX
tsmatrix * init_square_matrix(int size_row, int size_col) {
	tsmatrix * smatrix = malloc(sizeof(tsmatrix));					// (tsmatrix-malloc) create malloc of s-matrix
	
	smatrix->nrow = size_row;										// number of rows equal to row size
	smatrix->ncol = size_col;										// number of columns equal to column size
	smatrix->data = malloc(smatrix->nrow*sizeof(float *));			// (nrow*float_pointer-malloc) cr mall of row data
	
	for (int i = 0; i < smatrix->nrow; i++) {						// for each row
		smatrix->data[i] = malloc(smatrix->ncol*sizeof(float));		// (ncol*float-malloc) create malloc of columns data
	}
	
	return smatrix;
}

/* FLOAT VECTOR
typedef struct {
	float * data;		// (1D-array) vector data
	int length;			// vector length
} tvector;
*/

// INITIALIZE FLOAT VECTOR
tvector * init_vector(int size_vector) {
	tvector * vector = malloc(sizeof(tvector));				// (tvector-malloc) create malloc of tvector
	
	vector->length = size_vector;							// vector length equal to vector size
	vector->data = malloc(vector->length*sizeof(float));	// (length*float-malloc) create malloc of vector data

	return vector;
}

/*	SBASIS
	typedef struct {
		tsmatrix * state_matrix;		// (2D-array) basis states
		int nstates;					// number of states
		int length;						// chain length
	} tsbasis;
*/

// INITIALIZE SBASIS
tsbasis * init_sbasis(int chain_length) {
	tsbasis * sbasis = malloc(sizeof(tsbasis));
	
	sbasis->length = chain_length;
	sbasis->nstates = (int)(pow(2, chain_length)+1e-9);
	sbasis->data = malloc(sbasis->nstates*sizeof(tvector));
	for (int i = 0; i < sbasis->nstates; i++) {
		sbasis->data[i] = init_vector(chain_length);
	}
	
	return sbasis;
}

/*	VBASIS
	typedef struct {
		tsmatrix * vector_matrix;		// (2D-array) basis vectors
		int nvectors;					// number of vectors
		int * nup;						// number of spin-up
	} tvbasis;
*/

// INITIALIZE VBASIS
tvbasis * init_vbasis(int chain_length) {
	tvbasis * vbasis = malloc(sizeof(tvbasis));					// (tbasis-malloc) create malloc of tbasis
	
	vbasis->nvectors = (int)(pow(2, chain_length)+1e-9);
	vbasis->nup = malloc(vbasis->nvectors*sizeof(int));								// number of states equal to basis size
	vbasis->vector_matrix = init_square_matrix(vbasis->nvectors, vbasis->nvectors);		// (nstates*tstate-malloc) create malloc of states data
	
	return vbasis;
}

/* SPECTRUM
typedef struct {
	float * eigvals;			// (1D-array) eigenvalues
	tsmatrix * eigvectors;		// (tsmatrix-malloc) eigenvectors
	int length;					// (int) number of eigenvalues
} tspectrum;
*/

// INITIALIZE SPECTRUM
tspectrum * init_spectrum(int n_eigvals, int dim) {
	tspectrum * spectrum = malloc(sizeof(tspectrum));				// malloc of tspectrum
	
	spectrum->length = n_eigvals;									// length = number of eigenvalues
	spectrum->eigvals = malloc(spectrum->length*sizeof(float));		// (length*float) malloc of eigvals array
	spectrum->eigvectors = init_square_matrix(spectrum->length, dim);	// smatrix of eigvectors (sz =lngth)
	
	return spectrum;
}

tsmatrix * identity_matrix(int dim) {
	tsmatrix * identity = init_square_matrix(dim, dim);
	
	for (int i = 0; i < dim; i++) {
		identity->data[i][i] = 1;
	}
	
	return identity;
}

tsmatrix * tensor_matrix_product(tsmatrix *A, tsmatrix *B) {
	tsmatrix * C = init_square_matrix(A->nrow*B->nrow, A->ncol*B->ncol);
	for (int ia = 0; ia < A->nrow; ia++) {
		for (int ja = 0; ja < A->ncol; ja++) {
			for (int ib = 0; ib < B->nrow; ib++) {
				for (int jb = 0; jb < B->ncol; jb++) {
					C->data[ib+ia*B->nrow][jb+ja*B->ncol] = A->data[ia][ja] * B->data[ib][jb];
				}
			}
		}
	}
	
	return C;
}

void tensor_vector_product(tvector * vector_in, int v2_index) {
	
	tvector * vector_out = init_vector(2*vector_in->length);
	
	if (v2_index == 1) {
		for (int i = 0; i < vector_in->length; i++) {
			vector_out->data[2*i] = vector_in->data[i];
			vector_out->data[2*i+1] = 0;
		}
	} else {
		for (int i = 0; i < vector_in->length; i++) {
			vector_out->data[2*i+1] = vector_in->data[i];
			vector_out->data[2*i] = 0;
		}
	}
	
	for (int i = 0; i < vector_out->length; i++) {
		vector_in->data[i] = vector_out->data[i];
	}
	vector_in->length = vector_out->length;
}

// BASIS STATES

	void basis_states(tsbasis *sbasis, tvbasis *vbasis) {
		int n;
		
		for (int i = 0; i < sbasis->nstates; i++) {
			n = i;
			for (int j = 0; j < sbasis->length; j++) {
				sbasis->data[i]->data[j] = n%2;
				n /= 2;
			}
		}
		
		for (int i = 0; i < sbasis->nstates; i++) {
			tvector * aux_vector = init_vector(sbasis->nstates);
			
			int nup = 0;
			
			if (sbasis->data[i]->data[0] == 1) {
				aux_vector->data[0] = 1;
				nup = 1;
			} else {
				aux_vector->data[1] = 1;
			}
			aux_vector->length = 2;
			
			for (int j = 1; j < sbasis->length; j++) {
				tensor_vector_product(aux_vector, sbasis->data[i]->data[j]);
				nup += sbasis->data[i]->data[j];
			}
			
			int pos1;
			for (int j = 0; j < vbasis->nvectors; j++) {
				if (aux_vector->data[j] == 1) {
					pos1 = j;
					vbasis->nup[j] = nup;
				}
			}
			
			for (int j = 0; j < vbasis->nvectors; j++) {
				vbasis->vector_matrix->data[pos1][j] = aux_vector->data[j];
			}
			free(aux_vector);
		}
	}

void sum_matrix(tsmatrix *A, tsmatrix *B, float g) {
	if (A->nrow != B->nrow) {
		if (A->ncol != B->ncol) {
			printf("ERROR\n");
			exit(1);
		}
	}
	
	for (int i = 0; i < A->nrow; i++) {
		for (int j = 0; j < A->ncol; j++) {
			A->data[i][j] += g*B->data[i][j];
		}
	}
}

// MATRIX - ARRAY PRODUCT
float * matrix_array_product(tsmatrix *matrix, float *array) {
	float * result = malloc(matrix->nrow*sizeof(float));
	for (int i = 0; i < matrix->nrow; i++) {
		for (int j = 0; j < matrix->ncol; j++) {
			result[i] += matrix->data[i][j]*array[j];
		}
	}
	
	return result;
}

// SCALAR PRODUCT
float scalar_product(float *array_1, float *array_2, int length) {
	float result = 0;
	for (int i = 0; i < length; i++) {
		result += array_1[i]*array_2[i];
	}
	
	return result;
}

	tsmatrix * hamiltonian(int length, float g) {
		int a;
		int dim = (int)(pow(2, length)+1e-9);
		tsmatrix * Hmatrix = init_square_matrix(dim, dim);
		
		tsmatrix * sigma_zz = init_square_matrix(4,4);
		sigma_zz->data[0][0] = -1;
		sigma_zz->data[1][1] = 1;
		sigma_zz->data[2][2] = 1;
		sigma_zz->data[3][3] = -1;
		
		sum_matrix(Hmatrix, tensor_matrix_product(sigma_zz, identity_matrix((int)(pow(2, length-2)+1e-9))), 1);
		for (int i = 1; i < (length-2); i++) {
			sum_matrix(Hmatrix, tensor_matrix_product(
								tensor_matrix_product(identity_matrix((int)(pow(2, i)+1e-9)), sigma_zz),
								identity_matrix((int)(pow(2, length-i-2)+1e-9))),
								1);
		}
		sum_matrix(Hmatrix, tensor_matrix_product(identity_matrix((int)(pow(2, length-2)+1e-9)), sigma_zz), 1);

		tsmatrix * sigma_x = init_square_matrix(2,2);
		sigma_x->data[1][0] = -1;
		sigma_x->data[0][1] = -1;
		
		sum_matrix(Hmatrix, tensor_matrix_product(sigma_x, identity_matrix((int)(pow(2, length-1)+1e-9))), g);
		for (int i = 1; i < (length-1); i++) {
			sum_matrix(Hmatrix, tensor_matrix_product(
								tensor_matrix_product(identity_matrix((int)(pow(2, i)+1e-9)), sigma_x),
								identity_matrix((int)(pow(2, length-i-1)+1e-9))),
								g);
		}
		sum_matrix(Hmatrix, tensor_matrix_product(identity_matrix((int)(pow(2, length-1)+1e-9)), sigma_x), g);
		
		return Hmatrix;
	}

/*	for (int k = 0; k < (length-1); k++) {
		a = 0;
		for (int i = 0; i < (int)(pow(2, k)+1e-9); i++) {
			for (int j = 0; j < (int)(pow(2, length-k-2)+1e-9); j++) {
				Hmatrix->data[a][a] += -1;
				a++;
			}
			for (int j = 0; j < (int)(pow(2, length-k-1)+1e-9); j++) {
				Hmatrix->data[a][a] -= -1;
				a++;
			}
			for (int j = 0; j < (int)(pow(2, length-k-2)+1e-9); j++) {
				Hmatrix->data[a][a] += -1;
				a++;
			}
		}
	}*/
/*	int dimk, dimNk;
	for (int k = 0; k < length; k++) {
		dimk = (int)(pow(2, k)+1e-9);
		dimNk = (int)(pow(2,length-k-1)+1e-9);
		for (int j = 0; j < dimk; j++) {
			for (int i = 0; i < dimNk; i++) {
				Hmatrix->data[2*j*dimNk+i+dimNk][2*j*dimNk+i] += -1;
				Hmatrix->data[2*j*dimNk+i][2*j*dimNk+i+dimNk] += -1;
			}
		}
	}
*/

	// EIGENVALUES AND VECTORS OF ONE MATRIX
	
	tspectrum * eigenvalues(tsmatrix *Hmatrix) {
		int dim = Hmatrix->nrow;
		tspectrum * spectrum = init_spectrum(dim, dim);		// (tspectrum) initialize spectrum
		int INFO = 1;		// (int) out-parameter of LAPACKE_ssyev function
								// = 0: successful exit
								// < 0: if INFO = -i, the i-th argument had an illegal value
								// > 0: if INFO = i, the algorithm failed to converge; i off-diagonal elements of an
									// intermediate tridiagonal form did not converge to zero
		tvector * matrix_array = init_vector(Hmatrix->ncol*Hmatrix->nrow);// (ncol*nrow_Hmatrix vector) conversion
																				// from tsmatrix to float-array
		
		for (int i = 0; i < Hmatrix->nrow; i++) {								// for each row
			for (int j = 0; j < Hmatrix->ncol; j++) {							// and for each column
				matrix_array->data[i*Hmatrix->nrow+j] = Hmatrix->data[i][j];	// insert matrix elements into the
				//printf("%f	", matrix_array->data[i*Hmatrix->nrow+j]);			// float-array matrix_array->data
			}
		}
		
		//LAPACKE_ssyev(int matrix_layout, char jobz, char uplo, lapack_int n, float* a, lapack_int lda, float* w);
		INFO = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'V', 'L', dim, matrix_array->data, dim, spectrum->eigvals);
		if (INFO > 0) {
			printf("the algorithm failed to compute eigenvalues.\n");
			exit(1);
		}
		
		// copy eigenvectors from matrix_array->data[] to eigvector_matrix->data[][]
		for (int i = 0; i < Hmatrix->nrow; i++) {
			for (int j = 0; j < Hmatrix->ncol; j++) {
				spectrum->eigvectors->data[i][j] = matrix_array->data[i*Hmatrix->nrow+j];
			}
		}
		
		return spectrum;
	}

	// (int) matrix_layout: indicates whether the input and output matrices are stored in row mayor order or column
			// major order, where:
				// = LAPACK_ROW_MAJOR: the matrices are stored in row major order (C/C++ default)
				// = LAPACK_COL_MAJOR: the matrices are stored in column major order (Fortran default)
		// (char) jobz: indicates the type of computation to be performed, where:
				// = 'N': eigenvalues only are computed
				// = 'V': eigenvalues and eigenvectors are computed
		// (char) uplo: indicates whether the upper or lower part of matrix A is referenced, where:
				// = 'U': the upper triangular part is referenced
				// = 'L': the lower triangular part is referenced
		// (lapack_int) n: is the order of matrix A used in computation (n >= 0)
		// (float_pointer) a: is the real symmetric or complex Hermitian matrix A of order n
		// (lapack_int) lda: is the leading dimension of the array specified for A (lda >= n)
		// (float_pointer) w: is the vector w of length n, which contains the computed eigenvalues in acending order

	float magnetization_expectation_value(tvbasis *vbasis, float *state, int length) {
		float M = 0;
		
		for (int i = 0; i < vbasis->nvectors; i++) {
			M += pow(state[i], 2)*abs(2*vbasis->nup[i]- length);
		}
		
		M /= (float) length;
		
		return M;
	}

////////////////////////////////  FERMION SYSTEM  ////////////////////////////////

int compare_states(tvector *vi, tvector *vj) {
	int result = 1;
	for (int i = 0; i < vi->length; i++) {
		if (vi->data[i] != vj->data[i]) {
			result = 0;
		}
	}
	
	return result;
}

void copy_vector(tvector *v_out, tvector *v_in) {
	for (int i = 0; i < v_in->length; i++) {
		v_out->data[i] = v_in->data[i];
	}
}

tvector * Cact(tvector *state, int index, int act[4]) {
	tvector * result = init_vector(state->length);
	copy_vector(result, state);
	
	if (act[0] == 1) {
		result->data[index] -= 1;
	}
	if (act[1] == 1) {
		result->data[index] += 1;
	}
	if (act[2] == 1) {
		result->data[index+1] -= 1;
	}
	if (act[3] == 1) {
		result->data[index+1] += 1;
	}
	
	return result;
}

float Hact(tvector *vi, tvector *vj, float g) {
	float H = 0;
	
	int act1[4] = {1,1,0,0};
	int act2[4] = {1,0,1,0};
	int act3[4] = {0,1,0,1};
	int act4[4] = {0,1,1,0};
	int act5[4] = {1,0,0,1};
	
	for (int i = 0; i < (vi->length-1); i++) {
		H -= g;
		H -= -2*g*compare_states(Cact(vi, i, act1), vj);
		H -= compare_states(Cact(vi, i, act2), vj);
		H -= compare_states(Cact(vi, i, act3), vj);
		H -= compare_states(Cact(vi, i, act4), vj);
		H -= compare_states(Cact(vi, i, act5), vj);
	}
	H -= g;
	H -= -2*g*compare_states(Cact(vi, (vi->length-1), act1), vj);
	
	return H;
}

tsmatrix * H(tsbasis *basis, float g) {
	tsmatrix * Hmatrix = init_square_matrix(basis->nstates, basis->nstates);				// (tsmatrix) H matrix
	for (int i = 0; i < Hmatrix->nrow; i++) {												// for each row
		for (int j = 0; j < i; j++) {														// and for each column
			Hmatrix->data[i][j] = Hact(basis->data[i], basis->data[j], g);	// Hact on ch pair of states
			Hmatrix->data[j][i] = Hmatrix->data[i][j];
		}
		Hmatrix->data[i][i] = Hact(basis->data[i], basis->data[i], g);
	}
	
	return Hmatrix;
}

//tsmatrix * H_BdG(int length, float g) {
	
//}

/*tsmatrix * magnetization(int length) {
	int a;
	int dim = (int)(pow(2, length)+1e-9);
	tsmatrix * Mmatrix = init_square_matrix(dim, dim);
	
	for (int k = 0; k < length; k++) {
		a = 0;
		for (int i = 0; i < (int)(pow(2, k)+1e-9); i++) {
			for (int j = 0; j < (int)(pow(2, length-k-1)+1e-9); j++) {
				Mmatrix->data[a][a] += 1;
				a++;
			}
			for (int j = 0; j < (int)(pow(2, length-k-1)+1e-9); j++) {
				Mmatrix->data[a][a] -= 1;
				a++;
			}
		}
	}
	
	return Mmatrix;
}*/

/*void copy_matrix(tsmatrix *matrix_out, tsmatrix *matrix_in) {
	for (int i = 0; i < matrix_in->nrow; i++) {
		for (int j = 0; j < matrix_in->ncol; j++) {
			matrix_out->data[i][j] = matrix_in->data[i][j];
		}
	}
	matrix_out->nrow = matrix_in->nrow;
	matrix_out->ncol = matrix_in->ncol;
}

void JW_operators(tsmatrix *ci_matrix, tsmatrix *ci_plus_matrix, int length, int i) {
	int dimi = (int)(pow(2, i+1)+1e-9);
	
	tsmatrix * aux1_matrix = init_square_matrix(dimi, dimi);
	aux1_matrix->nrow = 2;
	aux1_matrix->ncol = 2;
	aux1_matrix->data[0][0] = 1;
	aux1_matrix->data[1][1] = -1;
	
	tsmatrix * sigma_plus = init_square_matrix(2,2);
	tsmatrix * sigma_minus = init_square_matrix(2,2);
	sigma_plus->data[0][1] = 1;
	sigma_minus->data[1][0] = 1;
	
	tsmatrix * aux_plus_matrix = init_square_matrix(dimi, dimi);
	tsmatrix * aux_minus_matrix = init_square_matrix(dimi, dimi);
	
	if (i == 0) {
		copy_matrix(aux_plus_matrix, sigma_plus);
		copy_matrix(aux_minus_matrix, sigma_minus);
	} else {
		for (int j = 1; j < i; j++) {
			tsmatrix * aux2_matrix = init_square_matrix(2*aux1_matrix->nrow, 2*aux1_matrix->ncol);
			aux2_matrix = tensor_matrix_product(aux1_matrix, sigma_z);
			copy_matrix(aux1_matrix, aux2_matrix);
			free(aux2_matrix);
		}
		aux_plus_matrix = tensor_matrix_product(aux1_matrix, sigma_plus);
		aux_minus_matrix = tensor_matrix_product(sigma_minus, aux1_matrix);
		
	}
	
	if (i < (length-1)) {
		tsmatrix * identity = identity_matrix((int)(pow(2, length-(i+1))+1e-9));
		ci_matrix = tensor_matrix_product(aux_plus_matrix, identity);
		ci_plus_matrix = tensor_matrix_product(aux_minus_matrix, identity);
	} else {
		copy_matrix(ci_matrix, aux_plus_matrix);
		copy_matrix(ci_plus_matrix, aux_minus_matrix);
	}
}


tsmatrix * JW_hamiltonian(int length, float g) {
	int dim = (int)(pow(2, length)+1e-9);
	tsmatrix * JW_Hmatrix = init_square_matrix(dim, dim);
	tsmatrix * ci_matrix = init_square_matrix(dim, dim);
	tsmatrix * ci_plus_matrix = init_square_matrix(dim, dim);
	tsmatrix * ci_1_matrix = init_square_matrix(dim, dim);
	tsmatrix * ci_1_plus_matrix = init_square_matrix(dim, dim);
	tsmatrix * identity = identity_matrix(dim);
	
	for (int i = 0; i < length; i++) {
		sum_matrix(JW_Hmatrix, identity, -g);
		JW_operators(ci_matrix, ci_plus_matrix, length, i);
		JW_operators(ci_1_matrix, ci_1_plus_matrix, length, i+1);
		sum_matrix(JW_Hmatrix, matrix_matrix_product(ci_matrix, ci_plus_matrix), 2*g);
		sum_matrix(JW_Hmatrix, matrix_matrix_product(ci_1_matrix, ci_matrix), -1);
		sum_matrix(JW_Hmatrix, matrix_matrix_product(ci_plus_matrix, ci_1_plus_matrix), -1);
		sum_matrix(JW_Hmatrix, matrix_matrix_product(ci_plus_matrix, ci_1_matrix), -1);
		sum_matrix(JW_Hmatrix, matrix_matrix_product(ci_1_plus_matrix, ci_matrix), -1);
	}
	
	return JW_Hmatrix;
}*/
