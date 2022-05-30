#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <lapacke.h>

// SQUARE MATRIX
typedef struct {
	float ** data;		// (2D-array) matrix data
	int nrow;			// number of rows
	int ncol;			// number of columns
} tsmatrix;

// FLOAT VECTOR
typedef struct {
	float * data;		// (1D-array) vector data
	int length;			// vector length
} tvector;

// SBASIS
typedef struct {
	tvector ** data;
	int nstates;
	int length;
} tsbasis;

// VBASIS
typedef struct {
	tsmatrix * vector_matrix;		// (2D-array) basis states
	int nvectors;		// number of states
	int * nup;
} tvbasis;

// SPECTRUM
typedef struct {
	float * eigvals;			// (1D-array) eigenvalues
	tsmatrix * eigvectors;		// (tsmatrix-malloc) eigenvectors
	int length;					// (int) number of eigenvalues
} tspectrum;

tsmatrix * init_square_matrix(int size_row, int size_col);
tvector * init_vector(int size_vector);
tsbasis * init_sbasis(int chain_length);
tvbasis * init_vbasis(int chain_length);
tspectrum * init_spectrum(int n_eigvals, int dim);
tsmatrix * identity_matrix(int dim);
tsmatrix * tensor_matrix_product(tsmatrix *A, tsmatrix *B);
void tensor_vector_product(tvector * vector_in, int v2_index);
void basis_states(tsbasis *sbasis, tvbasis *vbasis);
void sum_matrix(tsmatrix *A, tsmatrix *B, float g);
float * matrix_array_product(tsmatrix *matrix, float *array);
float scalar_product(float *array_1, float *array_2, int length);
tsmatrix * hamiltonian(int length, float g);
tspectrum * eigenvalues(tsmatrix *Hmatrix);
float magnetization_expectation_value(tvbasis *vbasis, float *state, int length);

int compare_states(tvector *vi, tvector *vj);
void copy_vector(tvector *v_out, tvector *v_in);
tvector * Cact(tvector *state, int index, int act[4]);
float Hact(tvector *vi, tvector *vj, float g);
tsmatrix * H(tsbasis *basis, float g);