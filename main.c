
	#include <stdio.h>
	#include <stdbool.h>
	#include <string.h>
	#include <math.h>
	#include <lapacke.h>
	#include "include/def.h"
	#include <time.h>

	int main(int argc, char **argv) {
		
		time_t t_0, t_1;
		t_0 = time(NULL);
		
		int length = 10;
		float g = 2;
		
		tsbasis * sbasis = init_sbasis(length);
		tvbasis * vbasis = init_vbasis(length);
		basis_states(sbasis, vbasis);

		tsmatrix * Hmatrix = hamiltonian(length, g);
		tspectrum * spectrum = eigenvalues(Hmatrix);
		float M = magnetization_expectation_value(vbasis,
							spectrum->eigvectors->data[0],
							length);

		t_1 = time(NULL);
		printf("Timelapse: %ld seconds\n", (t_1 - t_0));
		
		return 0;
	}
	
	int dim = (int)(pow(2, length)+1e-9);
	int numg = 100;
	float ming = 0;
	float maxg = 2;
	float stepg = (maxg - ming)/((float)numg);
	float * g_array = malloc((numg+1)*sizeof(float));*/
	float g = 2;
	
	FILE *fp;
	fp = fopen("data/data_fermion_g2.txt", "w");
	
	
	for (int i = 0; i < spectrum->length; i++) {
		fprintf(fp, "%f	%f\n", g, spectrum->eigvals[i]);
	}
	
/*	for (int i = 0; i < spectrum->length; i++) {
		for (int j = 0; j < dim; j++) {
			printf("%.0f	", Hmatrix->data[i][j]);
		}
		printf("\n");
		printf("%f\n", spectrum->eigvals[i]);
	}*/
	
	/*for (int j = 0; j < (numg + 1); j++) {
		g_array[j] = ming + stepg*j;
		tsmatrix * Hmatrix = H(sbasis, g_array[j]);
		tspectrum * spectrum = eigenvalues(Hmatrix);
		
		for (int i = 0; i < spectrum->length; i++) {
			//printf("%f	%f\n", M->data[i], g_array[i]);
			fprintf(fp, "%f	%f\n", g_array[j], spectrum->eigvals[i]);
		}
		free(Hmatrix);
		free(spectrum);
	}*/
	
	fclose(fp);
	
	t_1 = time(NULL);
	printf("Timelapse: %ld seconds\n", (t_1 - t_0));
	
	return 0;
}


/*	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < length; j++) {
			printf("%.0f	", sbasis->state_matrix->data[i][j]);
		}
		printf("\n");
	}
	printf("====================================\n");
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			printf("%.0f	", vbasis->vector_matrix->data[i][j]);
		}
		printf("\n");
	}
	
	for (int i = 0; i < dim; i++) {
		printf("%d	", vbasis->nup[i]);
	}
	printf("\n");
*/
/*	
	tspectrum * ground_spectrum = init_spectrum(numg+1, dim);
	tvector * M = init_vector(numg+1);
	tsmatrix ** Hmatrix_array = malloc(11*sizeof(tsmatrix));
	tspectrum ** spectrum_array = malloc(11*sizeof(tspectrum));
	
	for (int i = 0; i < (numg+1); i++) {
		g_array[i] = pow(10, ming + stepg*i);
		Hmatrix_array[i] = hamiltonian(length, g_array[i]);
		spectrum_array[i] = eigenvalues(Hmatrix_array[i]);
		ground_spectrum->eigvals[i] = spectrum_array[i]->eigvals[0];
		for (int j = 0; j < spectrum_array[i]->eigvectors->ncol; j++) {
			ground_spectrum->eigvectors->data[i][j] = spectrum_array[i]->eigvectors->data[0][j];
		}
		M->data[i] = magnetization_expectation_value(vbasis, ground_spectrum->eigvectors->data[i], length);
		//free(Hmatrix);
		//free(spectrum);
	}
	
	FILE *fp;
	fp = fopen("data/data6.txt", "w");
	
	for (int i = 0; i < (numg+1); i++) {
		printf("%.4f	%f\n", M->data[i], g_array[i]);
		fprintf(fp, "%f	%f\n", M->data[i], g_array[i]);
	}
	
	
	tsmatrix * Hmatrix;
	tspectrum * spectrum;
	Hmatrix = hamiltonian(length, g);
	spectrum = eigenvalues(Hmatrix);
	float Mground = magnetization_expectation_value(vbasis, spectrum->eigvectors->data[0], length);
	printf("=======\n%.4f	%.4f\n", Mground, g);
*/	
	//for (int i = 0; i < spectrum->length; i++) {
	//	printf("%.2f	", spectrum->eigvectors->data[0][i]);
	//}

/*	
	//for (int i = 0; i < spectrum->length; i++) {
		printf("%.4f	", spectrum->eigvals[0]);
	//}
	printf("\n==============\n");
	//for (int i = 0; i < spectrum->length; i++) {
		for (int j = 0; j < spectrum->length; j++) {
			printf("%.4f	", spectrum->eigvectors->data[0][j]);
		}
		printf("\n");
	//}
*/
	
	/*for (int i = 0; i < Hmatrix->nrow; i++) {
		for (int j = 0; j < Hmatrix->ncol; j++) {
			printf("%.2f	", Hmatrix->data[i][j]);
		}
		printf("\n");
	}*/

/*	tbasis * basis = init_basis(length);
	basis_states(basis);
	
	for (int i = 0; i < basis->nstates; i++) {
		for (int j = 0; j < basis->length; j++) {
			printf("%d	", basis->data[i]->spin[j]);
		}
		printf("\n");
	}
*/