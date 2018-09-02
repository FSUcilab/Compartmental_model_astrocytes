
void matrixMult(double* L, double *r, int n) {
	int i,j,k;
        //show_matrix(L);
  	//printf("\n");
  	//for (i=0; i < NPROC; i++) { printf("rinit: %f\n", r[i]); }


	double s;

	for ( i = 0; i < n; i++){	
		s = 0;
		for ( j = 0; j < (i+1); j++) {
			s = s + L[n*i+j]*r[j];
			//printf("j=%i,L[dd]=%f ",j,L[NPROC*i+j] );
		}
	r[i] = s;
	//printf("\n");
	}
  	//for (i=0; i < NPROC; i++) { printf("rfinal: %f\n", r[i]); }
}


