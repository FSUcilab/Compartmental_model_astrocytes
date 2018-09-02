
/******************************************************************************/ 
void astro_d ( double t, double y[], double yp[] )

/******************************************************************************/
/*
  Purpose:

    astro_d evaluates the derivative for the ODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 August 2016

  Author:

    John Burkardt

  Modifier: 

     Evan Cresswell

  Parameters:

    Input, double T, the value of the independent variable.

    Input, double Y(NEQN), the value of the dependent variable.

    Output double YP(NEQN), the value of the derivative dY(1:NEQN)/dT.
*/
{
  double c4 = 50.0;
  double r = .31;
  //double der = 0.001; // was .001 10/18/2016
  int ii ;
  double dr[NPROC];
  for(ii = 0;ii<NPROC; ii++) 
  {
	  dr[ii] = r4_normal_ab ( mu, sigma, &seed );
	  if(ii==0)
	  {
	  	dr[ii] *= addOn;
	  }
	if (dr[ii]>1e5)
	{
		printf("oh shit");
		exit(0);
	}
  }

  matrixMult(L, dr, NPROC); // correlation step

  //if (I_step % 10)
  {
    if (dr[0]>1e5)
    {
	printf("oh shit");
	exit(0);
    } 
    fwrite(dr,sizeof(double),NPROC,   PTR_MYFILE);
  } 
  double cytosol_diff = 0;
  double der_diff = 0;
#ifdef VD
  double Vd_surf = pow(Vd,2./3);
#endif
  //astrocyte compartment equations....

  
  // MODEL
  int i;
  printf("t: %f,(r,neuro_input, dr): ", t);
  for(i=0;i<NPROC;i++) {
  	 //if (i == 2) { r = 0.;} else { r = 0.31; }
  	 //r = 0.;  // there should be no soma spikes
  	 yp[2*i] = (neuro_input[i] * dr[i]/sqrt(dt) + r - c4 * f(y[2*i],y[2*i+1]) - y[2*i] + dc*(y[NEQN-2] - y[2*i])) / tau  ; //process i cytosol
	 cytosol_diff += dc*(y[2*i] - y[NEQN-2]); //process i cytosol diffusion term 
  	 yp[2*i+1] = (f(y[2*i],y[2*i+1]) + der*(y[NEQN-1] - y[2*i+1])) / (eps*tau); //process i er
	 der_diff += der*(y[2*i+1] - y[NEQN-1]); //process i er diffusioh term 

  	// In analysis code, do not forget to divide the current terms by eps when used for ER calcium balance
  	serca_current[i]  = fserca(y[2*i],y[2*i+1]);
  	cicr_current[i]   = fcicr(y[2*i],y[2*i+1]);
  	linear_current[i] = flinear(y[2*i],y[2*i+1]);
    ca_removal[i]     = -y[2*i];
  	diffc[i]   = dc*(y[NEQN-2]  - y[2*i]);
  	differ[i]  = der*(y[NEQN-1] - y[2*i+1]);
  }

	
  diffc[NGROUPS-1] = cytosol_diff;
  differ[NGROUPS-1] = der_diff;
  serca_current[NGROUPS-1]  = fserca(y[NEQN-2],y[NEQN-1]);
  cicr_current[NGROUPS-1]   = fcicr(y[NEQN-2],y[NEQN-1]);
  linear_current[NGROUPS-1] = flinear(y[NEQN-2],y[NEQN-1]);
  ca_removal[NGROUPS-1] = -y[NEQN-2];

  // GE: May 1, 2018: turn soma off past a certain time
  double rsoma = r;
  //rsoma = 0.;
  if (t > t_rsoma_turnoff) rsoma = 0.;
  //printf("t= %f\n", t);
  //printf("rsoma= %f\n", rsoma);exit(0);
  printf("t= %f, rsoma= %f\n", t, rsoma);

#ifdef VD
  yp[NEQN-2] = ( (rsoma -c4 * f(y[NEQN-2],y[NEQN-1]) - y[NEQN-2]) * Vd_surf + cytosol_diff) / (Vd*tau); //soma cytosol
  yp[NEQN-1] = (f(y[NEQN-2],y[NEQN-1]) * Vd_surf + der_diff) / (eps*tau*Vd); //soma er
#else
  yp[NEQN-2] = (rsoma - c4 * f(y[NEQN-2],y[NEQN-1]) - y[NEQN-2] + cytosol_diff) / tau; //soma cytosol
  yp[NEQN-1] = (f(y[NEQN-2],y[NEQN-1]) + der_diff) / (eps*tau); //soma er
  printf("  not computing with Vd though you should be   ");
#endif
//-----------------------------------




  return;
}

//----------------------------------------------------------------------
double *cholesky(double *A) {
    int i,j,k;
    double *L = (double*)calloc(NPROC * NPROC, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);
 
    for ( i = 0; i < NPROC; i++)
        for ( j = 0; j < (i+1); j++) {
            double s = 0;
            for ( k = 0; k < j; k++)
                s += L[i * NPROC + k] * L[j * NPROC + k];
            L[i * NPROC + j] = (i == j) ?
                           sqrt(A[i * NPROC + i] - s) :
                           (1.0 / L[j * NPROC + j] * (A[i * NPROC + j] - s));
        }
 
    return L;
}
 
void show_matrix(double *A) {
	int i,j;
    for (i = 0; i < NPROC; i++) {
        for (j = 0; j < NPROC; j++){
            printf("%2.5f ", A[i * NPROC + j]);
	}
        printf("\n");
    }

}

