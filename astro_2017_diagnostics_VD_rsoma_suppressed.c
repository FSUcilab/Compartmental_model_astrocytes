// GE:  Add diagnostics to track different components of right hand side of equations
// I will only track equation for the processes and the soma equation (for now)
// This is a diagnostic tool to understand where CER is going as one increases DER. 

// clamp soma signal at 0.4 (to prevent it from spiking)
// leave the soma ER as is. 

#define VD 1

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include "normal.h"
# include <unistd.h>

int I_step;


void euler(  double y[],double yp[],double t, int neqn,double dt);
void writeIC(double y[],double t);
double readIC(double y[],double t);
void show_matrix(double *A); 
double *cholesky(double *A); 
void matrixMult(double* L, double *r, int n);


void astro (int N_step,double dur);
void astro_d ( double t, double y[], double yp[] );
double f(double ca, double cer);
FILE *PTR_MYFILE;
FILE *PTR_BALANCE; //GE
double dt;
double der;
double mu;
double dc; //was .0025;
double sigma;
FILE *PTR_IC;
int seed;

//defined in preprocesses
# ifndef NPROC 
	# define NPROC 5
# endif
# define NEQN 2*NPROC+2 
# define NGROUPS NPROC+1 

// neuro_input[i]: multiplicative factor for neuronal input into processes. 
// This allows me to set one of them to zero.
double neuro_input[NPROC];




int main ( int argc, char *argv[] );
double diag = 1.0;
double correlation;
double C[NPROC*NPROC];
double *L; 
//double alpha;
double addOn;
double tau;
double eps;

#ifdef VD
double Vd;
#endif
double t_rsoma_turnoff;  // May 1, 2018 (time at which to turn r off in soma)


// Gordon Erlebacher (GE) to do budget analysis
// one group per cytosol/ER pair
// I will leave eps and tauc on left hand side. They have to be taken into accounting in budget balance
double serca_current[NGROUPS];
double cicr_current[NGROUPS];
double linear_current[NGROUPS];
double ca_removal[NGROUPS];
double diffc[NGROUPS];
double differ[NGROUPS];
double source[NGROUPS];


/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:


  Discussion:

  Modified:

    17 June 2006

  Author:

*/

{

  //alpha = .1666666666666666666666666;
  tau = 8;
  eps = .04;

	for (int i=0; i < NPROC; i++) {
		neuro_input[i] = 1.;
	}


  int i,j,k;
  /* Create File */
  int N_step;
  double dur;
  dur  = atof(argv[1]);

  N_step = atoi(argv[2]);
  printf("N_step=%i\n",N_step);

  mu = 0.0;
  der= atof(argv[3]);
  //der = der * eps* tau/alpha; // 7/17/2017

  sigma= atof(argv[4]);
  //printf("alpha %f\n\n",alpha);
  //sigma = sigma*sqrt(alpha);// 1/17/2017
  //printf("\nsigma is %f\n",sigma);

  dc = atof(argv[5]);
  //dc = dc * tau/alpha;// 7/17/2017

  correlation = atof(argv[6]);
  correlation = correlation > .99999 ? 0.99999 : correlation;
  printf("correlation is %f",correlation);


  addOn = atof(argv[7]);

  PTR_MYFILE = fopen(argv[8], "wb");
  PTR_BALANCE = fopen("concentration_balance.bin", "wb");
  printf("\n %s \n",argv[8]);

  seed = atoi(argv[9]);

  printf("\n\nder is %f\n\n", der);
  printf("\n\ndc is %f\n\n", dc);
  printf("\n\nsigma is %f\n\n", sigma);

// evan's code breaks if VD is 0. Mine doesn't. That is because if VD==0, the argument positions are incorrect
#ifdef VD
  Vd = atof(argv[10]);
  printf("Vd is %f\n",Vd);
#endif

  t_rsoma_turnoff = atof(argv[11]);
  printf("t_rsoma_turnoff= %f\n", t_rsoma_turnoff);

  neuro_input[0] = atof(argv[12]);
  printf("neuro_input[0]= %f\n", neuro_input[0]);


printf("NPROC is  %d\n",NPROC);
printf("NEQN is %d",NEQN);

  for(i=0;i<NPROC;i++){
	for(j=0;j<NPROC;j++){
		C[i*NPROC+j] = correlation;
		if(i==j){
			C[i*NPROC+j] = diag;
		}
	}
  }

  // Generate covariance matrix..
  //  double L[NPROC][NPROC];
  printf("\n");

  printf("\n");
  //free(L);
  L = cholesky(C);
  printf("\n");
	
  
  //show_matrix(C);

  //show_matrix(L);
  astro (N_step,dur);

  return 0;


}

void astro (int N_step,double dur )
/******************************************************************************/
{

  //int i_step;
  //double dur = 5000.0;
  double relerr;
  double t;
  double t_out;
  double t_start;
  double t_stop;
  double y[NEQN];
  double yp[NEQN];


  # define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))

  t_start = 0.0;
  t_stop = t_start + dur;

  //n_step = t_stop*5;
  dt = (t_stop - t_start)/N_step;
//  dt = dt*alpha; //7/18/2017

  t = 0.0;
  t_out = 0.0;

  if (access( "IC.bin", F_OK ) != -1 ) {
          printf("t=%f\n",t);
          t = readIC(y, t);
          printf("t=%f\n",t);         
          //update t_start and t_stop
          t_start = t; 
          printf("t_start is %f\n", t_start);
          t_stop = t_start + dur;
          printf("t_stop is %f\n", t_stop);
  } else {
	int i;
  	for(i=0;i<NPROC;i++)
	{
		y[2*i]=.31;
		y[2*i+1]=1.0;
	}
	y[NEQN-2]=.31;
	y[NEQN-1]=1.0; 
 }


  if(!PTR_MYFILE) {
	  printf("Unable to open file!\n");
	  return;
  } else { 
  	  printf("seed is: %i\n",seed);

	  for ( I_step = 1; I_step <= N_step; I_step++ ) {
	      t = ( ( double ) ( N_step - I_step + 1 ) * t_start 
		    + ( double ) (          I_step - 1 ) * t_stop ) 
		    / ( double ) ( N_step              );

		  //t = ((N_step-I_step+1.) * t_start + (I_step-1.)*t_stop) / N_step;

	      t_out = ( ( double ) ( N_step - I_step ) * t_start 
		    + ( double ) (	   I_step ) * t_stop ) 
		    / ( double ) ( N_step );
            
          fwrite (&t ,sizeof(double), 1, PTR_MYFILE);
          fwrite(y,sizeof(double), NEQN, PTR_MYFILE);

	      if (y[0]>1.0e5){
		    printf("step %d\n",I_step);
	        printf("y=%f\n",y[0]);
    	    printf("yp=%f\n",yp[0]);
		    exit(0);
          }

	     euler (  y, yp, t, NEQN,dt );
         if (I_step==N_step) {
	        /* Create Initial Condition File*/
	        writeIC(y,t);
         }

		 fwrite(&t, sizeof(double), 1, PTR_BALANCE);
		 fwrite(serca_current,  sizeof(double), NGROUPS, PTR_BALANCE);
		 fwrite(cicr_current,   sizeof(double), NGROUPS, PTR_BALANCE);
		 fwrite(linear_current, sizeof(double), NGROUPS, PTR_BALANCE);
		 fwrite(ca_removal, sizeof(double), NGROUPS, PTR_BALANCE);
		 fwrite(diffc,  sizeof(double), NGROUPS, PTR_BALANCE);
		 fwrite(differ, sizeof(double), NGROUPS, PTR_BALANCE);
	  }
        
      fclose(PTR_MYFILE);
      fclose(PTR_BALANCE);
      return;
 }
 
} /******************************************************************************/

/******************************************************************************/

void writeIC( double y[],double t )

/******************************************************************************/

{
  PTR_IC = fopen("IC.bin", "wb");
  fwrite(&t,sizeof(double), 1, PTR_IC);
  fwrite(y,sizeof(double), NEQN, PTR_IC);
  fwrite(&seed,sizeof(double), 1, PTR_IC);
  printf("writing initial condition: \n");
  //printf("times is:  %f\n",t);
  int i;
  for(i=0;i<NEQN;i++)
  {
    //printf("%f\n",y[i]) ;
  } 
  fclose(PTR_IC); 
}
/******************************************************************************/

/******************************************************************************/

double readIC( double y[],double t )

/******************************************************************************/

{
  PTR_IC = fopen("IC.bin", "rb");
  fread(&t,sizeof(double), 1, PTR_IC);
  fread(y,sizeof(double), NEQN,PTR_IC);
  fread(&seed,sizeof(int),1,PTR_IC);
  printf("reading initial condition: \n");
  //printf("times is:  %f\n",t);
  int i;
  for(i=0;i<NEQN;i++)
  {
    //printf("%f\n",y[i]) ;
  }
  fclose(PTR_IC);
  return t;
}
/******************************************************************************/

/******************************************************************************/

void euler(double y[], double yp[], double t, int neqn, double dt)

/******************************************************************************/

{

  astro_d(t,y,yp );
  int i;

  for (i = 0; i < neqn; i++ ){
  
    y[i] = y[i] + dt*yp[i] ;

	if(y[i]>10000)
	{
		printf("\nIn Euler\n");
		printf("I_step is  %d\n", I_step);
		printf("y[%d] is %f\n", i, y[i]);
		printf("yp[%d] is %f", i, yp[i]);
	//	exit(0);
	}
  }
  return ;

}

/******************************************************************************/

/******************************************************************************/

double f(double ca,double cer) 

/******************************************************************************/
{
  double c1 = .13;
  double c2= 1.0;
  double c3 = .004;
  
  return c1* pow(ca,2)/(1+pow(ca,2)) - (pow(cer,2)/(1+pow(cer,2)))*(pow(ca,4)/(pow(c2,4)+pow(ca,4))) - c3*cer;
}

double fserca(double ca,double cer) 
{
  double c1 = .13;
  return c1* pow(ca,2)/(1+pow(ca,2));
}
double fcicr(double ca,double cer) 
{
  double c2= 1.0;
  return -(pow(cer,2)/(1+pow(cer,2)))*(pow(ca,4)/(pow(c2,4)+pow(ca,4)));
}
double flinear(double ca, double cer)
{
  double c3 = .004;
  return -c3*cer;
}

/******************************************************************************/

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

  
  int i;
  //printf("t: %f,(r,neuro_input, dr): ", t);
  for(i=0;i<NPROC;i++) {
  	 //if (i == 2) { r = 0.;} else { r = 0.31; }
  	 //r = 0.;  // there should be no soma spikes
  	 yp[2*i] = (neuro_input[i] * dr[i]/sqrt(dt) + r - c4 * f(y[2*i],y[2*i+1]) - y[2*i] + dc*(y[NEQN-2] - y[2*i])) / tau  ; //process i cytosol
     //printf("(%e,%e,%e), ", r, neuro_input[i], dr[i]);
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
  //printf("\n");

	
  diffc[NGROUPS-1] = cytosol_diff;
  differ[NGROUPS-1] = der_diff;
  serca_current[NGROUPS-1]  = fserca(y[NEQN-2],y[NEQN-1]);
  cicr_current[NGROUPS-1]   = fcicr(y[NEQN-2],y[NEQN-1]);
  linear_current[NGROUPS-1] = flinear(y[NEQN-2],y[NEQN-1]);
  ca_removal[NGROUPS-1] = -y[NEQN-2];

  // GE: May 1, 2018: turn soma off past a certain time
  double rsoma = r; // 0.31 by default

  //if (t > t_rsoma_turnoff) rsoma = 0.;
  //printf("t= %f\n", t);
  //printf("rsoma= %f\n", rsoma);exit(0);
  //printf("t= %f, rsoma= %f\n", t, rsoma);

  // Turn off soma (cytool and CR) through its currents when Ca in the cytosol is above a certain threshold
  double f_factor = 1.;
  double soma_thresh = 0.3;

  if (y[NEQN-2] > soma_thresh) {
      f_factor = 0.;
  } else {
      f_factor = 1.;
  }
  printf("t= %f, f_factor= %f, soma_cytosol= %f\n", t, f_factor, yp[NEQN-2]);

#ifdef VD
  yp[NEQN-2] = ( (rsoma -c4 * f_factor * f(y[NEQN-2],y[NEQN-1]) - y[NEQN-2]) * Vd_surf + cytosol_diff) / (Vd*tau); //soma cytosol
  yp[NEQN-1] = (f_factor * f(y[NEQN-2],y[NEQN-1]) * Vd_surf + der_diff) / (eps*tau*Vd); //soma er
  //float soma_clamp = 0.4;
  //if (yp[NEQN-2] > soma_clamp) yp[NEQN-2] = 0.2; //soma_clamp;
  //if (yp[NEQN-2] > soma_clamp) rsoma = 0.0;
#else
  yp[NEQN-2] = (rsoma - c4 * f(y[NEQN-2],y[NEQN-1]) - y[NEQN-2] + cytosol_diff) / tau; //soma cytosol
  yp[NEQN-1] = (f(y[NEQN-2],y[NEQN-1]) + der_diff) / (eps*tau); //soma er
  printf("  not computing with Vd though you should be   ");
#endif




  return;
}

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

