#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MSOFT.h"


int main(int argc, char **argv) {

  int step; /* Simulation loop iteration index */

  FILE *f2=fopen("psisq.dat","w");
  FILE *f3=fopen("norm.dat","w");
  FILE *f4=fopen("potential_ad.dat","w");
  FILE *f5=fopen("population.dat","w");
  FILE *f6=fopen("potential_di.dat","w");
  //FILE *f7=fopen("test.dat","w");

  init_param();  /* Initialize input parameters */
  init_prop();   /* Initialize the kinetic & potential propagators */
  init_wavefn(); /* Initialize the electron wave function */
  generate_trajectory(); /* Generate the trajectory for Tully Surface hopping*/
  pop_tsh_state();
  print_pot_ad(f4);
  print_pot_di(f6);
  print_wavefn(0,f2,f3); /* print initial conditions */

  for (step=1; step<=NSTEP; step++) {
    single_step(step); /* Time propagation for one step, DT */
    tsh_single_step();
    pop_tsh_state();
    if (step%NECAL==0) {
      pop_states();
      print_pop(step,f5);
    }
    if (step%NNCAL==0) {
      calc_norm();
      print_wavefn(step,f2,f3);
    }
  }

  return 0;
}
/*------------------------------------------------------------------------------*/
void init_param() {
/*------------------------------------------------------------------------------
  Initializes parameters by reading them from standard input.
------------------------------------------------------------------------------*/
  /* Initialize control parameters */
  LX=50.0;
  DT=0.02;
  NECAL=1;
  NNCAL=100;

  /* Calculate the mesh size */
  dx = LX/NX;
}

/*----------------------------------------------------------------------------*/
void init_prop() {
/*------------------------------------------------------------------------------
  Initializes the kinetic & potential propagators.
------------------------------------------------------------------------------*/
  int i;
  double x, k;

  A=0.01;
  B=1.6;
  C=0.005;
  D=1.0;
  M=2000;

  D1=0.02278; /*Parameter for Morse potential*/
  B1=0.675;
  b1=2.0;
  E1=0;
  D2=0.01025;
  B2=0.953;
  b2=3.212;
  E2=0.0038;
  A12=0.006337;
  b12=0.56;
  Rx=2.744;

  /* Set up kinetic propagators */
  for (i=0; i<=NX; i++) {
    if (i < NX/2)
      k = 2*M_PI*i/LX;
    else
      k = 2*M_PI*(i-NX)/LX;

    /* kinetic operator */
    T[i] = k*k*0.5/M;
    /* kinetic propagator */
    t[i][0] = cos(-DT*T[i]);
    t[i][1] = sin(-DT*T[i]);
  }
  dp= 2*M_PI/LX;
  /* Set up potential propagator */
  for (i=1; i<=NX; i++) {
    x = -0.5*LX+ dx*i;
    /* Construct the matrix h */
    /*Tully 2 potential */
     h[i][0][0] = A*(1-exp(-B*fabs(x)))*x/fabs(x);
     h[i][0][1] = C*exp(-D*x*x);
     h[i][1][0] = h[i][0][1];
     h[i][1][1] = -h[i][0][0];
     if(i==NX/2) {
       h[i][0][0] = 0;
       h[i][0][1] = C;
       h[i][1][0] = C;
       h[i][1][1] = 0;
     }
    /*Morse potential*/
    /*h[i][0][0]=D1*(1-exp(-B1*(x-b1)))*(1-exp(-B1*(x-b1)))+E1;
    h[i][1][1]=D2*(1-exp(-B2*(x-b2)))*(1-exp(-B2*(x-b2)))+E2;
    h[i][1][0]=A12*exp(-b12*(x-Rx)*(x-Rx));
    h[i][0][1]=h[i][1][0];*/

    if (h[i][1][1]>h[i][0][0]){
      intersect=i;
    }


    /* calc eigenvalues of h */
    calc_eigenvalues(i);
    /* calc De and Det */
    calc_De_and_Det(i);
    /* Half-step diagonal propagator */
    /* 1st component */
    u[i][0][0] = cos(-0.5*DT*E[i][0]);
    u[i][0][1] = sin(-0.5*DT*E[i][0]);

    /* 2nd component */
    u[i][1][0] = cos(-0.5*DT*E[i][1]);
    u[i][1][1] = sin(-0.5*DT*E[i][1]);
  }
}

/*----------------------------------------------------------------------------*/
void init_wavefn() {
/*------------------------------------------------------------------------------
  Initializes the components of the wave function.
------------------------------------------------------------------------------*/
  int sx,s;
  double x,gauss,Csq,norm_fac;
  /*Parameter for Tully surface in comment*/
  X0=0; /*-5.8*/
  S0=2.5;    /*1*/
  P0=0;   /*55*/ //sb: 55

  /* Calculate the the wave function value mesh point-by-point */
  for (sx=1; sx<=NX; sx++) {
    x = -0.5*LX+ dx*sx;
    gauss = exp(-S0*(x-X0)*(x-X0));
    C1[sx][0] = gauss*cos(P0*(x-X0)); 	/* wf on surface 1 */
    C1[sx][1] = gauss*sin(P0*(x-X0));
    C2[sx][0] = 0;			/* wf on surface 2 */
    C2[sx][1] = 0;
  }//[0] represent the real part and [i] represent the imaginary part e(it)=cos(t)+isin(t)

  /* Normalize C1 */
  Csq=0.0;
  for (sx=1; sx<=NX; sx++)
    for (s=0; s<2; s++)
      Csq += C1[sx][s]*C1[sx][s]; //multiply
  Csq *= dx;
  norm_fac = 1.0/sqrt(Csq);
  for (sx=1; sx<=NX; sx++)
    for (s=0; s<2; s++)
      C1[sx][s] *= norm_fac;

  /* Normalize C2 */
/*    Csq=0.0;
    for (sx=1; sx<=NX; sx++)
    for (s=0; s<2; s++)
    Csq += C2[sx][s]*C2[sx][s];
    Csq *= dx;
    norm_fac = 1.0/sqrt(Csq);
    for (sx=1; sx<=NX; sx++)
      for (s=0; s<2; s++)
        C2[sx][s] *= norm_fac; */

  periodic_bc();
}

/*----------------------------------------------------------------------------*/
void generate_trajectory(){
    /*Generate trajectory according to a Gaussian distribution to mimick the shape of wavepacket*/
    int nb,sx,pop_x,nb_index;
    double x;
    double pop_size=1/nb_traj;
    double p0=55;
    nb=0;
    for (sx=1;sx<=NX;sx++){ /*May need some improvement*/
      x = -0.5*LX+ dx*sx;
      pop_x=(int)((C1[sx][0]*C1[sx][0]+C1[sx][1]+C1[sx][1])*dx);
      if ((int)(pop_x/pop_size)>1) {
          for (nb_index = nb; nb_index <= (nb + pop_x / pop_size); nb_index++) {
              traj[nb_index][1] = x;
              traj[nb_index][2]= p0;
              traj[nb_index][3]=0;
          }
          nb=nb_index+1;
      }
    }

}
/*----------------------------------------------------------------------------*/
void periodic_bc() {
/*------------------------------------------------------------------------------
  Applies the periodic boundary condition to wave function, by copying
  the boundary values to the auxiliary array positions at the other ends.
------------------------------------------------------------------------------*/
  int s;

  /* Copy boundary wave function values */
  for (s=0; s<=1; s++) {
    C1[0][s] = C1[NX][s];
    C1[NX+1][s] = C1[1][s];
    C2[0][s] = C2[NX][s];
    C2[NX+1][s] = C2[1][s];
  }
}

/*----------------------------------------------------------------------------*/
void four1(double data[], unsigned long nn, int isign)
/*******************************************************************************
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as
1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
if isign is input as -1.  data is a complex array of length nn or, equivalently,
a real array of length 2*nn.  nn MUST be an integer power of 2 (this is not
checked for!).
*******************************************************************************/
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) { /* This is the bit-reversal section of the routine. */
    if (j > i) {
      SWAP(data[j],data[i]); /* Exchange the two complex numbers. */
      SWAP(data[j+1],data[i+1]);
    }
    m=nn;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  mmax=2;
  while (n > mmax) { /* Outer loop executed log2 nn times. */
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax); /* Initialize the trigonometric recurrence. */
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) { /* Here are the two nested inner loops. */
      for (i=m;i<=n;i+=istep) {
        j=i+mmax; /* This is the Danielson-Lanczos formula. */
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
/*----------------------------------------------------------------------------*/
void create_C1f() {
  int sx;
  for (sx=1; sx <= NX+1; sx++) {
    Cf[2*sx-1] = C1[sx-1][0];
    Cf[2*sx] = C1[sx-1][1];
  }
}
/*----------------------------------------------------------------------------*/
void create_C2f() {
  int sx;
  for (sx=1; sx <= NX+1; sx++) {
    Cf[2*sx-1] = C2[sx-1][0];
    Cf[2*sx] = C2[sx-1][1];
  }
}
/*----------------------------------------------------------------------------*/
void update_C1() {
  int sx;
  for (sx=1; sx <= NX+1; sx++) {
    C1[sx-1][0] = Cf[2*sx-1];
    C1[sx-1][1] = Cf[2*sx];
  }
}
/*----------------------------------------------------------------------------*/
void update_C2() {
  int sx;
  for (sx=1; sx <= NX+1; sx++) {
    C2[sx-1][0] = Cf[2*sx-1];
    C2[sx-1][1] = Cf[2*sx];
  }
}
/*----------------------------------------------------------------------------*/
void calc_eigenvalues(int i) {
  //E1 = (h11+h22+sqrt((h11-h22)*(h11-h22)+4*h12*h21))/2;  //nb1//energy obtained by diagonalization of matrix => eigenvalue problem
  //E2 = (h11+h22-sqrt((h11-h22)*(h11-h22)+4*h12*h21))/2;
  E[i][0] = (h[i][0][0]+h[i][1][1]+sqrt((h[i][0][0]-h[i][1][1])*(h[i][0][0]-h[i][1][1])+4*h[i][0][1]*h[i][1][0]))/2;
  E[i][1] = (h[i][0][0]+h[i][1][1]-sqrt((h[i][0][0]-h[i][1][1])*(h[i][0][0]-h[i][1][1])+4*h[i][0][1]*h[i][1][0]))/2;
}
/*----------------------------------------------------------------------------*/
void calc_De_and_Det(int i) {
  double norm_fac,der_1,der_2,der_3,der_4,xpos;
  /* De */
  if(i<=intersect) {
    /* 1st column = 1st eigenvector */
    De[i][0][0] = h[i][0][1]/(E[i][0]-h[i][0][0]);
    De[i][1][0] = 1;
    norm_fac = sqrt(De[i][0][0]*De[i][0][0]+1);
    De[i][0][0] /= norm_fac;
    De[i][1][0] /= norm_fac;
    /* 2nd column = 2nd eigenvector */
    De[i][0][1] = -De[i][1][0];
    De[i][1][1] = De[i][0][0];
  }
  else {
    /* 2nd column = 2nd eigenvector */
    De[i][0][1] = h[i][0][1]/(E[i][1]-h[i][0][0]);
    De[i][1][1] = 1;
    norm_fac = sqrt(De[i][0][1]*De[i][0][1]+1);
    De[i][0][1] /= norm_fac;
    De[i][1][1] /= norm_fac;
    /* 1st column = 1st eigenvector */
    De[i][1][0] = -De[i][0][1];
    De[i][0][0] = De[i][1][1];
  }
  Det[i][0][0] = De[i][0][0];
  Det[i][1][1] = De[i][1][1];
  Det[i][1][0] = De[i][0][1];
  Det[i][0][1] = De[i][1][0];

    if (i>0){
        der_1=(De[i][0][0]-De[i-1][0][0])/dx;
        der_2=(De[i][1][0]-De[i-1][1][0])/dx;
        der_3=(De[i][0][1]-De[i-1][0][1])/dx;
        der_4=(De[i][1][1]-De[i-1][1][1])/dx;
        d12[i]=De[i][0][0]*der_3+De[i][1][0]*der_4;
        xpos=-0.5*LX+i*dx;
        printf("%15.10f %15.10f\n",xpos,d12[i]);
    }
    if (i==1){
        d12[0]=d12[1];
    }


}

/*----------------------------------------------------------------------------*/
void pot_prop() {
/*------------------------------------------------------------------------------
  Potential propagator for a half time step, DT/2.
------------------------------------------------------------------------------*/
  int sx;
  double wr1,wi1,wr2,wi2;

  for (sx=1; sx<=NX; sx++) {
    /* De<r|psi> */
    wr1=De[sx][0][0]*C1[sx][0]+De[sx][0][1]*C2[sx][0];
    wi1=De[sx][0][0]*C1[sx][1]+De[sx][0][1]*C2[sx][1];
    wr2=De[sx][1][0]*C1[sx][0]+De[sx][1][1]*C2[sx][0];
    wi2=De[sx][1][0]*C1[sx][1]+De[sx][1][1]*C2[sx][1];
    C1[sx][0]=wr1;
    C1[sx][1]=wi1;
    C2[sx][0]=wr2;
    C2[sx][1]=wi2;
    /* exp(-iE1dt)*C1, exp(-iE2dt)*C2 */
    /* 1st component */
    wr1=u[sx][0][0]*C1[sx][0]-u[sx][0][1]*C1[sx][1];
    wi1=u[sx][0][0]*C1[sx][1]+u[sx][0][1]*C1[sx][0];
    /* 2nd component */
    wr2=u[sx][1][0]*C2[sx][0]-u[sx][1][1]*C2[sx][1];
    wi2=u[sx][1][0]*C2[sx][1]+u[sx][1][1]*C2[sx][0];

    C1[sx][0]=wr1;
    C1[sx][1]=wi1;
    C2[sx][0]=wr2;
    C2[sx][1]=wi2;
    /* Det<r|fi> */
    wr1=Det[sx][0][0]*C1[sx][0]+Det[sx][0][1]*C2[sx][0];
    wi1=Det[sx][0][0]*C1[sx][1]+Det[sx][0][1]*C2[sx][1];
    wr2=Det[sx][1][0]*C1[sx][0]+Det[sx][1][1]*C2[sx][0];
    wi2=Det[sx][1][0]*C1[sx][1]+Det[sx][1][1]*C2[sx][1];
    C1[sx][0]=wr1;
    C1[sx][1]=wi1;
    C2[sx][0]=wr2;
    C2[sx][1]=wi2;
  }
  periodic_bc();
}
/*----------------------------------------------------------------------------*/
void kin_prop() {
/*------------------------------------------------------------------------------
  Kinetic propagation for t step.
-------------------------------------------------------------------------------*/
  int sx,s;
  double wr,wi;

  for (sx=1; sx<=NX; sx++) {
    /* 1st component C1 */
    wr=t[sx][0]*C1[sx][0]-t[sx][1]*C1[sx][1];
    wi=t[sx][0]*C1[sx][1]+t[sx][1]*C1[sx][0];
    C1[sx][0]=wr;
    C1[sx][1]=wi;
    /* 2nd component C2 */
    wr=t[sx][0]*C2[sx][0]-t[sx][1]*C2[sx][1];
    wi=t[sx][0]*C2[sx][1]+t[sx][1]*C2[sx][0];
    C2[sx][0]=wr;
    C2[sx][1]=wi;
  }
}
/*----------------------------------------------------------------------------*/
void single_step(int step) {
/*------------------------------------------------------------------------------
  Propagates the wave function for a unit time step, DT.
------------------------------------------------------------------------------*/
  int j, s, i, ibis;
  /* half step potential propagation */
  pot_prop();
  /* fft of the 2 components of <r|psi> */
  /* 1st component */
  create_C1f();
  four1(Cf, (unsigned long) NX, -1);
  for (j=0; j <= 2*(NX+1); j++)
    Cf[j] /= NX;
  update_C1();
  /* 2nd component */
  create_C2f();
  four1(Cf, (unsigned long) NX, -1);
  for (j=0; j <= 2*(NX+1); j++)
    Cf[j] /= NX;
  update_C2();
  /* step kinetic propagation   */
  kin_prop();
  /* fft^(-1) */
  /* 1st component */
  create_C1f();
  four1(Cf, (unsigned long) NX, 1);
  update_C1();
  /* 2nd component */
  create_C2f();
  four1(Cf, (unsigned long) NX, 1);
  update_C2();
  /* half step potential propagation */
  pot_prop();
}
/*----------------------------------------------------------------------------*/
void tsh_single_step(){
/*------------------------------------------------------------------------------
  Propagates the trajectory for a unit time step, DT.
------------------------------------------------------------------------------*/
    int nb,surf_1,surf_2,pos;
    double f_t,f_tdt,hop_prob,random_nb,new_energy;
    /*Computation of the classical dynamic*/
    for (nb=0;nb<nb_traj;nb++){
        /*Verley-Velocity algorithm or other scheme (Runge-Kutta-Gill)*/
        f_t=0 ;/*compute the force at position x(t) => Use Hellmann-Feynmann theorem to compute force*/
        traj[nb][1]=traj[nb][1]+(traj[nb][2]*DT/(M))+(f_t*DT*DT/(2*M)); /*Compute the new position*/
        f_tdt=0 ; /*compute the force at position x(t+dt)*/
        traj[nb][2]=traj[nb][2]+((f_tdt+f_t)*DT/(2)); /*Compute the new momentum*/
    }
    /*Compute the evolution of the quantum amplitude A */

    /*Compute the amplitude transfer by hop*/
    for (surf_1=0;surf_1<=1;surf_1++) {
        for(surf_2=0;surf_2<=1;surf_2++){
            /*Compute the hop probability*/
            if(surf_1!=surf_2) {
                hop_prob =0.5 ;
                for (nb = 0; nb < nb_traj; nb++) {
                    if (traj[nb][3] == surf_1) {
                        /*Generate a random number between 0 and 1*/
                        random_nb = (double)rand()/RAND_MAX;
                        /*Transfer of population if condition satisfied*/
                        if(random_nb<=hop_prob){
                            pos=(int)(traj[nb][1]+0.2*LX/dx);/*Compute the index of position*/
                            new_energy=(traj[nb][2]*traj[nb][2]/(2*M))+E[pos][surf_1]-E[pos][surf_2]; /*!!!Need a check for energy conservation*/
                            if (new_energy>0){ /*condition to ignore frustrated hop*/
                                traj[nb][3]=surf_2;/*Transfer of population*/
                                traj[nb][2]=sqrt(new_energy*2*M);/*Change of momentum for energy conservation*/
                            }
                        }
                    }
                }
            }
        }
    }



}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*Differential equation*/
/*double coupled_channel_R(double t, double c,int i,int pos){
    double result=E[pos][i]*c;
    return result;
}
double coupled_channel_Im(double t, double c,int i,int p){
    double result=-p*d_ij;
    return result;
}

double newton_equation(){

}*/
/*----------------------------------------------------------------------------*/
void pop_states() {
  int i;
  P1=0.0;
  P2=0.0;
  for (i=1; i<=NX;i++) {
    P1 += C1[i][0]*C1[i][0]+C1[i][1]*C1[i][1];
    P2 += C2[i][0]*C2[i][0]+C2[i][1]*C2[i][1];
  }
  P1 *= dx;
  P2 *= dx;
}
/*----------------------------------------------------------------------------*/
void pop_tsh_state(){
    double pop_1,pop_2;
    int nb;

    pop_1=0;
    pop_2=0;
    for (nb=0;nb<nb_traj;nb++){
        if (traj[nb][3]==0){
            pop_1++;
        }
        else{
            pop_2++;
        }
    }
    //printf("%15.10f %15.10f\n",pop_1/nb_traj,pop_2/nb_traj);
}
/*----------------------------------------------------------------------------*/

void print_pop(int step, FILE *f5) {
  int i;

  fprintf(f5,"%8i %15.10f %15.10f\n",step,P1,P2);
}

/*-------------------------------------------------------------------------------
  Print the potential
-------------------------------------------------------------------------------*/

void print_pot_ad(FILE *f4) {

  int i;
  double x;

  for (i=1; i<=NX; i++) {
    x = dx*i;
    fprintf(f4,"%8i %15.10f %15.10f %15.10f\n",i,x,E[i][1],E[i][0]);  //nb1
  }
}
/*-------------------------------------------------------------------------------*/
void print_pot_di(FILE *f6) {

  int i;
  double x;

  for (i=1; i<=NX; i++) {
    x = dx*i;
    fprintf(f6,"%8i %15.10f %15.10f %15.10f %15.10f %15.10f\n"
            ,i,x,h[i][0][0],h[i][1][1],h[i][0][1],h[i][1][0]);
  }
}
/*-------------------------------------------------------------------------------*/
void calc_norm() {
/*------------------------------------------------------------------------------
  Calculate the norm
-------------------------------------------------------------------------------*/

  int sx;
  double psisq,psisq2;

  norm=0.0;

  for (sx=1; sx<=NX; sx++)                              // domod from 1 to NX (not NX+1)
  {
    psisq=C1[sx][0]*C1[sx][0]+C1[sx][1]*C1[sx][1]+
          C2[sx][0]*C2[sx][0]+C2[sx][1]*C2[sx][1];
    psisq2=C1[sx+1][0]*C1[sx+1][0]+C1[sx+1][1]*C1[sx+1][1]+
           C2[sx+1][0]*C2[sx+1][0]+C2[sx+1][1]*C2[sx+1][1];
    norm=norm+dx*((psisq2+psisq)/2.0);
  }
  norm=sqrt(norm);
}

/*-------------------------------------------------------------------------------*/
void print_wavefn(int step, FILE *f2, FILE *f3) {

  int sx;
  double x;

  fprintf(f2,"\n");
  fprintf(f2,"\n");
  for (sx=1; sx<=NX; sx++)
  {
    x=dx*sx;
    fprintf(f2,"%8i %15.10f %15.10f %15.10f\n",sx,x,
            ((C1[sx][0]*C1[sx][0]+C1[sx][1]*C1[sx][1])/100), // "/100" for visualization purpose
            ((C2[sx][0]*C2[sx][0]+C2[sx][1]*C2[sx][1])/100));  //nb1
  }

  if (step>0)
  {
    fprintf(f3,"%8i %15.10f\n",step,norm);
  }
}
/*----------------------------------------------------------------------------*/
void print_avg(int step, FILE *f7) {
  fprintf(f7, "%8i %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n",step,DT*step,X1[step],X2[step],X_avg[step],MP_1[step],MP_2[step],MP_avg[step]);
}
/*----------------------------------------------------------------------------*/
