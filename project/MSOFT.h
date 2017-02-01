/*******************************************************************************
File MSOFT.h is a header file for program MSOFT.c.
*******************************************************************************/
#define NX 1024  /* Number of mesh points */
#define NSTEP 10000 /* number of simulation step*/
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define nb_traj 1000 /* Number of trajectory for TSH*/
#define n_t 10 /*scale number of dt with respect to DT dt=DT/n_t*/
#include "complex.h"


/* Function prototypes ********************************************************/
void init_param();
void init_prop();
void init_wavefn();
void single_step(int step);
void pot_prop();
void kin_prop();
void periodic_bc();
void four1(double data[], unsigned long nn, int isign);
void create_C1f();
void create_C2f();
void update_C1();
void update_C2();
void calc_norm();
void calc_eigenvalues(int i);
void calc_De_and_Det(int i);
void pop_states();
void print_pop(int step, FILE *f5);
void print_pot_ad(FILE *f4);
void print_pot_di(FILE *f6);
void print_wavefn(int step, FILE *f2, FILE *f3);
void print_avg(int step, FILE *f7);

/* Function prototype for Tully surface hopping*/
void tsh_single_step(int step);
void generate_trajectory();
void pop_tsh_state(FILE *f7,int step);
double rkg4(double (*f)(double, double, int), double y0, double x0,double h);
double energy_diabatic(double pos,int i_1,int i_2);
double eigenvalue_calc(double h_dia[2][2],double surf);
double calc_eigenvector(double x, int i_1, int i_2,double h_dia[2][2],double h_adia[2]);

double box_muller(double m, double s);

void test();

/* Input parameters ***********************************************************/
double LX;       /* Simulation box length */
double DT;       /* Time discretization unit */
double M;	 /* Mass of the system */
/*int NSTEP;*/       /* Number of simulation steps */
int NECAL;       /* Interval to calculate energies */
int NNCAL;       /* Interval to calculate norm     */
double X0,S0,P0; /* Center-of-mass, spread & momentum of initial wave packet */
double A,B,C,D;  /* Parameters of Tully potential */
double D1,B1,b1,E1,D2,B2,b2,E2,A12,b12,Rx; /*Parameters of Morse potential*/
double E0;
double intercept_cont;
/* Arrays **********************************************************************
C1[NX+2][2]:     C1[i][0|1] is the real|imaginary part of the first component of the
		 wave function on mesh point i
C2[NX+2][2]:     C2[i][0|1] is the real|imaginary part of the second component of the
		 wave function on mesh point i
Cf[(NX+1)*2]:    C[2i|2i+1] is the real|imaginary part of the wave function
		 on mesh point i
T[NX+2]:	 T[i] is the kinetic energy at mesh point i
t[NX+2][2]:	 t[i][] is the kinetic propagator on i (real|imaginary part)
u[NX+2][2][2]:   u[i][j][] is the jth component of the diagonal potential propagator on i (real|imaginary part)
h[NX+2][2][2]:   h11[i][j][k] is the element j,k of the matrix h at mesh point i
E[NX+2][2]:      E[i][j] is the j-th eigenvalue of the matrix at mesh point i
De[NX+2][2][2]:  De[i][j][k] is the element j,k of the matrix De, which has the
		 eigenvectors (of the matrix h) as columns, at mesh point i
Det[NX+2][2][2]: Det[i][j][k] is the element j,k of the matrix Det, which is the
		 transposed of De, at mesh point i

P_1,P_2,X1,X2      variable for calculating the average momentum and position of state 1 and 2
*******************************************************************************/
double C1[NX+2][2];
double C2[NX+2][2];
double Cf[(NX+2)*2];
double T[NX+2];
double t[NX+2][2];
double u[NX+2][2][2];
double h[NX+2][2][2];
double E[NX+2][2];
double De[NX+2][2][2];
double Det[NX+2][2][2];

double MP_1[NSTEP+1];
double MP_2[NSTEP+1];
double X1[NSTEP+1];
double X2[NSTEP+1];
double MP_avg[NSTEP+1];
double X_avg[NSTEP+1];

/*variable for trajectory
 * parameter 0=position
 * parameter 1=momentum
 * parameter 2=electronic surface
 * *****************************/
double traj[nb_traj][3];
double d12[NX+2];
double complex a[nb_traj][2][2];//for each trajectory [nb_traj] and for each state k,l [2][2]: complex number [2]
double complex b[nb_traj][2][2];//for each trajectory [nb_traj] and for each state k,l [2][2]: real number
/* Variables *******************************************************************
dx   = Mesh spacing
P1 = Population of state 1
P2 = Population of state 2
*******************************************************************************/
double dx;
double dp;
double norm;
double P1,P2,d1,d2;
int intersect;
