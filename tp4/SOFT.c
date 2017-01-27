/*******************************************************************************
Quantum dynamics (QD) simulation, SOFT algorithm.
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SOFT.h"

#define FOLDER "C:\\Users\\Ni\\Desktop\\Uni\\9.Semester\\Computational Methods in Molecular Quantum Mechanics\\tp4\\output\\tun1_"

int main(int argc, char **argv) {

    int step; /* Simulation loop iteration index */

    FILE *f1 = fopen(FOLDER "energy.dat", "w");
    FILE *f2 = fopen(FOLDER "psi.dat", "w");
    FILE *f3 = fopen(FOLDER "psisq.dat", "w");
    FILE *f4 = fopen(FOLDER "norm.dat", "w");
    FILE *f5 = fopen(FOLDER "potential.dat", "w");

    init_param();      /* Initialize input parameters */
    init_prop_tun();   /* Initialize the kinetic & potential propagators (choice: "box" or "mol")*/
    init_wavefn(f2, f3, f4);     /* Initialize the electron wave function */

    print_pot(f5);
    print_wavefn(0, f2, f3, f4); /* Daniele: starting wf */

    for (step = 1; step <= NSTEP; step++) {
        single_step(step); /* Time propagation for one step, DT */
        if (step % NECAL == 0) {
            calc_energy();
            print_energy(step, f1);
        }
        if (step % NNCAL == 0) {
            calc_norm();
            print_wavefn(step, f2, f3, f4);
        }
        if(step % 1000 == 0){
            printf("%d %\n",step/1000);
        }
    }

    return 0;
}

/*------------------------------------------------------------------------------*/
void init_param() {
/*------------------------------------------------------------------------------
  Initializes parameters.
------------------------------------------------------------------------------*/
    /* initialize control parameters */
    LX = 50.0;
    DT = 0.001;
    NSTEP = 300000;//= 100000;
    NECAL = 1;    /* print energy */
    NNCAL = 100;  /* print psi, psisq, norm */

    /* Calculate the mesh size */
    dx = LX / NX;
}

/*----------------------------------------------------------------------------*/
void init_prop_box() {
/*------------------------------------------------------------------------------
  Initializes the kinetic & potential propagators for an electron in a box.
------------------------------------------------------------------------------*/
    int i;
    double x, k;

    BH = 0.8;   /* height of central barrier */
    BW = 1.0;   /* width  of central barrier */
    EH = 50.0;  /* height of edge barrier    */

    X0 = 10;    /* starting pos                      */
    E0 = 9.0;   /* height of the initial gaussian WF */
    S0 = 3.0;   /* width  of the initial gaussian WF */

    /* Set up kinetic propagators */
    for (i = 0; i <= NX; i++) {
        if (i < NX / 2)
            k = 2 * M_PI * i / LX;
        else
            k = 2 * M_PI * (i - NX) / LX;

        /* kinetic operator */
        T[i] = 0.5 * k * k;
        /* kinetic propagator */
        t[i][0] = cos(-DT * T[i]);
        t[i][1] = sin(-DT * T[i]);
    }

    /* Set up potential propagator */
    for (i = 1; i <= NX; i++) {
        x = dx * i;
        /* Construct the edge potential */
        if (i == 1 || i == NX)
            v[i] = EH;
            /* Construct the barrier potential */
        else if (0.5 * (LX - BW) < x && x < 0.5 * (LX + BW))
            v[i] = BH;
        else
            v[i] = 0.0;
        /* Half-step potential propagator */
        u[i][0] = cos(-0.5 * DT * v[i]);
        u[i][1] = sin(-0.5 * DT * v[i]);
    }
}


/*----------------------------------------------------------------------------*/
void init_prop_tun() {
/*------------------------------------------------------------------------------
  Initializes the kinetic & potential propagators for an electron in a box.
------------------------------------------------------------------------------*/
    int i;
    double x, k;

    BH = 8;   /* height of central barrier */
    BW = 1.0;   /* width  of central barrier */
    EH = 50.0;  /* height of edge barrier    */

    X0 = 10;    /* starting pos                      */
    E0 = 4;   /* height of the initial gaussian WF */
    S0 = 3.0;   /* width  of the initial gaussian WF */
    double MASS = 1;

    /* Set up kinetic propagators */
    for (i = 0; i <= NX; i++) {
        if (i < NX / 2)
            k = 2 * M_PI * i / LX;
        else
            k = 2 * M_PI * (i - NX) / LX;

        /* kinetic operator */
        T[i] = 0.5 * k * k / MASS;
        /* kinetic propagator */
        t[i][0] = cos(-DT * T[i]);
        t[i][1] = sin(-DT * T[i]);
    }

    /* Set up potential propagator */
    for (i = 1; i <= NX; i++) {
        x = dx * i;
        /* Construct the edge potential */
        if (i == 1 || i == NX)
            v[i] = EH;
            /* Construct the barrier potential */
        else if (0.5 * (LX - BW) < x && x < 0.5 * (LX + BW))
            v[i] = BH;
        else
            v[i] = 0.0;
        /* Half-step potential propagator */
        u[i][0] = cos(-0.5 * DT * v[i]);
        u[i][1] = sin(-0.5 * DT * v[i]);
    }
}

/*----------------------------------------------------------------------------*/
void init_prop_mol() {
/*------------------------------------------------------------------------------
  Initializes the kinetic & potential propagators for a morse potential.
------------------------------------------------------------------------------*/
    int i;
    double x, k;
    double b, D;
    b = 0.3;
    D = 1.0;

    X0 = 15.0;
    E0 = 0.0;
    S0 = 1.0;

    /* Set up kinetic propagators */
    for (i = 0; i <= NX; i++) {
        if (i < NX / 2)
            k = 2 * M_PI * i / LX;
        else
            k = 2 * M_PI * (i - NX) / LX;

        /* kinetic operator */
        T[i] = 0.5 * k * k;
        /* kinetic propagator */
        t[i][0] = cos(-DT * T[i]);
        t[i][1] = sin(-DT * T[i]);
    }

    /* Set up potential propagator */
    for (i = 1; i <= NX; i++) {
        x = dx * i;
        /* potential operator */
        v[i] = D * (exp(-2 * b * (x - 10)) - 2 * exp(-b * (x - 10))) + D;
        /* Half-step potential propagator */
        u[i][0] = cos(-0.5 * DT * v[i]);
        u[i][1] = sin(-0.5 * DT * v[i]);
    }
}

/*----------------------------------------------------------------------------*/
void init_wavefn() {
/*------------------------------------------------------------------------------
  Initializes the wave function as a traveling Gaussian wave packet.
------------------------------------------------------------------------------*/
    int sx, s;
    double x, gauss, psisq, norm_fac;


    /* Calculate the the wave function value mesh point-by-point */
    for (sx = 1; sx <= NX; sx++) {
        x = dx * sx - X0;
        gauss = exp(-0.25 * x * x / (S0 * S0));
        psi[sx][0] = gauss * cos(-sqrt(2.0 * E0) * x);
        psi[sx][1] = gauss * sin(-sqrt(2.0 * E0) * x);
    }

    /* Normalize the wave function */
    psisq = 0.0;
    for (sx = 1; sx <= NX; sx++)
        for (s = 0; s < 2; s++)
            psisq += psi[sx][s] * psi[sx][s];
    psisq *= dx;
    norm_fac = 1.0 / sqrt(psisq);
    for (sx = 1; sx <= NX; sx++)
        for (s = 0; s < 2; s++)
            psi[sx][s] *= norm_fac;
    periodic_bc();

}

/*----------------------------------------------------------------------------*/
void single_step(int step) {
/*------------------------------------------------------------------------------
  Propagates the electron wave function for a unit time step, DT.
------------------------------------------------------------------------------*/
    int j;

    pot_prop();                         /* half step potential propagation */
    create_psif();                      //psi -> psif
    four1(psif, (unsigned long) NX, -1);//fast fourier transform -> result in psif (-1 -> forward)
    for (j = 0; j <= 2 * (NX + 1); j++)
        psif[j] /= NX;                  //normalize
    update_psi();                       //psif -> psi
    kin_prop();                         /* step kinetic propagation   */
    if (step % NECAL == 0)
        calc_ekin();                    //psi has to be in p-basis!
    create_psif();                      //psi -> psif
    four1(psif, (unsigned long) NX, 1); //fft (this time backwards)
    update_psi();                       //psif -> psi
    pot_prop();                         /* half step potential propagation */
    if (step % NECAL == 0)
        calc_epot();                    //psi has to be in x-basis!
}

/*----------------------------------------------------------------------------*/
void pot_prop() {
/*------------------------------------------------------------------------------
  Potential propagator for a half time step, DT/2.
------------------------------------------------------------------------------*/
    int sx;
    double wr, wi;

    for (sx = 1; sx <= NX; sx++) {
        wr = u[sx][0] * psi[sx][0] - u[sx][1] * psi[sx][1];
        wi = u[sx][0] * psi[sx][1] + u[sx][1] * psi[sx][0];
        psi[sx][0] = wr;
        psi[sx][1] = wi;
    }
    periodic_bc();
}

/*----------------------------------------------------------------------------*/
void kin_prop() {
/*------------------------------------------------------------------------------
  Kinetic propagation for t step.
-------------------------------------------------------------------------------*/
    int sx, s;
    double wr, wi;

    for (sx = 1; sx <= NX; sx++) {
        wr = t[sx][0] * psi[sx][0] - t[sx][1] * psi[sx][1];
        wi = t[sx][0] * psi[sx][1] + t[sx][1] * psi[sx][0];
        psi[sx][0] = wr;
        psi[sx][1] = wi;
    }
}

/*----------------------------------------------------------------------------*/
void periodic_bc() {
/*------------------------------------------------------------------------------
  Applies the periodic boundary condition to wave function PSI, by copying
  the boundary values to the auxiliary array positions at the other ends.
------------------------------------------------------------------------------*/
    int s;

    /* Copy boundary wave function values */
    for (s = 0; s <= 1; s++) {
        psi[0][s] = psi[NX][s];
        psi[NX + 1][s] = psi[1][s];
    }
}

/*----------------------------------------------------------------------------*/
void calc_energy() {
/*------------------------------------------------------------------------------
  Calculates the total energy ETOT.
------------------------------------------------------------------------------*/
    /* Total energy */
    etot = ekin + epot;
}

/*----------------------------------------------------------------------------*/
void calc_epot() {
    int sx;
    /* Potential energy */
    epot = 0.0;
    for (sx = 1; sx <= NX; sx++) {
        epot += v[sx] * (psi[sx][0] * psi[sx][0] + psi[sx][1] * psi[sx][1]);
    }
    epot *= dx;
}

/*----------------------------------------------------------------------------*/
void calc_ekin() {
    int sx;
    double k;
    /* Kinetic energy */
    ekin = 0.0;
    for (sx = 1; sx <= NX; sx++) {
        if (sx < NX / 2)
            k = 2 * M_PI * sx / LX;
        else
            k = 2 * M_PI * (sx - NX) / LX;
        ekin += 0.5 * k * k * (psi[sx][0] * psi[sx][0] + psi[sx][1] * psi[sx][1]);
    }
    ekin *= dx;
    ekin *= NX;
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
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) { /* This is the bit-reversal section of the routine. */
        if (j > i) {
            SWAP(data[j], data[i]); /* Exchange the two complex numbers. */
            SWAP(data[j + 1], data[i + 1]);
        }
        m = nn;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    mmax = 2;
    while (n > mmax) { /* Outer loop executed log2 nn times. */
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax); /* Initialize the trigonometric recurrence. */
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) { /* Here are the two nested inner loops. */
            for (i = m; i <= n; i += istep) {
                j = i + mmax; /* This is the Danielson-Lanczos formula. */
                tempr = wr * data[j] - wi * data[j + 1];
                tempi = wr * data[j + 1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j + 1] = data[i + 1] - tempi;
                data[i] += tempr;
                data[i + 1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr; /* Trigonometric recurrence. */
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}

/*----------------------------------------------------------------------------*/
void create_psif() {
/*------------------------------------------------------------------------------
  Create an array for the fourier transform.
-------------------------------------------------------------------------------*/
    int sx;
    for (sx = 1; sx <= NX + 1; sx++) {
        psif[2 * sx - 1] = psi[sx - 1][0];
        psif[2 * sx] = psi[sx - 1][1];
    }
}

/*----------------------------------------------------------------------------*/
void update_psi() {
/*------------------------------------------------------------------------------
  Update psi with its fourier transform.
-------------------------------------------------------------------------------*/
    int sx;
    for (sx = 1; sx <= NX + 1; sx++) {
        psi[sx - 1][0] = psif[2 * sx - 1];
        psi[sx - 1][1] = psif[2 * sx];
    }
}

/*----------------------------------------------------------------------------*/
void print_energy(int step, FILE *f1) {
/*------------------------------------------------------------------------------
  Print energy values in file energy.dat
-------------------------------------------------------------------------------*/

    fprintf(f1, "%8i %15.10f %15.10f %15.10f\n", step, ekin, epot, etot);

}

/*----------------------------------------------------------------------------*/
void print_wavefn(int step, FILE *f2, FILE *f3, FILE *f4) {
/*------------------------------------------------------------------------------
  Print wf, squared wf and norm 
-------------------------------------------------------------------------------*/

    int sx;
    double x;

    fprintf(f2, "\n");
    fprintf(f2, "\n");
    fprintf(f3, "\n");
    fprintf(f3, "\n");
    for (sx = 1; sx <= NX; sx++) {
        x = dx * sx;
        fprintf(f2, "%8i %15.10f %15.10f %15.10f\n", sx, x, psi[sx][0], psi[sx][1]);
        fprintf(f3, "%8i %15.10f %15.10f\n", sx, x, psi[sx][0] * psi[sx][0] + psi[sx][1] * psi[sx][1]);
    }
    if (step > 0) {
        fprintf(f4, "%8i %15.10f\n", step, norm);
    }
}

/*----------------------------------------------------------------------------*/
void calc_norm() {
/*------------------------------------------------------------------------------
  Calculate the norm                        
-------------------------------------------------------------------------------*/

    int sx;
    double psisq, psisq2;

    norm = 0.0;

    for (sx = 1; sx <= NX + 1; sx++) {
        psisq = psi[sx][0] * psi[sx][0] + psi[sx][1] * psi[sx][1];
        psisq2 = psi[sx + 1][0] * psi[sx + 1][0] + psi[sx + 1][1] * psi[sx + 1][1];
        norm = norm + dx * ((psisq2 + psisq) / 2.0);
    }
    norm = sqrt(norm);
}

/*-------------------------------------------------------------------------------
  Print the potential
-------------------------------------------------------------------------------*/

void print_pot(FILE *f5) {

    int sx;
    double x;

    for (sx = 1; sx <= NX; sx++) {
        x = dx * sx;
        fprintf(f5, "%8i %15.10f %15.10f\n", sx, x, v[sx]);
    }
}
/*-------------------------------------------------------------------------------*/
