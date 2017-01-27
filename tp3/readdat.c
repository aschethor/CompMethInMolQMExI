#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include "system.h"

// read in system information
void Readdat(void) {
    FILE *FilePtr;

    FilePtr = fopen("C:\\Users\\Ni\\ClionProjects\\CompMethInMolQMExI\\run", "r");
    //FilePtr=fopen("input","r");
    fscanf(FilePtr, "%lf %d %d %lf %lf %d\n",
           &Box, &NumberOfParticles, &NumberOfSteps, &Temperature,
           &Deltat, &NumberOfInitializationSteps);

    /*
    while (!feof(FilePtr))printf("%c", fgetc(FilePtr));
    printf("-------------------\n");
    printf("%f \n %d \n %d \n %f \n %f \n %d",
           Box, NumberOfParticles, NumberOfSteps, Temperature, Deltat, NumberOfInitializationSteps);
    */

    fclose(FilePtr);

    //Because fscanf didn't work...
    Box = 5.0;
    NumberOfParticles = 100;    //variate that
    NumberOfSteps  = 50000;
    Temperature = 0.5;          //variate that...
    Deltat = 0.005;             //and that -> U_tot,U_kin,U_pot
    NumberOfInitializationSteps  = 5000;

    if (NumberOfParticles > MAXPART) {
        printf("Maximum number of particles is : %d\n", MAXPART);
        exit(1);
    }

    // Calculate Some Parameters
    CutOff = 0.49999 * Box;
    Ecut = 4.0 * (pow(SQR(CutOff), -6.0) - pow(SQR(CutOff), -3.0));

    // print information to the screen

    printf("Molecular Dynamics Program\n");
    printf("\n");
    printf("Number of particles   : %d\n", NumberOfParticles);
    printf("Boxlength             : %f\n", Box);
    printf("Density               : %f\n", NumberOfParticles / (Box * Box * Box));
    printf("Temperature           : %f\n", Temperature);
    printf("Cut-Off radius        : %f\n", CutOff);
    printf("Cut-Off energy        : %f\n", Ecut);
    printf("Number of uteps       : %d\n", NumberOfSteps);
    printf("Number of init steps  : %d\n", NumberOfInitializationSteps);
    printf("Timestep              : %f\n", Deltat);
}
