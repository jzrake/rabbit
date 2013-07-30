
/*------------------------------------------------------------------------------
 * FILE: euler.c
 *
 * Copyright (C) 2011 Jonathan Zrake, NYU CCPP: zrake <at> nyu <dot> edu
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.

 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 *
 *
 * Feel free to modify this code or use it for any research / educational
 * purposes. Just please use good manners when distributing it!
 *
 *
 *
 * DESCRIPTION:
 *
 * This code was written to demonstrate several numerical algorithms used in the
 * solution of the relativistic Euler equations. These equations describe the
 * motion of an inviscid fluid when its bulk or thermal velocity is near the
 * speed of light. All of the algorithms demonstrated here are of the
 * finite-volume type, which means that the Euler equations are solved in
 * conservative form, by explicitly passing conserved quantites (mass, momentum,
 * total energy) between neighboring volumes. These fluxes are computed from an
 * approximate solution of the Riemann problem posed by the discontinuous states
 * of adjacent zones. The algorithms demonstrated in this program are:
 *
 *         (1) The HLL approximate Riemann solver
 *         (2) Piecewise linear reconstruction (2nd order in space)
 *         (3) Runge-Kutta time integration (1st, 2nd, 3rd order in time)
 *         (4) Primitive variable recovery via Newton-Rapheson solver
 *
 * The code outputs an ASCII table named 'srhd.dat' containing the primitive
 * variables where the columns are
 *
 *                 [x, Density, Pressure, Vx, Vy, Vz]
 *
 * It also generates a minimal gnuplot script which visualizes the data, by
 * running
 *
 *                 $> gnuplot plot.gp
 *
 *
 * The code contained in this file was adapted largely from the Mara code, which
 * is a relativistic magnetohydrodynamic (RMHD) code written for the study of
 * astrophysical turbulence for the graduate work of J. Zrake under
 * A. MacFadyen. The routines contained here describe the full 3d form of the
 * Euler equations in cartesian coordinates. Therefore although the calling
 * functions assume a 1d problem, the code may easily be adapted to build a
 * fully 3d and parallel (MPI) simulation. Also, adding test problems and
 * integration algorithms should be straightforward as the code provides a
 * simple model for modularity. To add new test problems, place a function with
 * the following signature
 *
 *                  Primitive InitialCondition_xyz(double x);
 *
 * at the very bottom with the others, provide a prototype in the section below,
 * 'Test problem setup functions', and add it to the array TestProblems.
 *
 *
 * REFERENCES:
 *
 * The code is based largely on the 'How to Write a Hydro Code' document by
 * Weiqun Zhang
 *
 *            http://zrake.webfactional.com/static/notes/hydro_code.pdf
 *
 *
 *------------------------------------------------------------------------------
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "euler1d.h"

void SetBoundaryConditions(Conserved *U);
void TimeDerivative(const Conserved *U, Conserved *L);

double maxval3(double a, double b, double c);
double minval3(double a, double b, double c);
double sign(double x);
double minmod3(double ul, double u0, double ur);
char *ProgressBar(double frac);
void RunUserMenu();


// Test problem setup functions
// -----------------------------------------------------------------------------
Primitive InitialCondition_Shocktube1(double x);
Primitive InitialCondition_Shocktube2(double x);
Primitive InitialCondition_Blastwave1(double x);
Primitive InitialCondition_Blastwave2(double x);

struct TestProblem
{
  char Name[256];
  Primitive (*InitialCondition)(double x);
} TestProblems[5] =
  { {"Shocktube1", InitialCondition_Shocktube1 },
    {"Shocktube2", InitialCondition_Shocktube2 },
    {"Blastwave1", InitialCondition_Blastwave1 },
    {"Blastwave2", InitialCondition_Blastwave2 }, { "", NULL }
  } ;


// Time integration options
// -----------------------------------------------------------------------------
void AdvanceSolution_FwdEuler(Conserved *U, double dt);
void AdvanceSolution_Midpoint(Conserved *U, double dt);
void AdvanceSolution_ShuOsher(Conserved *U, double dt);

struct IntegrationMethod
{
  char Name[256];
  void (*AdvanceSolution)(Conserved *U, double dt);
} IntegrationMethods[4] =
  { { "FwdEuler", AdvanceSolution_FwdEuler },
    { "Midpoint", AdvanceSolution_Midpoint },
    { "ShuOsher", AdvanceSolution_ShuOsher }, { "", NULL }
  } ;



// Global variables go here
// =============================================================================
static double PlmTheta = 1.5;
static double AdiabaticGamma = 1.4;
static double CFL = 0.4;
static double OutputTime = 0.4;
static double MaxWavespeed = 0.0;
static int Nx = 400;
static int Ng = 2; // Number of ghost zones needed
static struct TestProblem *testProblem = &TestProblems[0];
static struct IntegrationMethod *integrationMethod = &IntegrationMethods[2];
// =============================================================================


#ifdef EULER1D_MAIN

int main(int argc, char **argv)
// =============================================================================
//                                 Main program
// =============================================================================
{
  RunUserMenu();

  FILE *OutputFile = fopen("euler.dat", "w");

  double *x = (double*) malloc((Nx+2*Ng)*sizeof(double));
  Conserved *U = (Conserved*) malloc((Nx+2*Ng)*sizeof(Conserved));
  int i;

  // ---------------------------------------------------------------------------
  // Initial conditions setup.
  for (i=0; i<Nx+2*Ng; ++i) {

    // Here we initialize the array of x coordinates, and set the initial
    // conditions on the conserved array.

    x[i] = -0.5 + (i-Ng+0.5) / Nx;

    Primitive P = testProblem->InitialCondition(x[i]);
    U[i] = PrimToCons(&P);
  }

  int CycleCounter = 0;
  double CurrentTime = 0.0;
  double EndTime = OutputTime;
  double dt = 0.0;
  double dx = 1.0 / Nx;

  // ---------------------------------------------------------------------------
  // This is the main iteration loop.
  while (CurrentTime < EndTime) {

    MaxWavespeed = 0.0;
    integrationMethod->AdvanceSolution(U, dt);

    CurrentTime += dt;
    CycleCounter++;
    dt = CFL * dx / MaxWavespeed;

    printf("\r[%s]", ProgressBar(CurrentTime/EndTime));
    fflush(stdout);
  }

  // ---------------------------------------------------------------------------
  // Output solution.
  for (i=Ng; i<Nx+Ng; ++i) {

    // This loop prints an ASCII table of x-coordinate values with their
    // associated primitive variables.

    Primitive P = ConsToPrim(&U[i]);
    fprintf(OutputFile, "%+8.6e %+8.6e %+8.6e %+8.6e %+8.6e %+8.6e\n",
            x[i], P.rho, P.p, P.vx, P.vy, P.vz);
  }

  free(x);
  free(U);
  fclose(OutputFile);

  // ---------------------------------------------------------------------------
  // Create a simple GNUplot script.
  FILE *gnuplot = fopen("plot.gp", "w");
  fprintf(gnuplot, "set title \"%s %s t=%3.2f PLM=%3.2f\"\n",
          testProblem->Name, integrationMethod->Name, OutputTime, PlmTheta);
  fprintf(gnuplot, "set xlabel \"x\"\n");
  fprintf(gnuplot, "plot \\\n");
  fprintf(gnuplot, "\"euler.dat\" u 1:2 title \"Density \",\\\n");
  fprintf(gnuplot, "\"euler.dat\" u 1:3 title \"Pressure\",\\\n");
  fprintf(gnuplot, "\"euler.dat\" u 1:4 title \"Velocity\"   \n");
  fprintf(gnuplot, "pause(-1)\\\n");
  fclose(gnuplot);

  printf("\n\n");
  printf("All done. Data is in euler.dat, but if you have gnuplot simply run\n"
	 "$> gnuplot plot.gp\n\n");

  return 0;
}

#endif // EULER1D_MAIN

void RunUserMenu()
{
  char input[200];
  int n;

  printf("Enter number of x-zones: [%d] ", Nx);
  fgets(input, 200, stdin);
  if (strlen(input) != 1) {
    Nx = atoi(input);
  }

  printf("Enter the output time: [%3.2f] ", OutputTime);
  fgets(input, 200, stdin);
  if (strlen(input) != 1) {
    OutputTime = atof(input);
  }

  printf("Enter PLM theta value: [%3.2f] ", PlmTheta);
  fgets(input, 200, stdin);
  if (strlen(input) != 1) {
    PlmTheta = atof(input);
  }

  printf("Enter Courant condition (CFL number): [%3.2f] ", CFL);
  fgets(input, 200, stdin);
  if (strlen(input) != 1) {
    CFL = atof(input);
  }

  // Choose the test problem to run.
  // ---------------------------------------------------------------------------
  struct TestProblem *tp = TestProblems;
  n = 0;
  printf("\n");
  while (tp->InitialCondition) {
    printf("%d: %s\n", n++, (tp++)->Name);
  }
  printf("Enter the problem setup to run: [%s] ", testProblem->Name);

  fgets(input, 200, stdin);
  if (strlen(input) != 1) {
    int choice = atoi(input);
    if (choice < n) {
      testProblem = &TestProblems[choice];
    }
    else {
      printf("That's not a choice.\n");
      exit(0);
    }
  }

  // Choose the integration method to use.
  // ---------------------------------------------------------------------------
  struct IntegrationMethod *ti = IntegrationMethods;
  n = 0;
  printf("\n");
  while (ti->AdvanceSolution) {
    printf("%d: %s\n", n++, (ti++)->Name);
  }
  printf("Enter the time integration to use: [%s] ", integrationMethod->Name);

  fgets(input, 200, stdin);
  if (strlen(input) != 1) {
    int choice = atoi(input);
    if (choice < n) {
      integrationMethod = &IntegrationMethods[choice];
    }
    else {
      printf("That's not a choice.\n");
      exit(0);
    }
  }
  printf("\n");
}


void SetBoundaryConditions(Conserved *U)
// -----------------------------------------------------------------------------
// This function sets outflow boundary conditions, also known as
// zero-gradient. The outermost Ng (number of ghost) zones are modified to take
// on the value of the nearest interior zone.
// -----------------------------------------------------------------------------
{
  int i;
  for (i=0; i<Ng; ++i) {
    U[i] = U[Ng];
  }
  for (i=Nx+Ng; i<Nx+2*Ng; ++i) {
    U[i] = U[Nx-1];
  }
}

Conserved RiemannSolver(const Primitive *Pl,
                        const Primitive *Pr, int dimension)
{
  Conserved Ul = PrimToCons(Pl);
  Conserved Ur = PrimToCons(Pr);

  Conserved Fl = FluxFunction(Pl, dimension);
  Conserved Fr = FluxFunction(Pr, dimension);

  double lamL[5]={0,0,0,0,0}, lamR[5]={0,0,0,0,0};
  Eigenvalues(Pl, dimension, lamL);
  Eigenvalues(Pr, dimension, lamR);

  const double am = minval3(lamL[0], lamR[0], 0.0); // left going wave speed
  const double ap = maxval3(lamL[4], lamR[4], 0.0); // right
  MaxWavespeed = maxval3(fabs(ap), fabs(am), MaxWavespeed);

  Conserved F;

  F.D   = (ap*Fl.D   - am*Fr.D   + ap*am*(Ur.D   - Ul.D  )) / (ap - am);
  F.tau = (ap*Fl.tau - am*Fr.tau + ap*am*(Ur.tau - Ul.tau)) / (ap - am);
  F.Sx  = (ap*Fl.Sx  - am*Fr.Sx  + ap*am*(Ur.Sx  - Ul.Sx )) / (ap - am);
  F.Sy  = (ap*Fl.Sy  - am*Fr.Sy  + ap*am*(Ur.Sy  - Ul.Sy )) / (ap - am);
  F.Sz  = (ap*Fl.Sz  - am*Fr.Sz  + ap*am*(Ur.Sz  - Ul.Sz )) / (ap - am);

  return F;
}

void ReconstructStates(const Primitive *P, Primitive *Pl, Primitive *Pr)
// -----------------------------------------------------------------------------
//The input value to this function is the primitive array, centered at zone i
// Its outputs are the extrapolated primitive states to either side of the i+1/2
// zone interface, which are used as input values for the Riemann solver by the
// calling function. The derivaties dP/dx are evaluated based on the generalized
// minmod slope limiter.
//
//               Pr_{i+1/2} = P_{i+1} - 0.5 * dx * (dP/dx)_{i+1}
//               Pl_{i+1/2} = P_{i}   + 0.5 * dx * (dP/dx)_{i}
//
// -----------------------------------------------------------------------------
{
  Pr->p = P[1].p - 0.5 * minmod3(P[ 0].p, P[1].p, P[2].p);
  Pl->p = P[0].p + 0.5 * minmod3(P[-1].p, P[0].p, P[1].p);

  Pr->rho = P[1].rho - 0.5 * minmod3(P[ 0].rho, P[1].rho, P[2].rho);
  Pl->rho = P[0].rho + 0.5 * minmod3(P[-1].rho, P[0].rho, P[1].rho);

  Pr->vx = P[1].vx - 0.5 * minmod3(P[ 0].vx, P[1].vx, P[2].vx);
  Pl->vx = P[0].vx + 0.5 * minmod3(P[-1].vx, P[0].vx, P[1].vx);

  Pr->vy = P[1].vy - 0.5 * minmod3(P[ 0].vy, P[1].vy, P[2].vy);
  Pl->vy = P[0].vy + 0.5 * minmod3(P[-1].vy, P[0].vy, P[1].vy);

  Pr->vz = P[1].vz - 0.5 * minmod3(P[ 0].vz, P[1].vz, P[2].vz);
  Pl->vz = P[0].vz + 0.5 * minmod3(P[-1].vz, P[0].vz, P[1].vz);
}

void TimeDerivative(const Conserved *U, Conserved *L)
// -----------------------------------------------------------------------------
// This function recieved the conserved quantities U, and returns their time
// derivative L.
// -----------------------------------------------------------------------------
{
  Primitive *P = (Primitive*) malloc((Nx+2*Ng)*sizeof(Primitive));
  Conserved *F = (Conserved*) malloc((Nx+2*Ng)*sizeof(Conserved));

  int i;
  const double dx = 1.0 / Nx;

  // (1) Convert the conserved array to primitives.
  // ---------------------------------------------------------------------------
  for (i=0; i<Nx+2*Ng; ++i) {
    P[i] = ConsToPrim(&U[i]);
  }

  // (2) Obtain the intercell fluxes of conserved variables with the HLL solver.
  // ---------------------------------------------------------------------------
  for (i=Ng-1; i<Nx+Ng; ++i) {

    Primitive Pl, Pr;
    ReconstructStates(&P[i], &Pl, &Pr);

    F[i] = RiemannSolver(&Pl, &Pr, 0); // F[i] := F^{HLL}_{i+1/2}
  }

  // (3) Difference the fluxes to obtain the time derivative, L.
  // ---------------------------------------------------------------------------
  for (i=1; i<Nx+2*Ng; ++i) {
    L[i].D   = -(F[i].D   - F[i-1].D  ) / dx;
    L[i].tau = -(F[i].tau - F[i-1].tau) / dx;
    L[i].Sx  = -(F[i].Sx  - F[i-1].Sx ) / dx;
    L[i].Sy  = -(F[i].Sy  - F[i-1].Sy ) / dx;
    L[i].Sz  = -(F[i].Sz  - F[i-1].Sz ) / dx;
  }

  free(P);
  free(F);
}

void AdvanceSolution_FwdEuler(Conserved *U, double dt)
{
  int i;

  Conserved *L  = (Conserved*) malloc((Nx+2*Ng)*sizeof(Conserved));

  SetBoundaryConditions(U);
  TimeDerivative(U, L);

  for (i=0; i<Nx+2*Ng; ++i) {
    U[i].D   += L[i].D   * dt;
    U[i].tau += L[i].tau * dt;
    U[i].Sx  += L[i].Sx  * dt;
    U[i].Sy  += L[i].Sy  * dt;
    U[i].Sz  += L[i].Sz  * dt;
  }

  free(L);
}
void AdvanceSolution_Midpoint(Conserved *U, double dt)
{
  int i;

  Conserved *L  = (Conserved*) malloc((Nx+2*Ng)*sizeof(Conserved));
  Conserved *U1 = (Conserved*) malloc((Nx+2*Ng)*sizeof(Conserved));

  SetBoundaryConditions(U);
  TimeDerivative(U, L);

  for (i=0; i<Nx+2*Ng; ++i) {
    U1[i].D   = U[i].D   + L[i].D   * 0.5 * dt;
    U1[i].tau = U[i].tau + L[i].tau * 0.5 * dt;
    U1[i].Sx  = U[i].Sx  + L[i].Sx  * 0.5 * dt;
    U1[i].Sy  = U[i].Sy  + L[i].Sy  * 0.5 * dt;
    U1[i].Sz  = U[i].Sz  + L[i].Sz  * 0.5 * dt;
  }

  SetBoundaryConditions(U1);
  TimeDerivative(U1, L);

  for (i=0; i<Nx+2*Ng; ++i) {
    U[i].D   = U[i].D   + L[i].D   * dt;
    U[i].tau = U[i].tau + L[i].tau * dt;
    U[i].Sx  = U[i].Sx  + L[i].Sx  * dt;
    U[i].Sy  = U[i].Sy  + L[i].Sy  * dt;
    U[i].Sz  = U[i].Sz  + L[i].Sz  * dt;
  }

  free(L);
  free(U1);
}

void AdvanceSolution_ShuOsher(Conserved *U, double dt)
{
  int i;

  Conserved *L  = (Conserved*) malloc((Nx+2*Ng)*sizeof(Conserved));
  Conserved *U1 = (Conserved*) malloc((Nx+2*Ng)*sizeof(Conserved));

  SetBoundaryConditions(U);
  TimeDerivative(U, L);

  for (i=0; i<Nx+2*Ng; ++i) {
    U1[i].D   = U[i].D   + L[i].D   * dt;
    U1[i].tau = U[i].tau + L[i].tau * dt;
    U1[i].Sx  = U[i].Sx  + L[i].Sx  * dt;
    U1[i].Sy  = U[i].Sy  + L[i].Sy  * dt;
    U1[i].Sz  = U[i].Sz  + L[i].Sz  * dt;
  }

  SetBoundaryConditions(U1);
  TimeDerivative(U1, L);

  for (i=0; i<Nx+2*Ng; ++i) {
    U1[i].D   = 3./4. * U[i].D   + 1./4. * U1[i].D   + 1./4. * L[i].D   * dt;
    U1[i].tau = 3./4. * U[i].tau + 1./4. * U1[i].tau + 1./4. * L[i].tau * dt;
    U1[i].Sx  = 3./4. * U[i].Sx  + 1./4. * U1[i].Sx  + 1./4. * L[i].Sx  * dt;
    U1[i].Sy  = 3./4. * U[i].Sy  + 1./4. * U1[i].Sy  + 1./4. * L[i].Sy  * dt;
    U1[i].Sz  = 3./4. * U[i].Sz  + 1./4. * U1[i].Sz  + 1./4. * L[i].Sz  * dt;
  }

  SetBoundaryConditions(U1);
  TimeDerivative(U1, L);

  for (i=0; i<Nx+2*Ng; ++i) {
    U[i].D   = 1./3. * U[i].D   + 2./3. * U1[i].D   + 2./3. * L[i].D   * dt;
    U[i].tau = 1./3. * U[i].tau + 2./3. * U1[i].tau + 2./3. * L[i].tau * dt;
    U[i].Sx  = 1./3. * U[i].Sx  + 2./3. * U1[i].Sx  + 2./3. * L[i].Sx  * dt;
    U[i].Sy  = 1./3. * U[i].Sy  + 2./3. * U1[i].Sy  + 2./3. * L[i].Sy  * dt;
    U[i].Sz  = 1./3. * U[i].Sz  + 2./3. * U1[i].Sz  + 2./3. * L[i].Sz  * dt;
  }

  free(L);
  free(U1);
}



Conserved PrimToCons(const Primitive *P)
// -----------------------------------------------------------------------------
// This function transforms the primitive variables into the conserved
// ones. There is no difficulty here.
// -----------------------------------------------------------------------------
{
  const double gm1 =  AdiabaticGamma - 1.0;
  const double v2  =  P->vx*P->vx + P->vy*P->vy + P->vz*P->vz;
  const double e   =  P->p / (P->rho * gm1); // specific internal energy

  Conserved U;

  U.D   = P->rho;
  U.tau = P->rho * (0.5*v2 + e);
  U.Sx  = P->rho * P->vx;
  U.Sy  = P->rho * P->vy;
  U.Sz  = P->rho * P->vz;

  return U;
}


Primitive ConsToPrim(const Conserved *U)
// -----------------------------------------------------------------------------
// This function transforms the conserved variables into the primitive
// ones. The inversion is trivial in the NR case.
// -----------------------------------------------------------------------------
{
  const double gm1 =  AdiabaticGamma - 1.0;
  const double S2  =  U->Sx*U->Sx + U->Sy*U->Sy + U->Sz*U->Sz;

  Primitive P;

  P.p   = (U->tau - 0.5*S2/U->D)*gm1;
  P.rho =  U->D;
  P.vx  =  U->Sx / U->D;
  P.vy  =  U->Sy / U->D;
  P.vz  =  U->Sz / U->D;

  if (P.p < 0.0) {
    printf("ConsToPrim got a negative pressure. Exiting.\n");
    exit(1);
  }
  if (P.rho < 0.0) {
    printf("ConsToPrim got a negative density. Exiting.\n");
    exit(1);
  }

  return P;
}



Conserved FluxFunction(const Primitive *P, int dimension)
// -----------------------------------------------------------------------------
// This function returns the flux of the conserved variables in the cartesian
// direction specified by the input parameter 'dimension' = 0,1,2
// -----------------------------------------------------------------------------
{
  Conserved U = PrimToCons(P);
  Conserved F = { 0,0,0,0,0 };

  switch (dimension) {
  case 0:
    F.D   = U.D  * P->vx;
    F.Sx  = U.Sx * P->vx + P->p;
    F.Sy  = U.Sy * P->vx;
    F.Sz  = U.Sz * P->vx;
    F.tau = (U.tau + P->p) * P->vx;
    break;
  case 1:
    F.D   = U.D  * P->vy ;
    F.Sx  = U.Sx * P->vy ;
    F.Sy  = U.Sy * P->vy + P->p;
    F.Sz  = U.Sz * P->vy ;
    F.tau = (U.tau + P->p) * P->vy;
    break;
  case 2:
    F.D   = U.D  * P->vz;
    F.Sx  = U.Sx * P->vz;
    F.Sy  = U.Sy * P->vz;
    F.Sz  = U.Sz * P->vz + P->p;
    F.tau = (U.tau + P->p) * P->vz;
    break;
  }
  return F;
}


void Eigenvalues(const Primitive *P, int dimension, double *evals)
// -----------------------------------------------------------------------------
// This function obtains the eigenvalues of the flux jacobian in the cartesian
// direction specified by the input parameter 'dimension' = 0,1,2. The result is
// placed in the variable 'evals' which must be an array of length 5.
// -----------------------------------------------------------------------------
{
  const double cs   =  sqrt(AdiabaticGamma * P->p / P->rho); // sound speed

  switch (dimension) {
  case 0:
    evals[0] =  P->vx - cs;
    evals[1] =  P->vx;
    evals[2] =  P->vx;
    evals[3] =  P->vx;
    evals[4] =  P->vx + cs;
    break;
  case 1:
    evals[0] =  P->vy - cs;
    evals[1] =  P->vy;
    evals[2] =  P->vy;
    evals[3] =  P->vy;
    evals[4] =  P->vy + cs;
    break;
  case 2:
    evals[0] =  P->vz - cs;
    evals[1] =  P->vz;
    evals[2] =  P->vz;
    evals[3] =  P->vz;
    evals[4] =  P->vz + cs;
    break;
  }
}




double maxval3(double a, double b, double c)
{
  return (a>b) ? ((a>c) ? a : c) : ((b>c) ? b : c);
}
double minval3(double a, double b, double c)
{
  return (a<b) ? ((a<c) ? a : c) : ((b<c) ? b : c);
}
double sign(double x)
{
  return (x > 0.0) - (x < 0.0);
}
double minmod3(double ul, double u0, double ur)
{
  const double a = PlmTheta * (u0 - ul);
  const double b =     0.5  * (ur - ul);
  const double c = PlmTheta * (ur - u0);
  const double sa = sign(a), sb = sign(b), sc = sign(c);
  const double fa = fabs(a), fb = fabs(b), fc = fabs(c);
  return 0.25*fabs(sa + sb)*(sa + sc)*minval3(fa, fb, fc);
}
char *ProgressBar(double frac)
{
  static const char ch[] =
    "=========================================================================";
  static char ret[100];
  const int count = frac * 72;
  sprintf(ret, "%2.1f%% %.*s", frac*100, count, ch);
  return ret;
}








// =============================================================================
// This is where new problem setups may be added.
// =============================================================================
Primitive InitialCondition_Shocktube1(double x)
{
  Primitive P;
  if (x < 0.0) {
    P.p   = 1.0;
    P.rho = 1.0;
    P.vx  = 0.0;
    P.vy  = 0.0;
    P.vz  = 0.0;
  }
  else {
    P.p   = 0.1;
    P.rho = 0.125;
    P.vx  = 0.0;
    P.vy  = 0.0;
    P.vz  = 0.0;
  }
  return P;
}

Primitive InitialCondition_Shocktube2(double x)
{
  Primitive P;
  if (x < 0.0) {
    P.p   = 0.95;
    P.rho = 1.08;
    P.vx  = 0.4;
    P.vy  = 0.3;
    P.vz  = 0.2;
  }
  else {
    P.p   = 1.0;
    P.rho = 1.0;
    P.vx  =-0.45;
    P.vy  =-0.20;
    P.vz  = 0.20;
  }
  return P;
}

Primitive InitialCondition_Blastwave1(double x) // Marti & Muller section 6.2
{
  Primitive P;
  if (x < 0.0) {
    P.p   = 13.33;
    P.rho = 10.0;
    P.vx  = 0.0;
    P.vy  = 0.0;
    P.vz  = 0.0;
  }
  else {
    P.p   = 0.01;
    P.rho = 1.0;
    P.vx  = 0.0;
    P.vy  = 0.0;
    P.vz  = 0.0;
  }
  return P;
}

Primitive InitialCondition_Blastwave2(double x) // Marti & Muller section 6.2
{
  Primitive P;
  if (x < 0.0) {
    P.p   = 1000.0;
    P.rho = 1.0;
    P.vx  = 0.0;
    P.vy  = 0.0;
    P.vz  = 0.0;
  }
  else {
    P.p   = 0.01;
    P.rho = 1.0;
    P.vx  = 0.0;
    P.vy  = 0.0;
    P.vz  = 0.0;
  }
  return P;
}
