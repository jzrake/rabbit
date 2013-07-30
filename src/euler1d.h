#ifndef EULER1D_H
#define EULER1D_H

typedef struct PrimitiveVariables
// -----------------------------------------------------------------------------
// p    : Gas pressure
// rho  : Mass density (fluid rest frame)
// v    : Fluid 3-velocity
// -----------------------------------------------------------------------------
{
  double p, rho, vx, vy, vz;
} Primitive;


typedef struct ConservedVariables
// -----------------------------------------------------------------------------
// D    : Mass density (lab frame)
// tau  : Total energy
// S    : Momentum vector
// -----------------------------------------------------------------------------
{
  double D, tau, Sx, Sy, Sz;
} Conserved;

// Prototypes of functions which are contained in this file
// -----------------------------------------------------------------------------
void ReconstructStates(const Primitive *P, Primitive *Pl, Primitive *Pr);
void Eigenvalues(const Primitive *P, int dimension, double *evals);
Primitive ConsToPrim(const Conserved *U);
Conserved PrimToCons(const Primitive *P);
Conserved FluxFunction(const Primitive *P, int dimension);
Conserved RiemannSolver(const Primitive *Pl, const Primitive *Pr, int dimension);

#endif // EULER1D_H
