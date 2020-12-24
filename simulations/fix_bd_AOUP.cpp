/* ----------------------------------------------------------------------
   This fix requires the custom atom style "active" which includes
   a 3D per-atom "noise" variable and a 2D per atom "v_a" variable.
   These styles are implemented in atom.cpp and atom.h.

   For a 2D simulation just set x[i][2] += 0.
------------------------------------------------------------------------- */

#include <iostream>
#include "stdio.h"
#include "string.h"
#include "fix_bd_AOUP.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "math.h"
#include "math_const.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
/* ---------------------------------------------------------------------- */

// example command
// fix bd_active TEMP GAMMA AMP SEED


FixBD_AOUP::FixBD_AOUP(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"bd_AOUP") != 0 && narg < 7)
    error->all(FLERR,"Illegal fix bd_AOUP command");
  t_target = force->numeric(FLERR,arg[3]);
  gamma = force->numeric(FLERR,arg[4]);
  tau = force->numeric(FLERR,arg[5]);
  Pe = force->numeric(FLERR,arg[6]);
  frac = force->numeric(FLERR,arg[7]);
  int seed = force->inumeric(FLERR,arg[8]);

  // allocate the random number generator
  random = new RanPark(lmp,seed); //Park-Miller random number generator, should return RN between -1 and 1
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixBD_AOUP::setmask()
{
  int mask = 0;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBD_AOUP::init()
{
  dt = update->dt;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

// variables:
// frac =  fraction of active particles. Set frac = 1 for all active particles.
//  This code can be used to get a mixture of passive and active particles.
// t_target = target temperature
// fd_term = standard deviation of particle noise
// fda_term = standard deviation of active noise


void FixBD_AOUP::final_integrate()
{
  // friction coefficient, this taken to be a property of the solvent
  // so here gamma_i is gamma / m_i
  double fd_term = 0.;
  double fda_term = 0.;
  double **noise = atom->noise;
  double **x = atom->x; //** because it is a 2D array
  double **f = atom->f;
  double **v_a = atom->v_a;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  int npassive=(int)(floor(nlocal*frac));

  for (int i = 0; i < npassive; i++) {
    if (mask[i] & groupbit) {
      gamma_i = gamma / mass[type[i]];
      // in LJ units, t_target is given in kbT/epsilon. Thus D = t_target * epsilon/gamma_i, but we set epsilon to 1.
      fd_term = sqrt(2 * dt * t_target / (gamma_i));
      noise[i][0] = fd_term * random->gaussian();
      noise[i][1] = fd_term * random->gaussian();
      v_a[i][0] += 0;
      v_a[i][1] += 0;
      x[i][0] += (dt / gamma_i) * (f[i][0] + v_a[i][0]) + noise[i][0]; //x-coords of particle
      x[i][1] += (dt / gamma_i) * (f[i][1] + v_a[i][1]) + noise[i][1]; //y-coords
    }
  }
  for (int i = npassive; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      gamma_i = gamma / mass[type[i]];
      // in LJ units, t_target is given in kbT/epsilon. Thus D = t_target * epsilon/gamma_i, but we set epsilon to 1.
      fd_term = sqrt(2 * dt * t_target / (gamma_i));
      fda_term = sqrt(2 * dt * Pe * Pe / tau);    // active noise. Scales like f^2 and 1/tau
      noise[i][0] = fd_term * random->gaussian();
      noise[i][1] = fd_term * random->gaussian();
      x[i][0] += (dt / gamma_i) * (f[i][0] + v_a[i][0]) + noise[i][0]; //x-coords of particle
      x[i][1] += (dt / gamma_i) * (f[i][1] + v_a[i][1]) + noise[i][1]; //y-coords
      v_a[i][0] += -(dt / tau)*v_a[i][0] + fda_term * random->gaussian();
      v_a[i][1] += -(dt / tau)*v_a[i][1] + fda_term * random->gaussian();
     }
  }
}
