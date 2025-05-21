#include "fargo3d.h"

void _Init(int index) {
  
  int i,j,k;
  
  real *sigma  = Density->field_cpu;
  real *vphi   = Vx->field_cpu;
  real *vr     = Vy->field_cpu;
  real *cs     = Energy->field_cpu;

  real eta = 0.75*ASPECTRATIO*ASPECTRATIO;
  
  i = j = k = 0;
  
  for (j=0; j<Ny+2*NGHY; j++) {
    for (i=0; i<Nx+2*NGHX; i++) {
      
      real vkmed = sqrt(G*MSTAR/ymed(j));
      real vkmin = sqrt(G*MSTAR/ymin(j));

      if (index == 0) {
        sigma[l] = SIGMA0*pow(ymed(j)/R0,-0.5);
        vr[l]    = 2*EPSILON1*TS1/(pow(TS1,2) + pow(1+EPSILON1,2))*eta*vkmin;
        vphi[l]  = vkmed;
        vphi[l] -= (pow(TS1,2) + 1 + EPSILON1)/(pow(TS1,2)+pow(1+EPSILON1,2))*eta*vkmed;
        cs[l]    = vkmed*ASPECTRATIO;
      }
      if (index == 1) {
        sigma[l] = EPSILON1*SIGMA0*pow(ymed(j)/R0,-0.5);
        vr[l]    = -2*TS1/(pow(TS1,2) + pow(1+EPSILON1,2))*eta*vkmin;
        vphi[l]  = vkmed;
        vphi[l] -= (1 + EPSILON1)/(pow(TS1,2)+pow(1+EPSILON1,2))*eta*vkmed;
        cs[l]    = 0.0; // Pressureless fluid
      }
      if (index == 2) {
        sigma[l] = EPSILON2*SIGMA0*pow(ymed(j)/R0,-0.5);
        vr[l]    = -2*TS2/(pow(TS2,2) + pow(1+EPSILON2,2))*eta*vkmin;
        vphi[l]  = vkmed;
        vphi[l] -= (1 + EPSILON2)/(pow(TS2,2)+pow(1+EPSILON2,2))*eta*vkmed;
        cs[l]    = 0.0; // Pressureless fluid
      }
    }
  } 
}

void CondInit() {
   Fluids[0] = CreateFluid("gas",GAS);
   Fluids[1] = CreateFluid("dust1",DUST);
   Fluids[2] = CreateFluid("dust2",DUST);
   
   SelectFluid(0);
   _Init(0);
   
   SelectFluid(1);
   _Init(1);

   SelectFluid(2);
   _Init(2);
     
   ColRate(1.0/TS1, 0, 1, YES); //Fill the collision matrix with the
   ColRate(1.0/TS2, 0, 2, YES); 
                               //inverse of the Stokes Number
}