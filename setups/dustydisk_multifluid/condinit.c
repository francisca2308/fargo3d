#include "fargo3d.h"

void _CondInit(int id) {
  
  int i,j,k;
    int n;

  real *sigma  = Density->field_cpu;
  real *vphi   = Vx->field_cpu;
  real *vr     = Vy->field_cpu;
  real *cs     = Energy->field_cpu;

  real dlogPdlogr = -1.5;
  real eta = 0.5*ASPECTRATIO*ASPECTRATIO*dlogPdlogr;
  real beta = pow(1+2*eta,0.5);
  real xi   = 0.5*beta;
  real psi;
  real Sn, Qn;
  
  real stokes_plus[NFLUIDS];
  real stokes[NFLUIDS-1];
  real epsilons[NFLUIDS-1];
  real vrgas;
  real vphigas;
  real vrdust;
  real vphidust;
  real sq = SQ;
  real slope;
  // Stokes numbers
  real stmax = TSMAX;
  real stmin = TSMIN;
  real dst   = (log(stmax)-log(stmin))/(NFLUIDS-1);
  for(n=0;n<NFLUIDS;n++){
    stokes_plus[n] = stmin*exp(dst*n);
  }

  //Density distribution
  slope = 4-sq;
  for(n=0; n<NFLUIDS-1; n++){
    if(slope != 0.) {
      epsilons[n]  = pow(stokes_plus[n+1],slope) - pow(stokes_plus[n],slope) ;
      epsilons[n] *= EPSILON/(pow(stmax, slope) - pow(stmin,slope));
    }
    else{
      epsilons[n]  = log(stokes_plus[n+1]/stokes_plus[n]);
      epsilons[n] *= EPSILON/log(stmax/stmin);
    }
    stokes[n] = stokes_plus[n+1];
    if( NFLUIDS == 2) stokes[n] = TSMAX;
    if(CPU_Master) printf("%1.16f \t %1.16f \n", stokes[n], epsilons[n]);
  }
#ifdef DRAGFORCE
  if(id > 0) {
#ifdef STOKESNUMBER
   Coeffval[0]   = 1.0/stokes[id-1];
#endif
#ifdef DUSTSIZE
    Coeffval[1]   = 1.0/(stokes[id-1]*R0/R0_CGS);    
    Coeffval[2]   = RHOSOLID/(MSTAR_CGS/(R0_CGS*R0_CGS*R0_CGS))*(MSTAR/(R0*R0*R0));
#endif
    if(CPU_Master) printf("Ts %1.16f \t eps %1.16f \n", stokes[id-1], epsilons[id-1]);
  }
#endif



  //summation
  Sn = 0.0;
  Qn = 0.0;

  for(n=0; n<NFLUIDS-1; n++){
    Qn += epsilons[n]*stokes[n]/(1. + stokes[n]*stokes[n]);
    Sn += epsilons[n]/(1. + stokes[n]*stokes[n]);
  }
  psi = 1./( (Sn+beta)*(Sn+2.*xi) + Qn*Qn);


  i = j = k = 0;
  
  for (j=0; j<Ny+2*NGHY; j++) {
    for (i=0; i<Nx+2*NGHX; i++) {
      
      real vkmed = sqrt(G*MSTAR/ymed(j));
      real vkmin = sqrt(G*MSTAR/ymin(j));

     
      vphigas = -((Sn+2*xi)*Sn+Qn*Qn)*psi*(beta-1)*vkmed;
      vrgas   = -2*beta*psi*(beta-1)*Qn*vkmin;

      vrdust  = 2*stokes[id-1]*(beta-1)*vkmed;
      vrdust += vrgas + 2*stokes[id-1]*vphigas;
      vrdust /= (1+stokes[id-1]*stokes[id-1]);

      vphidust  = 2.*(beta-1)*vkmed;
      vphidust += 2*vphigas - stokes[id-1]*vrgas;
      vphidust /= 2*(1.+stokes[id-1]*stokes[id-1]);

      if (id == 0) {
        sigma[l] = SIGMA0*pow(ymed(j)/R0,-0.5);
        vr[l]    = vrgas;
        vphi[l]  = beta*vkmed;
        vphi[l] += vphigas;
        cs[l]    = vkmed*ASPECTRATIO;
      }
      if (id > 0) {
        sigma[l] = epsilons[id-1]*SIGMA0*pow(ymed(j)/R0,-0.5);
        vr[l]    = vrdust;
        vphi[l]  = vkmed;
        vphi[l] += vphidust;
        cs[l]    = 0.0; // Pressureless fluid
      }
      vphi[l] -= OMEGAFRAME*ymed(j);
    }
  } 
  #ifdef COLLISIONS
   if(id>0)  ColRate(1.0/stokes[id-1], 0, id, YES);
  #endif
}

void CondInit() {
  int id_gas = 0;
  int feedback = YES;
  //We first create the gaseous fluid and store it in the array Fluids[]
  Fluids[id_gas] = CreateFluid("gas",GAS);

  //We now select the fluid
  SelectFluid(id_gas);

  //and fill its fields
  _CondInit(id_gas);

  //We repeat the process for the dust fluids
  char dust_name[MAXNAMELENGTH];
  int id_dust;

  for(id_dust = 1; id_dust<NFLUIDS; id_dust++) {
    sprintf(dust_name,"dust%d",id_dust); //We assign different names to the dust fluids

    Fluids[id_dust]  = CreateFluid(dust_name, DUST);
    SelectFluid(id_dust);
    _CondInit(id_dust);

  }
}
