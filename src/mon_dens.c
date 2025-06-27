//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void mon_dens_cpu() {

//<USER_DEFINED>
  INPUT(Density);
  OUTPUT(Slope);
  real xp = Xplanet; 
  real yp = Yplanet;
  real mass = MplanetVirtual;
//<\USER_DEFINED>


//<EXTERNAL>// real MplanetVirtual(1);
  real* dens = Density->field_cpu;
  real* interm = Slope->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real rh = pow(mass/3./MSTAR, 1./3.)*sqrt(xp*xp+yp*yp);
//<\EXTERNAL>

//<INTERNAL>/home/fangulo/myforkfargo3d/fargo3d/src/mon_dens.c
  int i;
  int j;
  int k;
  int ll;
  real dist2;
  real dx;
  real dy;
  real dz;
  real planet_distance;
//<\INTERNAL>

//<CONSTANT>
// real Sxi(Nx);
// real Syk(Nz+2*NGHZ);
// real InvVj(Ny+2*NGHY);
// real Xplanet(1);
// real Yplanet(1);
// real MplanetVirtual(1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
#ifdef CYLINDRICAL
	dx = ymed(j)*cos(xmed(i))-xp;
	dy = ymed(j)*sin(xmed(i))-yp;
#endif
	dist2 = dx*dx+dy*dy+dz*dz;
  planet_distance=sqrt(dist2);

	ll = l;
  if (planet_distance > rh)
	  interm[ll] = 0.0;
  else{
	  interm[ll] = dens[ll]*Vol(i,j,k);
  }
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}