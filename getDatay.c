/* Title: Getting y data from Basilisk file
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"
char filename[80];
face vector shear[], D2p[];

event init(t = 0)
{
  restore (file = filename);
  boundary(all);
  FILE * fp = ferr;
  int LEVEL = 8;

  foreach_face(y){
    shear.y[] = (u.x[0,0] - u.x[0,-1])/(Delta);
  }
  boundary ((scalar *){shear});

  foreach_face(y) {
    double D11 = (u.y[] - u.y[0,-1]);
    double D22 = ((u.x[1,0]-u.x[-1,0])+(u.x[1,-1]-u.x[-1,-1]))/4.0;
    double D12 = 0.5*(((u.y[1,0]-u.y[-1,0])+(u.y[1,-1]-u.y[-1,-1]))/4.0 + (u.x[] - u.x[0,-1]));
    double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12))/(Delta);
    // D2 /= sqrt(2.0); // note that this is $\|D_{ij}\|$
    D2p.y[] = D2;
  }
  boundary ((scalar *){D2p});

  for (double y = 0.; y < 1.0; y += 1./pow(2.,LEVEL)){
    fprintf(ferr, "%g %g %g %g\n", y, interpolate(u.x, L0/2.0, y),
            interpolate(shear.y, L0/2.0, y), interpolate(D2p.y, L0/2.0, y));
  }
  fflush (fp);
  fclose (fp);
}

int main(int a, char const *arguments[])
{
  // Boundary condition: periodic right - left
  periodic (right);
  u.t[top] = neumann(0);
  u.n[top] = neumann(0);
  uf.n[top] = neumann(0);
  // bottom is moving wall
  u.n[bottom] = dirichlet(0);
  uf.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  p[top] = neumann(0);
  pf[top] = neumann(0);
  p[bottom] = neumann(0);
  pf[bottom] = neumann(0);
  sprintf (filename, "%s", arguments[1]);
  run();
}
