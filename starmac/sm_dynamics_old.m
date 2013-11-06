function [ f, dfdx, dfdu, dfdt ] = sm_dynamics( x, u, t, user )

L = user.starmac_params.length;
M = user.starmac_params.mass;
I = user.starmac_params.Iy;
Fvel = user.starmac_params.friction_velocity;
Fangle = user.starmac_params.friction_angle;
g = user.starmac_params.gravity_accel;

px = user.idxs.px;
pz = user.idxs.pz;
theta = user.idxs.theta;
dotx = user.idxs.dotx;
dotz = user.idxs.dotz;
dottheta = user.idxs.dottheta;
T1 = user.idxs.T1;
T2 = user.idxs.T2;

f = zeros( 6, 1 );
f( px ) = x(dotx);
f( pz ) = x(dotz);
f( theta ) = x(dottheta);
f( dotx ) = 1 / M * ( - Fvel * x(dotx) + sin( x(theta) ) * ( u(T1) + u(T2) ) );
f( dotz ) = 1 / M * ( - Fvel * x(dotz) + cos( x(theta) ) * ( u(T1) + u(T2) ) ) - g;
f( dottheta ) = 1 / I * ( - Fangle * x(dottheta) + L * ( - u(T1) + u(T2) ) );


if( nargout > 1 )
  dfdx = zeros( 6, 6 );
  dfdx( px, dotx ) = 1;
  dfdx( pz, dotz ) = 1;
  dfdx( theta, dottheta ) = 1;
  dfdx( dotx, dotx ) = - Fvel / M;
  dfdx( dotx, theta ) = 1 / M * cos( x(theta) ) * sum(u);
  dfdx( dotz, dotz ) = - Fvel / M;
  dfdx( dotz, theta ) = - 1 / M * sin( x(theta) ) * sum(u);
  dfdx( dottheta, dottheta ) = - Fangle / I;

  dfdu = zeros( 6, 2 );
  dfdu( dotx, T1 ) = 1 / M * sin( x(theta) );
  dfdu( dotx, T2 ) = 1 / M * sin( x(theta) );
  dfdu( dotz, T1 ) = 1 / M * cos( x(theta) );
  dfdu( dotz, T2 ) = 1 / M * cos( x(theta) );
  dfdu( dottheta, T1 ) = - L / I;
  dfdu( dottheta, T2 ) = L / I;

  dfdt = zeros( 6, 1 );
end
