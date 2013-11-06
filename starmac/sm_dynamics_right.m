function [ f, dfdx, dfdu ] = sm_dynamics_right( x, u, t, user )

Ninputs = user.Ninputs;

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

f = zeros( 6, 1 );
f( px ) = x(dotx);
f( pz ) = x(dotz);
f( theta ) = x(dottheta);
f( dotx ) = 1 / M * ( - Fvel * x(dotx) + sin( x(theta) ) * M * g );
f( dotz ) = 1 / M * ( - Fvel * x(dotz) + cos( x(theta) ) * M * g ) - g;
% f( dotx ) = 1 / M * ( - Fvel * x(dotx) + x(theta) * M * g );
% f( dotz ) = 1 / M * ( - Fvel * x(dotz) + x(theta) * M * g ) - g;
f( dottheta ) = 1 / I * ( - Fangle * x(dottheta) - L * u( 1 ) );


if( nargout > 1 )
  dfdx = zeros( 6, 6 );
  dfdx( px, dotx ) = 1;
  dfdx( pz, dotz ) = 1;
  dfdx( theta, dottheta ) = 1;
  dfdx( dotx, dotx ) = - Fvel / M;
  dfdx( dotx, theta ) = 1 / M * cos( x(theta) ) * M * g;
%   dfdx( dotx, theta ) = g;
  dfdx( dotz, dotz ) = - Fvel / M;
  dfdx( dotz, theta ) = - 1 / M * sin( x(theta) ) * M * g;
%   dfdx( dotz, theta ) = g;
  dfdx( dottheta, dottheta ) = - Fangle / I;

  dfdu = zeros( 6, Ninputs );
  dfdu( dottheta )  = - L / I;

end
