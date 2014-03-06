function [ f, dfdx, dfdu ] = sm_dynamics( x, u, t, user )

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
T3 = user.idxs.T3;

f = [ x(dotx)                              ; ...
      x(dotz)                              ; ...
      x(dottheta)                          ; ...
      1/M * sin( x(theta) ) * sum(u)       ; ...
      1/M * cos( x(theta) ) * sum(u) - g ; ...
      L/I * ( u(T1) - u(T3) )              ];

if( nargout >= 3 )
  dfdx = [  0  0   0                              1  0  0 ; ...
            0  0   0                              0  1  0 ; ...
            0  0   0                              0  0  1 ; ...
            0  0   1/M * cos( x(theta) ) * sum(u) 0  0  0 ; ...
            0  0  -1/M * sin( x(theta) ) * sum(u) 0  0  0 ; ...
            0  0   0                              0  0  0 ];

  dfdu = [ 0                 0                  0                 ; ...
           0                 0                  0                 ; ...
           0                 0                  0                 ; ...
           1/M*sin(x(theta)) 1/M*sin(x(theta))  1/M*sin(x(theta)) ; ...
           1/M*cos(x(theta)) 1/M*cos(x(theta))  1/M*cos(x(theta)) ; ...
           L/I               0                 -L/I               ];
       
end