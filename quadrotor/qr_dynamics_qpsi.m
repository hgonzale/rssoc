function [ f, dfdx, dfdu ] = qr_dynamics_qpsi( x, u, t, user )
%This function calculates the dynnamics equations and its derivatives with
%respect to the state and input for the mode of operation that corresponds
%to movement in the qpsi direction 

Nstates = user.Nstates;  
Ninputs = user.Ninputs;

Lr = user.quadrotor_params.length;
m = user.quadrotor_params.mass;
M = user.quadrotor_params.Iy;
g = user.quadrotor_params.gravity_accel;
Kf = user.quadrotor_params.motor_const; 
Ki = user.quadrotor_params.air_const; 
b = user.quadrotor_params.frict_coeff; 

qx = user.idxs.qx;
qy = user.idxs.qy;
qz = user.idxs.qz;
vx = user.idxs.vx;
vy = user.idxs.vy;
vz = user.idxs.vz;
qtheta = user.idxs.qtheta;
qphi = user.idxs.qphi;
qpsi = user.idxs.qpsi;
vtheta = user.idxs.vtheta;
vphi = user.idxs.vphi;
vpsi = user.idxs.vpsi;

f = zeros( Nstates, 1 );
f( qx ) = x(vx);
f( qy ) = x(vy);
f( qz ) = x(vz);
f( vx ) = -b / m * x(vx) + g * tan(x(qphi)) / cos(x(qtheta)); 
f( vy ) = -b / m * x(vy) - g * tan(x(qtheta)); 
f( vz ) = -b / m * x(vz);
f( qtheta ) = x(vtheta); 
f( qphi ) = x(vphi); 
f( qpsi ) = x(vpsi); 
f( vtheta ) = -x(vphi) * x(vpsi); 
f( vphi ) = x(vtheta) * x(vpsi); 
f( vpsi ) = u; 

if( nargout > 1 )
    dfdx = zeros( Nstates, Nstates );
    dfdx( qx, vx ) = 1;
    dfdx( qy, vy ) = 1;
    dfdx( qz, vz ) = 1;
    dfdx( vx, vx ) = -b/m; 
    dfdx( vy, vy ) = -b/m; 
    dfdx( vz, vz ) = -b/m; 
    dfdx( vx, qtheta ) = (g * tan(x(qtheta)) * tan(x(qphi))) / cos(x(qtheta));  
    dfdx( vx, qphi ) = g / ( cos(x(qtheta)) * (cos(x(qphi)))^2); 
    dfdx( vy, qtheta ) = -g / (cos(x(qtheta)))^2; 
    dfdx( qtheta, vtheta ) = 1; 
    dfdx( qphi, vphi ) = 1; 
    dfdx( qpsi, vpsi ) = 1; 
    dfdx( vtheta, vphi ) = -x(vpsi); 
    dfdx( vtheta, vpsi ) = -x(vphi); 
    dfdx( vphi, vtheta ) = x(vpsi); 
    dfdx( vphi, vpsi) = x(vtheta); 
    
    dfdu = zeros( Nstates, Ninputs );
    dfdu( vpsi ) = 1; 

end