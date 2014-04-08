function [ fx, dfdx ] = qr_instant_cons( x, user )
%This function calculates the constraints and its corresponding derivatives

Nstates = user.Nstates; 
len_cons = user.len_cons; 
x_min = user.x_min; 
x_max = user.x_max; 

qx = user.idxs.qx;
qy = user.idxs.qy;
qz = user.idxs.qz;
qtheta = user.idxs.qtheta;
qphi = user.idxs.qphi;

assert( len_cons == 10 );

fx = zeros( len_cons, 1 );
fx(1) = -x(qx) + x_min(qx); %lower constraint on qx 
fx(2) = x(qx) - x_max(qx); %upper constraint on qx 
fx(3) = -x(qy) + x_min(qy); 
fx(4) = x(qy) - x_max(qy); 
fx(5) = -x(qz) + x_min(qz); 
fx(6) = x(qz) - x_max(qz);
fx(7) = -x(qtheta) + x_min(qtheta); 
fx(8) = x(qtheta) - x_max(qtheta);
fx(9) = -x(qphi) + x_min(qphi); 
fx(10) = x(qphi) - x_max(qphi); 

if( nargout >= 2 )
  dfdx = zeros( len_cons, Nstates );
  dfdx(1,qx) = -1;
  dfdx(2,qx) = 1;
  dfdx(3,qy) = -1;
  dfdx(4,qy) = 1;
  dfdx(5,qz) = -1;
  dfdx(6,qz) = 1;
  dfdx(7,qtheta) = -1;
  dfdx(8,qtheta) = 1;
  dfdx(9, qphi) = -1;
  dfdx(10, qphi) = 1;
end


end

