function [ phi, dphidx ] = qr_terminal_cost( x, user )
%This function calculates the terminal cost, Phi, and its derivative with
%respect to the states, x
%note: the derivative with respect to the input u is not necessary because
%u does not appear in the equation for Phi

Nstates = user.Nstates;

qxd = user.cost.terminal_cost.qxd;
qyd = user.cost.terminal_cost.qyd; 
qzd = user.cost.terminal_cost.qzd;
qpsid = user.cost.terminal_cost.qpsid;
K = user.cost.terminal_cost.K;

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



phi = ( x(qx) - qxd )^2 + ( x(qy) - qyd )^2 + ( x(qz) - qzd )^2 + ( x(qpsi) - qpsid )^2 + ...
  x(vx)^2 + x(vy)^2 + x(vz)^2 + x(qtheta)^2 + x(qphi)^2 + x(vtheta)^2 + x(vphi)^2 + x(vpsi)^2;
phi = K * phi;

if( nargout >= 2 )
  dphidx = zeros( 1, Nstates );
  dphidx(qx) = 2 * ( x(qx) - qxd ); 
  dphidx(qy) = 2 * ( x(qy) - qyd );
  dphidx(qz) = 2 * ( x(qz) - qzd );
  dphidx(qpsi) = 2 * ( x(qpsi) - qpsid );
  dphidx(vx) = 2 * x(vx); 
  dphidx(vy) = 2 * x(vy); 
  dphidx(vz) = 2 * x(vz);
  dphidx(qtheta) = 2 * x(qtheta);
  dphidx(qphi) = 2 * x(qphi);
  dphidx(vtheta) = 2 * x(vtheta);
  dphidx(vphi) = 2 * x(vphi); 
  dphidx(vpsi) = 2 * x(vpsi);
  dphidx = K * dphidx;
end

end

