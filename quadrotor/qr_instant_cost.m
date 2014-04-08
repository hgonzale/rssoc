function [ L, dLdx, dLdu ] = qr_instant_cost( x, u, t, user )
%This function calculates the cost function, L, and its derivatives with
%respect to the states x and the input u 

Nstates = user.Nstates;
Ninputs = user.Ninputs;

Qvx = user.cost.instant_cost.Qvx; 
Qvy = user.cost.instant_cost.Qvy; 
Qvz = user.cost.instant_cost.Qvz;
Qqtheta = user.cost.instant_cost.Qqtheta; 
Qqphi = user.cost.instant_cost.Qqphi;
Qvtheta = user.cost.instant_cost.Qvtheta;
Qvphi = user.cost.instant_cost.Qvphi;
Qvpsi = user.cost.instant_cost.Qvpsi;
R = user.cost.instant_cost.R; 

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

Q = zeros(Nstates,Nstates); 
Q(vx,vx) = Qvx; 
Q(vy,vy) = Qvy; 
Q(vz,vz) = Qvz; 
Q(qtheta,qtheta) = Qqtheta; 
Q(qphi,qphi) = Qqphi; 
Q(vtheta,vtheta) = Qvtheta; 
Q(vphi,vphi) = Qvphi; 
Q(vpsi, vpsi) = Qvpsi; 

L = x' * Q * x + u' * R * u; 
 
if( nargout >= 2 )
    dLdx = 2 * x' * Q;  
    dLdu = 2 * R * u;
end

end

