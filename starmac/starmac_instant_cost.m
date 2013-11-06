function [ L, dLdx, dLdu ] = starmac_instant_cost( x, u, t, user )

Nstates = user.Nstates;
Ninputs = user.Ninputs;

Ktime = user.cost.instant_cost.Ktime;
Q = user.cost.instant_cost.Q;
Qdottheta = user.cost.instant_cost.Qdottheta;
zground = user.cost.instant_cost.zground;
roll = user.cost.instant_cost.roll;
t1 = user.cost.instant_cost.t1;
t2 = user.cost.instant_cost.t2;
t3 = user.cost.instant_cost.t3;
kx = user.cost.instant_cost.kx;
kz = user.cost.instant_cost.kz;
kt = user.cost.instant_cost.kt;

if( t < t1 )
  xd = t;
  zd = zground;
  thetad = roll;
elseif( t < t1 + t2 )
  xd = t1 + cos( ( t - t1 ) - pi/2 );
  zd = zground + 2 + 2 * sin( ( t - t1 ) - pi/2 );
  thetad = ( t - t1 );
elseif( t < t1 + t2 + t3 )
  xd = t - t2;
  zd = zground;
  thetad = roll;
else
  xd = t1 + t3;
  zd = zground;
  thetad = roll;
end

px = user.idxs.px;
pz = user.idxs.pz;
theta = user.idxs.theta;
dottheta = user.idxs.dottheta;

L = 0.5 * kx * ( x(px) - xd )^2 + 0.5 * kz * ( x(pz) - zd )^2 + kt * sin(0.5 * ( x(theta) - thetad ) )^2 + Ktime + 0.5 * u' * Q * u + 0.5 * x(dottheta)' * Qdottheta * x(dottheta);

if( nargout >= 2 )
  dLdx = zeros( 1, Nstates );
  dLdx(px) = kx * ( x(px) - xd );
  dLdx(pz) = kz * ( x(pz) - zd );
  dLdx(theta) = kt * sin( 0.5 * ( x(theta) - thetad ) ) * cos( 0.5 * ( x(theta) - thetad ) );
  dLdx(dottheta) = x(dottheta)' * Qdottheta;
  dLdu = u' * Q;
end
