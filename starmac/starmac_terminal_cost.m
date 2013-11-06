function [ phi, dphidx ] = starmac_terminal_cost( x, waypoint, user )

Nstates = user.Nstates;
pos = user.idxs.pos;
theta = user.idxs.theta;
K = user.cost.terminal_cost.K;
Ktheta = user.cost.terminal_cost.Ktheta;

dist = x(pos) - waypoint(pos);
phi = 0.5 * ( K * dist' * dist ) + Ktheta * sin( 0.5 * ( x(theta) - waypoint(theta) ) )^2;

if( nargout >= 2 )
  dphidx = zeros( 1, Nstates );
  dphidx(pos) = K * dist;
  dphidx(theta) = Ktheta * sin( 0.5 * ( x(theta) - waypoint(theta) ) ) * cos( 0.5 * ( x(theta) - waypoint(theta) ) );
end
