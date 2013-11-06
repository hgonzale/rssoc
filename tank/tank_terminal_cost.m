function [ phi, dphidx, dphidt ] = tank_terminal_cost( x, waypoint, user )

Nstates = user.Nstates;
pos = user.idxs.pos;
K = user.cost.terminal_cost.K;

dist = x(pos) - waypoint(pos);
phi = 0.5 * ( K * dist' * dist );

if( nargout >= 2 )
  dphidx = zeros( 1, Nstates );
  dphidx(pos) = K * dist;
  dphidt = 0;
end
