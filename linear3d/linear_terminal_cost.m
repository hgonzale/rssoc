function [ phi, dphidx ] = linear_terminal_cost( x, waypoint, user )

Qx = user.cost.terminal_cost.Qx;

dist = x - waypoint;
phi = 0.5 * dist' * Qx * dist;

if( nargout > 1 )
  dphidx = dist' * Qx;
end
