function [ L, dLdx, dLdu ] = linear_quad_cost( x, u, t, user )

Ktime = user.cost.instant_cost.Ktime;
Qx = user.cost.instant_cost.Qx;
Qu = user.cost.instant_cost.Qu;

L = 0.5 * ( x' * Qx * x + u' * Qu * u ) + Ktime;

if( nargout >= 2 )
  dLdx = x' * Qx;
  dLdu = u' * Qu;
end
