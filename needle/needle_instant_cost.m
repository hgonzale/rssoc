function [ L, dLdx, dLdu, dLdt ] = needle_instant_cost( x, u, t, user )

Nstates = user.Nstates;
Ninputs = user.Ninputs;

Ktime = user.cost.instant_cost.Ktime;
Q = user.cost.instant_cost.Q;

L = 0.5 * u' * Q * u + Ktime;

if( nargout >= 2 )
  dLdx = zeros( 1, Nstates );
  dLdu = u' * Q;
  dLdt = 0;
end
