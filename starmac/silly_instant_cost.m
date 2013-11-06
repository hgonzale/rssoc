function [ L, dLdx, dLdu ] = silly_instant_cost( x, u, t, user )

Nstates = user.Nstates;
Q = user.cost.instant_cost.Q;

L = 0.5 * u' * Q * u;

if( nargout >= 2 )
    dLdx = zeros( 1, Nstates );
    dLdu = u' * Q;
end
