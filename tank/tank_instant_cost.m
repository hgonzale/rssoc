function [ L, dLdx, dLdu, dLdt ] = tank_instant_cost( x, u, t, user )

Nstates = user.Nstates;
Ninputs = user.Ninputs;

L = 2 * ( x( 2 ) - 3 )^2;

if( nargout >= 2 )
  dLdx = zeros( 1, Nstates );
  dLdx( 1, 2 ) = 4 * ( x( 2 ) - 3 );
  
  dLdu = zeros( 1, Ninputs );
  dLdt = 0;
end
