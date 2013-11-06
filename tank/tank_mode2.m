function [ f, dfdx, dfdu, dfdt ] = tank_mode2( x, u, t, user )

Nstates = user.Nstates;
Ninputs = user.Ninputs;

f = [ 2 - sqrt( x( 1 ) ); sqrt( x( 1 ) ) - sqrt( x( 2 ) ) ];

if( nargout >= 2 )
  dfdx = zeros( Nstates, Nstates );
  
  dfdx( 1, 1 ) = -0.5 * x( 1 ).^( -0.5 );
  dfdx( 2, 1 ) = 0.5 * x( 1 ).^( -0.5 );
  dfdx( 2, 2 ) = -0.5 * x( 2 ).^( -0.5 );
      
  dfdu = zeros( Nstates, Ninputs );
  dfdt = zeros( Nstates, 1 );
end
