function [ f, dfdx, dfdu, dfdt ] = needle_dynamics_turn( x, u, t, user )

f = [  0; ...
       0; ...
       0; ...
       0; ...
       0; ...
       u(2); ];

if( nargout >= 2 )
  dfdx = zeros( 6, 6 );
      
  dfdu = zeros( 6, 2 );
  dfdu(6,2) = 1;

  dfdt = zeros( 6, 1 );
end
