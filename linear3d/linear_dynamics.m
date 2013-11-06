function [ f, dfdx, dfdu ] = linear_dynamics( x, u, t, A, B, user )

f = A*x + B*u;

if( nargout > 1 )
  dfdx = A;
  dfdu = B;
end
