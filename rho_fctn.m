function [ domnew, unew, dnew ] = rho_fctn( user, N, tau, u, d )
Nmodes = user.Nmodes;
Ninputs = user.Ninputs;

[ ~,   uwav ] = wavelet_transform( N, tau, u );
[ dom, dwav ] = wavelet_transform( N, tau, d );

Nnew = pow2( N );
domnew = zeros( 1, Nnew * Nmodes + 1 );
unew = zeros( Ninputs, Nnew * Nmodes );
dnew = zeros( Nmodes, Nnew * Nmodes );

dt = dom( end ) - dom( 1 );
assert( length( dom ) == Nnew + 1 );
for k = 1:Nnew
  accum = 0;
  for j = 1:Nmodes
    domnew( ( k - 1 ) * Nmodes + j ) = dom( k ) + accum;
    accum = accum + dwav( j, k ) * dt / Nnew;
    unew( :, ( k - 1 ) * Nmodes + j ) = uwav( :, k );
    dnew( j, ( k - 1 ) * Nmodes + j ) = 1;
  end
end
domnew( Nnew * Nmodes + 1 ) = tau(end);
