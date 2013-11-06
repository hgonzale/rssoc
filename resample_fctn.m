function [ domnew, unew, dnew ] = resample_fctn( user, tau, u, d )

Nmodes = user.Nmodes;
Ninputs = user.Ninputs;
tf = user.tf;

N = length( tau ) - 1;
Nnew = N * 2;
domnew = zeros( 1, Nnew + 1 );
unew = zeros( Ninputs, Nnew );
dnew = zeros( Nmodes, Nnew );

difftau = diff( tau );
for k = 1:N
  domnew( 2 * ( k - 1 ) + 1 ) = tau( k );
  domnew( 2 * ( k - 1 ) + 2 ) = tau( k ) + 0.5 * difftau( k );
  unew( :, 2 * ( k - 1 ) + 1 ) = u( :, k );
  unew( :, 2 * ( k - 1 ) + 2 ) = u( :, k );
  dnew( :, 2 * ( k - 1 ) + 1 ) = d( :, k );
  dnew( :, 2 * ( k - 1 ) + 2 ) = d( :, k );
end
domnew( Nnew + 1 ) = tf;
