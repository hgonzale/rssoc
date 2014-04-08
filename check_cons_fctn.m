function [ grad, num_grad ] = check_cons_fctn( user, tau, u, d )

if( nargin == 1 )
  x0 = [user.data.u(:); user.data.d(:)];
else
  x0 = [ u(:); d(:) ];
  user.data.tau = tau;
end

[ grad, num_grad ] = check_gradient( @my_f, @my_g, x0, user );

end

function out = my_f( s, Prob )

Ninputs = Prob.user.Ninputs;
Nmodes = Prob.user.Nmodes;
Nsamples = Prob.user.Nsamples;
instant_cons = Prob.user.functions.instant_cons;
len_cineq = Prob.user.len_cineq;
len_cons = Prob.user.len_cons;

u = reshape( s( 1:( Ninputs * Nsamples ) ), Ninputs, Nsamples );
d = reshape( s( ( Ninputs * Nsamples + 1 ):end ), Nmodes, Nsamples );

x = fwd_euler( Prob.user.x0, Prob.user.data.tau, u, d, Prob.user );

out = zeros( len_cineq - 1, 1 ); % no dJ and no zeta ( hence its -1 in row and column )

for k = 2:( Nsamples + 1 )
  idx = ( k - 2 ) * len_cons + 1;
  out( idx:( idx + len_cons - 1 ) ) = instant_cons( x( :, k ), Prob.user );
end

end

function out = my_g( x, Prob )

user = Prob.user;
Ninputs = user.Ninputs;
Nmodes = user.Nmodes;
Nsamples = user.Nsamples;

user.data.u = reshape( x( 1:( Ninputs * Nsamples ) ), Ninputs, Nsamples );
user.data.d = reshape( x( ( Ninputs * Nsamples + 1 ):end ), Nmodes, Nsamples );

user = state_variation_precompute( user );
[ ~, ~, A, ~, ~ ] = encode_LP_cons( user );

out = A( Nsamples + 2:end, 2:end ); % Ram says this is the derivative.

end