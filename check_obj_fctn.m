function [ grad, num_grad ] = check_obj_fctn( user, tau, u, d )

if( nargin == 1 )
  x0 = [user.data.u(:); user.data.d(:)];
else
  x0 = [ u(:); d(:) ];
  user.data.tau = tau;
end

[ grad, num_grad ] = check_gradient( @my_f, @my_g, x0, user );

end

function out = my_f( x, Prob )

Ninputs = Prob.user.Ninputs;
Nmodes = Prob.user.Nmodes;
Nsamples = Prob.user.Nsamples;

u = reshape( x( 1:( Ninputs * Nsamples ) ), Ninputs, Nsamples );
d = reshape( x( ( Ninputs * Nsamples + 1 ):end ), Nmodes, Nsamples );

out = obj_fctn( Prob.user, Prob.user.data.tau, u, d );

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

out = A( Nsamples + 1, 2:end ); % Ram says this is the derivative.

end