function [ zeta, u, d ] = state_decode( user, s )

Nsamples = user.Nsamples;
Nmodes = user.Nmodes;
Ninputs = user.Ninputs;

zidx = user.idxs.zeta; % RV: This is a dummy variable, do not confuse with \zeta(\xi,\eta) in the paper.
uidx = user.idxs.u;
didx = user.idxs.d;

zeta = reshape( s( zidx(:) ), 1, 1 );
u = reshape( s( uidx(:) ), [ Ninputs, Nsamples ] );
d = reshape( s( didx(:) ), [ Nmodes, Nsamples ] );
  
