function [ dfdx, dfdu, dfdt ] = linear_instant_cons_pattern( user )

Nstates = user.Nstates;
Ninputs = user.Ninputs;
len_cons = user.len_cons;
pos = user.idxs.pos;

dfdx = ones( len_cons, Nstates );

dfdu = zeros( len_cons, Ninputs );
dfdt = zeros( len_cons, 1 );