function [ dfdx, dfdu, dfdt ] = tank_instant_cons_pattern( user )

Nstates = user.Nstates;
Ninputs = user.Ninputs;
len_cons = user.len_cons;
pos = user.idxs.pos;

dfdx = zeros( len_cons, Nstates );
dfdx( :, pos ) = ones( len_cons, length( pos ) );

dfdu = zeros( len_cons, Ninputs );
dfdt = zeros( len_cons, 1 );