function user = create_idxs_struct( user )

Nsamples = user.Nsamples;
Nstates = user.Nstates;
Nmodes = user.Nmodes;
Ninputs = user.Ninputs;

user.idxs.states = 1:Nstates;
user.idxs.cost = Nstates + 1;

offset = 0;

user.idxs.zeta = offset + 1;
offset = offset + length( user.idxs.zeta(:) ); % RV: This is a dummy variable, do not confuse with \zeta(\xi,\eta) in the paper.

user.idxs.u = offset + ...
    reshape( 1:( Ninputs * ( Nsamples ) ), [ Ninputs, Nsamples ] );
offset = offset + length( user.idxs.u(:) );

user.idxs.d = offset + ...
    reshape( 1:( Nmodes * ( Nsamples ) ), [ Nmodes, Nsamples ] );
