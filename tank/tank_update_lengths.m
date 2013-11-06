function user = tank_update_lengths( user )

Nstates = user.Nstates;
obstacle = user.obstacle;

dummy = 1:Nstates;
user.bc_idxsmax = dummy( isfinite( user.x_max( 1:Nstates ) ) );
user.bc_idxsmin = dummy( isfinite( user.x_min( 1:Nstates ) ) );
user.len_cons = length( obstacle ) + length( user.bc_idxsmax ) + length( user.bc_idxsmin );
