function user = update_user( user, iter_Nsamples, iter_tau, iter_u, iter_d )

if ( nargin > 1 )
    user.Nsamples = iter_Nsamples;
    user.data.tau = iter_tau;
    user.data.u = iter_u;
    user.data.d = iter_d;
end

Nsamples = user.Nsamples;
Nmodes = user.Nmodes;
Ninputs = user.Ninputs;
update_lengths = user.functions.update_lengths;

% updates the number of user constraints at each instance in time (i.e.
% user.len_cons keep in mind that the state constraints are included
% in this count)
user = update_lengths( user );
% the size of the optimization variables ( discrete inputs, continuous
% inputs and the minmax )
user.len_s = ( Nsamples ) * ( Nmodes + Ninputs ) + 1 ;
% the number of inequality constraints ( each of the constraints, without counting the initial condition, and the
% cost )
user.len_cineq = Nsamples * user.len_cons + 1;

user = create_idxs_struct( user ); % Create idxs struct

assert( length( user.idxs.zeta(:) ) + length( user.idxs.d(:) ) + ...
    length( user.idxs.u(:) ) == user.len_s );

% Modify 'rectangle' obstacles to make them compatible as polytopes.
for k = 1:length( user.obstacle )
    if( strcmpi( user.obstacle{k}.type, 'rectangle' ) == 1 )
        x = user.obstacle{k}.x;
        y = user.obstacle{k}.y;
        width = user.obstacle{k}.width;
        height = user.obstacle{k}.height;
        
        user.obstacle{k}.A = [ 0, 1, 0, -1; 1, 0, -1, 0 ];
        user.obstacle{k}.b = [ y + height/2; x + width/2; - y + height/2; - x + width/2 ];
        user.obstacle{k}.type_idx = 1;
    elseif( strcmpi( user.obstacle{k}.type, 'polytope' ) == 1 )
        user.obstacle{k}.type_idx = 1;
    elseif( strcmpi( user.obstacle{k}.type, 'sphere' ) == 1 )
        user.obstacle{k}.type_idx = 2;
    else
        warning( 'relax:update_user', 'Unknown obstacle type (%s). The obstacle is left without a type index.', user.obstacle{k}.type );
    end
end
