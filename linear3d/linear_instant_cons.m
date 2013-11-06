function [ fx, dfdx ] = linear_instant_cons( x, user )

Nstates = user.Nstates;
obstacle = user.obstacle;
len_cons = user.len_cons;
len_obs = length( obstacle );
x_min = user.x_min;
x_max = user.x_max;
bc_idxsmin = user.bc_idxsmin;
bc_idxsmax = user.bc_idxsmax;
len_boxmax = length( bc_idxsmax );
len_boxmin = length( bc_idxsmin );

assert( len_cons == ( len_obs + len_boxmax + len_boxmin ) );

fx = zeros( len_cons, 1 );

if( nargout == 1 )
    for k = 1:len_obs
        fx( k ) = obstacle_cons( x, obstacle{k}, user );
    end
    for k = 1:len_boxmax
        fx( len_obs + k ) = x( bc_idxsmax( k ) ) - x_max( bc_idxsmax( k ) );
    end
    for k = 1:len_boxmin
        fx( len_obs + len_boxmax + k ) =  - x( bc_idxsmin( k ) ) + x_min( bc_idxsmin( k ) );
    end
else
    dfdx = zeros( len_cons, Nstates );
    for k = 1:len_obs
        [ fx( k ), dfdx( k, : ) ] = obstacle_cons( x, obstacle{k}, user );
    end
    for k = 1:len_boxmax
        fx( len_obs + k ) = x( bc_idxsmax( k ) ) - x_max( bc_idxsmax( k ) );
        dfdx( len_obs + k, bc_idxsmax( k ) ) = 1;
    end
    for k = 1:len_boxmin
        fx( len_obs + len_boxmax + k ) =  - x( bc_idxsmin( k ) ) + x_min( bc_idxsmin( k ) );
        dfdx( len_obs + len_boxmax + k, bc_idxsmin( k ) ) = -1;
    end
end


