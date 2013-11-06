function cost = obj_fctn( user, tau, u, d, x )

x0 = user.x0;
terminal_cost = user.functions.terminal_cost;
states_idxs = user.idxs.states;
cost_idx = user.idxs.cost;

N = length( d(1,:) );

if( nargin == 4 )
  x = fwd_euler( x0, tau, u, d, user, N );
elseif( nargin ~=5 )
  error( 'Number of arguments is incorrect.' );
end

cost = terminal_cost( x( states_idxs, N + 1 ), user ) + x( cost_idx, N + 1 );
