function x = fwd_euler( x0, tau, u, d, user, N )

if( nargin < 6 )
    N = user.Nsamples;
end

Nmodes = user.Nmodes;
Nstates = user.Nstates;

sys_model = user.functions.sys_model;
instant_cost = user.functions.instant_cost;

x = zeros( Nstates + 1, N + 1 );
x( :, 1 ) = [ x0; 0 ];

for k = 1:N
    t = tau( k );
    dt = tau( k + 1 ) - tau( k );
    for j = 1:Nmodes;
        x( :, k + 1 ) = x( :, k + 1 ) + dt * d( j, k ) * ...
            [ sys_model{j}( x( 1:Nstates, k ), u( :, k ), t, user ); instant_cost{j}( x( 1:Nstates, k ), u( :, k ), t, user ) ];
    end
    x( : , k + 1 ) = x( :, k + 1 ) + x( :, k );
end    

% [ state, cost ] = multistep_integ( user, tau, u, d );
% 
% x = [ state; cost ];