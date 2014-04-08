function user = state_variation_precompute( user )

Nsamples = user.Nsamples;
Nstates = user.Nstates;
Ninputs = user.Ninputs;
Nmodes = user.Nmodes;
sys_model = user.functions.sys_model;
instant_cost = user.functions.instant_cost;
terminal_cost = user.functions.terminal_cost;
instant_cons = user.functions.instant_cons;
len_cons = user.len_cons;
higher_order = user.higher_order;

tau = user.data.tau;
u = user.data.u;
d = user.data.d;

% construct all of the data structures that we are planning on saving
x = zeros( Nstates, Nsamples + 1 );
cost = zeros( 1, Nsamples + 1 );
dcost = zeros( 1, Nstates + 1 );
dPhi = zeros( Nstates + 1, Nstates + 1, Nsamples + 1, Nsamples + 1 );
dMdu = zeros( Nstates + 1, Ninputs, Nsamples );
M = zeros( Nstates + 1, Nmodes, Nsamples );
if ( len_cons > 0 )
  Psi = -Inf;
  dhdx = zeros( len_cons, Nstates, Nsamples + 1 );
end

% construct the state and cost for this particular iteration
if( higher_order )
  [ x, cost ] = multistep_integ( user, tau, u, d );
else
  x_cost = fwd_euler( user.x0, tau, u, d, user );
  x( :, : ) = x_cost( 1:Nstates, : );
  cost( :, : ) = x_cost( end, : );
end
[ ~, dphidx ] = terminal_cost( x( 1:Nstates, end ), user );
dcost( 1, 1:Nstates ) = dphidx;
dcost( 1, end ) = 1;

% evaluate the constraint
if ( len_cons > 0 )
  for k = 1:( Nsamples + 1 )
    [ hx, dhdx( :, :, k ) ] = instant_cons( x( 1:Nstates, k ), user );
    Psi = max( [ hx; Psi ] );
  end
end

% precomputing the vector fields used over and over...
dMdx = zeros( Nstates + 1, Nstates + 1, Nsamples + 1 );
for k = 1:Nsamples
  t = tau( k );
  for i = 1:Nmodes
    [ f, dfdx, dfdu ] = sys_model{i}( x( :, k ), u( :, k ), t, user );
    [ L, dLdx, dLdu ] = instant_cost{i}( x( :, k ), u( :, k ), t, user );
    M( :, i, k ) = [ f; L ];
    dMdx( :, 1:Nstates, k ) = dMdx( :, 1:Nstates, k ) + d( i, k ) * [ dfdx; dLdx ];
    dMdu( :, :, k ) = dMdu( :, :, k ) + d( i, k ) * [ dfdu; dLdu ];
  end
end

% compute the STM
for ell = 1:( Nsamples + 1 )
  dPhi( :, :, ell, ell ) = eye( Nstates + 1 );
  for k = ell:Nsamples
    dt = tau( k + 1 ) - tau( k );
    dPhi( :, :, k + 1, ell ) = dPhi( :, :, k, ell ) + dt * dMdx( :, :, k ) * dPhi( :, :, k, ell );
  end
end

% saving the stuff we need to use later
user.data.x = x;
user.data.cost = cost;
user.data.dcost = dcost;
user.data.dPhi = dPhi;
user.data.dMdu = dMdu;
user.data.M = M;
if ( len_cons > 0 )
  user.data.Psi = Psi;
  user.data.dhdx = dhdx;
end
