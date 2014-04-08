function out = optfctn_cons_g( s, Prob )

user = Prob.user;
Nsamples = user.Nsamples;
Nstates = user.Nstates;
len_s = user.len_s;
len_cons = user.len_cons;
len_cineq = user.len_cineq;
uidx = user.idxs.u; % [ Ninputs, Nsamples ]
didx = user.idxs.d; % [ Nmodes, Nsamples ]

% precomputed data
tau = user.data.tau;
u = user.data.u;
d = user.data.d;
dcost = user.data.dcost;
dPhi = user.data.dPhi;
dMdu = user.data.dMdu;
M = user.data.M;

if( len_cons > 0 )
  dhdx = user.data.dhdx;
end

% the vector over which we are optimizing
[ ~, u_p, d_p ] = state_decode( user, s ); % RV: This is a dummy variable, do not confuse with \zeta(\xi,\eta) in the paper.

out = zeros( len_cineq, len_s );

% computing the variation of the state at each instant in time
[ ~, dDxdp ] = state_variation( dPhi, dMdu, M, u_p, d_p, user );

% zeta is present as an upper in all the constraints.
out( :, 1) = repmat( -1, len_cineq, 1);

% derivative of \xi' - \xi function
for k = 1:Nsamples
  dt = tau( k + 1 ) - tau( k );
  out( :, uidx( :, k ) ) = repmat( 2 * dt * ( u_p( :, k ) - u( :, k ) )', len_cineq, 1 );
  out( :, didx( :, k ) ) = repmat( 2 * dt * ( d_p( :, k ) - d( :, k ) )', len_cineq, 1 );
end

% derivative of dJ
out( 1, : ) = out( 1, : ) + dcost * dDxdp( :, :, end );

% derivative of each of the dpsi_{j,t}
if ( len_cons > 0 )
  for k = 2:( Nsamples + 1 )
    for j = 1:len_cons
      idx = ( k - 2 ) * len_cons + j + 1;
      out( idx, : ) = out( idx, : ) + dhdx( j, :, k ) * dDxdp( 1:Nstates, :, k );
    end
  end
end