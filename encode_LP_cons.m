function [ F, c, A, b_L, b_U ] = encode_LP_cons( user )

Nsamples = user.Nsamples;
Nstates = user.Nstates;
Ninputs = user.Ninputs;
Nmodes = user.Nmodes;
len_s = user.len_s;
len_cons = user.len_cons;
len_cineq = user.len_cineq;
delta = user.optfctn.delta;
gamma = user.optfctn.gamma;
uidx = user.idxs.u; % [ Ninputs, Nsamples + 1 ]
didx = user.idxs.d; % [ Nmodes, Nsamples + 1 ]

% precomputed data
tau = user.data.tau;
u = user.data.u;
d = user.data.d;
dcost = user.data.dcost;
dPhi = user.data.dPhi;
dMdu = user.data.dMdu;
M = user.data.M;

if ( len_cons > 0 )
    Psi = user.data.Psi;
    dhdx = user.data.dhdx;
end

% computing a large array that we use to fill in everything ( compare with
% the construction in state_variation.m )
deltaTd = zeros( length( didx( : ) ), 1 );
deltaTu = zeros( length( uidx( : ) ), 1 );
dummy1 = 1:Nmodes;
dummy2 = 1:Ninputs;
for k = 1:Nsamples
    dt = tau( k + 1 ) - tau( k );
    deltaTd( dummy1 ) = dt;
    deltaTu( dummy2 ) = dt;
    dummy1 = dummy1 + Nmodes;
    dummy2 = dummy2 + Ninputs;
end

dDxdp = zeros( Nstates + 1, user.len_s, Nsamples + 1 );
for k = 1:( Nsamples + 1 )
    for l = 1:( k - 1 )
        dt = tau( l + 1 ) - tau( l );
        for i = 1:Nmodes
            dDxdp( :, didx( i, l ), k ) = dt * dPhi( :, :, k, l + 1 ) * M( :, i, l );
        end
        dDxdp( :, uidx( :, l ), k ) = dt * dPhi( :, :, k, l + 1 ) * dMdu( :, :, l );
    end
end

% objective function
dummy = zeros( len_s, 1 );
dummy( [ uidx( : ); didx( : ) ] ) = [ deltaTu; deltaTd ];
F = 2 * delta * diag( dummy ); % [ len_s len_s ] (multiplied by two because qpsolve assumes 1/2 x' * F * x )

c = zeros( len_s, 1 );
c( 1 ) = 1; % slack variable

% constraint functions
A = zeros( Nsamples + len_cineq, len_s );
b_L = -Inf * ones( Nsamples + len_cineq, 1 );
b_U = zeros( Nsamples + len_cineq, 1 );

% discrete input constraints
[ A( 1:Nsamples, : ), ~, ~ ] = dinput_encode_lin_cons( user );
b_L( 1:Nsamples ) = zeros( Nsamples, 1 );

% portion of the constraint that is identical for all the remaining
% constraints
for idx = ( Nsamples + 1 ):( len_cineq + Nsamples )
    A( idx, 1 ) = -1; % this is the dummy variable zeta
end

% constraint due to DJ
idx = Nsamples + 1;
A( idx, : ) = A( idx, : ) + ( dcost * dDxdp( :, :, end ) ); % the portion actually due to DJ multiplied by up and dp
if ( len_cons > 0 )
    if( Psi > 0 )
        b_U( idx ) = b_U( idx ) + Psi;
    end
end

% constraint due to each of the dPsi
if ( len_cons > 0 )
    for k = 1:( Nsamples + 1 )
        for j = 1:len_cons
            idx = ( k - 1 ) * len_cons + j + Nsamples + 1;
            A( idx, : ) = A( idx, : ) + ( dhdx( j, :, k ) * dDxdp( 1:Nstates, :, k ) ); % the portion actually due to Dpsi DJ multiplied by up and dp
            if ( Psi <= 0 )
                b_U( idx ) = b_U( idx ) - ( gamma * Psi );
            end
        end
    end
end
