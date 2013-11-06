function  [ F, c, qc ] = encode_QP_cons( user )

Nsamples = user.Nsamples;
Nstates = user.Nstates;
Ninputs = user.Ninputs;
Nmodes = user.Nmodes;
len_s = user.len_s;
len_cons = user.len_cons;
len_cineq = user.len_cineq;
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
sep = 0;
dummy1 = 1:Nmodes;
dummy2 = 1:Ninputs;
for k = 1:Nsamples
    dt = tau( k + 1 ) - tau( k );
    deltaTd( dummy1 ) = dt;
    deltaTu( dummy2 ) = dt;
    sep = sep + dt * ( u( :, k )' * u( :, k ) + d( :, k )' * d( :, k ) );
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

% objective function ( there is a linear component and thats it )
F = zeros( len_s, len_s );
c = zeros( len_s, 1 );
c( 1 ) = 1;

% inequality constraints
qc = struct;

% all of the inequality constraints include ( xi_p )^2
dummy = zeros( len_s, 1 );
dummy( [ uidx( : ); didx( : ) ] ) = [ deltaTu; deltaTd ];
for i = 1:len_cineq
    qc( i ).Q = sparse( diag( dummy ) );
    qc( i ).a = zeros( len_s, 1 ); % just initializing all of them...
    qc( i ).a( 1 ) = -1; % this is the dummy variable zeta ...
    qc( i ).a( didx( : ) ) = -2 * d( : )' .* deltaTd'; % the portion due to the 2 * d
    qc( i ).a( uidx( : ) ) = -2 * u( : )' .* deltaTu'; % the portion due to the 2 * u
end

% the constraint due to DJ
% xi_p^2 + DJ - 2xi'*xi_p - zeta <= xi^2
% remember the DJ component multiplies (xi_p - xi)
qc( 1 ).a = qc( 1 ).a + ( dcost * dDxdp( :, :, end ) )'; % the portion actually due to DJ multiplied by up and dp
qc( 1 ).r_U = dcost * dDxdp( :, [ uidx( : ); didx( : ) ], end ) * [ u( : ); d( : ) ] - sep;
if ( len_cons > 0 )
    if( Psi > 0 )
        qc( 1 ).r_U = qc( 1 ).r_U + Psi;
    end
end

% the constraint due to each of the Dpsi
% xi_p^2 + Dpsi - 2xi'*xi_p - zeta <= xi^2
% remember the Dpsi component multiplies (xi_p - xi)
if ( len_cons > 0 )
    for k = 1:( Nsamples + 1 )
        for j = 1:len_cons
            idx = ( k - 1 ) * len_cons + j + 1;
            qc( idx ).a = qc( idx ).a + ( dhdx( j, :, k ) * dDxdp( 1:Nstates, :, k ) )'; % the portion actually due to Dpsi multiplied by up and dp
            qc( idx ).r_U = dhdx( j, :, k ) * dDxdp( 1:Nstates, [ uidx( : ); didx( : ) ], k ) * [ u( : ); d( : ) ] - sep;
            if ( Psi <= 0 )
                qc( idx ).r_U = qc( idx ).r_U - ( gamma * Psi );
            end
        end
    end
end

% error check
assert( length( qc ) == len_cineq );

