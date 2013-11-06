function out = optfctn_cons_f( s, Prob )

user = Prob.user;
Nsamples = user.Nsamples;
Nstates = user.Nstates;
len_cons = user.len_cons;
len_cineq = user.len_cineq;
gamma = user.optfctn.gamma;

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

% the vector over which we are optimizing
[ zeta, u_p, d_p ] = state_decode( user, s ); % RV: This is a dummy variable, do not confuse with \zeta(\xi,\eta) in the paper.

out = zeros( len_cineq, 1 );

% the distance between \xi' and \xi
sep = 0;
for k = 1:Nsamples
    dt = tau( k + 1 ) - tau( k );
    sep = sep + dt * ( ( u_p( :, k ) - u( :, k ) )' * ( u_p( :, k ) - u( :, k ) ) + ...
        ( d_p( :, k ) - d( :, k ) )' * ( d_p( :, k ) - d( :, k ) ) );
end

% computing the variation of the state at each instant in time
Dx = zeros( Nstates + 1, Nsamples + 1 );
Dx( : , : ) = state_variation( dPhi, dMdu, M, u_p, d_p, user );

% DJ
if ( len_cons > 0 )
    if ( Psi > 0 )
        out( 1 ) = -Psi;
    end
end
out( 1 ) = out( 1 ) + dcost * Dx( :, end ) + sep - zeta;

% each of the d Psis
if ( len_cons > 0 )
    for k = 1:( Nsamples + 1 )
        for j = 1:len_cons
            idx = ( k - 1 ) * len_cons + j + 1;
            if ( Psi <= 0 )
                out( idx ) = gamma * Psi;
            end
            out( idx ) = out( idx ) + dhdx( j, :, k ) * Dx( 1:Nstates, k ) + sep - zeta;
        end
    end
end
