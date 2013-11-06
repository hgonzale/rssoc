function [ Dx, dDxdp ] = state_variation( dPhi, dMdu, M, u_p, d_p, user )

Nsamples = user.Nsamples;
Nstates = user.Nstates;
Nmodes = user.Nmodes;

tau = user.data.tau;
u = user.data.u;
d = user.data.d;

Dx = zeros( Nstates + 1, Nsamples + 1 );

% computing dxt
for k = 1:( Nsamples + 1 )
    for l = 1:( k - 1 )
        dt = tau( l + 1 ) - tau( l );
        for i = 1:Nmodes
            Dx( :, k ) = Dx( :, k ) + dt * dPhi( :, :, k, l + 1 ) * M ( :, i, l ) * ( d_p( i, l ) - d( i, l ) );
        end
        Dx( :, k ) = Dx( :, k ) + dt * dPhi( :, :, k, l + 1 ) * dMdu( :, :, l ) * ( u_p( :, l ) - u( :, l ) );
    end
end

if( nargout > 1 )
    dDxdp = zeros( Nstates + 1, user.len_s, Nsamples + 1 );
    uidx = user.idxs.u;
    didx = user.idxs.d;
    
    for k = 1:( Nsamples + 1 )
        for l = 1:( k - 1 )
            dt = tau( l + 1 ) - tau( l );
            for i = 1:Nmodes
                dDxdp( :, didx( i, l ), k ) = dt * dPhi( :, :, k, l + 1 ) * M( :, i, l );
            end
            dDxdp( :, uidx( :, l ), k ) = dt * dPhi( :, :, k, l + 1 ) * dMdu( :, :, l ); 
        end
    end
    
end


