function s = state_encode( user, zeta, u, d )

Nsamples = user.Nsamples;
Nmodes = user.Nmodes;
Ninputs = user.Ninputs;

zidx = user.idxs.zeta; % [ 1 ] % RV: This is a dummy variable, do not confuse with \zeta(\xi,\eta) in the paper.
uidx = user.idxs.u; % [ Ninputs, Nsamples + 1 ]
didx = user.idxs.d; % [ Nmodes, Nsamples + 1 ]
len_s = user.len_s;

s = zeros( 1, len_s );

% Encoding m
s( zidx ) = zeta;

% Encoding u
sizeu = size( u );
if( length( sizeu ) == 2 && all( sizeu >= [ Ninputs, Nsamples ] ) )
    s( uidx ) =  u( 1:Ninputs, 1:( Nsamples ) );
elseif( sizeu( 1 ) >= Ninputs )
    s( uidx ) = repmat( u( 1:Ninputs, 1 ), [ 1, Nsamples ] );
else
    error( 'Unknown state format.' );
end

% Encoding d
sized = size( d );
if( length( sized ) == 2 && all( sized >= [ Nmodes, Nsamples ] ) )
    s( didx ) =  d( 1:Nmodes, 1:( Nsamples ) );
elseif( sized( 1 ) >= Nmodes )
    s( didx ) = repmat( d( 1:Nmodes, 1 ), [ 1, Nsamples ] );
else
    error( 'Unknown state format.' );
end

% one last check...
assert( length( s ) == len_s );
