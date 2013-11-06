function Psi = max_cons_fctn( user, tau, u, d, x )

x0 = user.x0;
len_cons = user.len_cons;
instant_cons = user.functions.instant_cons;
states_idxs = user.idxs.states;

N = length( d(1,:) );

if( nargin == 4 )
  x = fwd_euler( x0, tau, u, d, user, N );
elseif( nargin ~=5 )
  error( 'Number of arguments is incorrect.' );
end

if ( len_cons > 0 )
  hx = zeros( len_cons * ( N + 1 ), 1 );

  for k = 1:( N + 1 )
    arr = ( 1:len_cons ) + (k-1) * len_cons;
    hx( arr ) = instant_cons( x( states_idxs, k ), user );
  end
  Psi = max( hx );
else
  warning( 'len_cons equals zero, so you should not be calling this function.' );
  Psi = 0;
end
