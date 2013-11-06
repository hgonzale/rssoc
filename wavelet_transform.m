function [ dom, val ] = wavelet_transform( N, tau, signal ) 

sigsize = size( signal );

assert( sigsize(2) > 1 ); % Make sure we don't have a transposed vector.
assert( length( tau ) == sigsize(2) + 1 );
% assert( tau(1) == 0 );
% assert( tau(end) == 1 );

t0 = tau(1);
tf = tau(end);
tau = ( tau - t0 ) / ( tf - t0 );

Nsamples = pow2( N );
dom = linspace( 0, 1, Nsamples + 1 );
val = zeros( sigsize(1), Nsamples );

for i = 1:sigsize(1)
  val( i, : ) = dot( diff( tau ), signal( i, : ) );
end

for k = 0:( N - 1 )
  kpow = pow2(k);
  for j = 0:( kpow - 1 )
    posinterval = [ j / kpow, ( j + 0.5 ) / kpow ];
    [ mydom, myval ] = eval_signal( posinterval, tau, signal );
    posarr = ( posinterval(1) * Nsamples + 1 ):( posinterval(2) * Nsamples );
    posval = kpow * myval * diff( mydom )';

    neginterval = [ ( j + 0.5 ) / kpow, ( j + 1 ) / kpow ];
    [ mydom, myval ] = eval_signal( neginterval, tau, signal );
    negarr = ( neginterval(1) * Nsamples + 1 ):( neginterval(2) * Nsamples );
    negval = kpow * myval * diff( mydom )';
    
    tmp =  repmat( posval - negval, 1, length( posarr ) );
    val( :, posarr ) = val( :, posarr ) + tmp; 
    val( :, negarr ) = val( :, negarr ) - tmp;
  end
end

dom = dom * ( tf - t0 ) + t0;

end


function [ dom, val ] = eval_signal( interval, tau, signal )
% assert( length( interval ) == 2 );
% assert( interval(1) < interval(2) );
% assert( interval(1) >= 0 );
% assert( interval(2) <= 1 );

idxs_before = ( tau < interval(1) );
idxs_in = ( interval(1) <= tau ) & ( tau < interval(2) );

if( ~any( idxs_in ) )
  k = find( idxs_before, 1, 'last' );
  % assert( ~isempty(k) ); % HG: I really hope that this assert is never triggered.
  
  dom = interval;
  val = signal( :, k );
else % idxs_in has something inside.
  k = find( idxs_in, 1, 'first' );
  if( interval(1) < tau(k) )
    k = find( idxs_before, 1, 'last' );
    % assert( ~isempty(k) ); % HG: I really hope that this assert is never triggered.

    dom = [ interval(1), tau( idxs_in ), interval(2) ];
    val = [ signal( :, k ), signal( :, idxs_in ) ];
  else
    dom = [ tau( idxs_in ), interval(2) ];
    val = [ signal( :, idxs_in ) ];
  end
end

end
