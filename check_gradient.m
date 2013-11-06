function [ grad, num_grad ] = check_gradient( f, g, x, st )

if( isfield( st, 'Nsamples' ) && ~isfield( st, 'user' ) ) % i.e. if the struct is actually 'user'
  Prob = struct;
  Prob.user = st;
elseif( isfield( st, 'user' ) ) % i.e. if the struct is 'Prob'
  Prob = st;
else
  error( 'Invalid struct.' );
end


test = f( x, Prob );
len_f = length( test );
len_x = length( x );
h = 1e-6;

num_grad = zeros( len_f, len_x );

for k = 1:len_x,
  hvec = zeros( size(x) );
  hvec(k) = h;
  fph = f( x + hvec, Prob );
  fmh = f( x - hvec, Prob );
  if( any( isnan( fph ) ) || any( isnan( fmh ) ) )
    error( 'NaN found at k = %d', num2str( k ) );
  end
  num_grad(:,k) = ( fph - fmh ) / ( 2 * h );
end

grad = g( x, Prob );

if( ~all( size( grad ) == [ len_f, len_x ] ) )
  warning( 'relax:check_gradient', ...
    'Matrix sizes do not match.\n\tsize( grad ) = [ %d %d ]\n\tsize( num_grad ) = [ %d %d ].', ...
    length( grad(:,1) ), length( grad(1,:) ), ...
    length( num_grad(:,1) ), length( num_grad(1,:) ) );
  return;
end

fprintf( 1, 'Maximum error = %d\n', max( abs( grad(:) - num_grad(:) ) ) );
if( max( abs( grad(:) - num_grad(:) ) ) > 1e-4 )
  spy( abs( grad - num_grad ) > 1e-4, 5 );
end
