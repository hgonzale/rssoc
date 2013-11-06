function [ h, dhdx ] = obstacle_cons( x, obstacle, user )

Nstates = user.Nstates;
pos = user.idxs.pos;

h = 0;
if( nargout == 2 )
  dhdx = zeros( 1, Nstates );
end
xpos = x(pos);

% if( strcmpi( obstacle.type, 'polytope' ) || strcmpi( obstacle.type, 'rectangle' ) ) % Polytope
if( obstacle.type_idx == 1 )
  A = obstacle.A;
  b = obstacle.b;
  p = obstacle.p;
  maxval = 20/p;
  exparg = A' * xpos - b;
  res = zeros( length( exparg ), 1 );
  for k = 1:length( exparg )
    if( exparg(k) <= maxval )
      res(k) = exp( exparg(k) * p );
    else
      % HG: This equation was found by Maryam and Ram, blame them if it doesn't work :)
      c = exp( maxval ) * A(:,k);
      d = exp( maxval ) * ( 1 - maxval - b(k) );
      res(k) = ( c' * xpos + d )^p;
    end
  end
  h = sum( res );

  if( nargout == 2 )
    for k = 1:length( exparg )
      if( exparg(k) <= maxval )
        dhdx(pos) = dhdx(pos) + A(:,k)' * res(k);
      else
        c = exp( maxval ) * A(:,k);
        d = exp( maxval ) * ( 1 - maxval - b(k) );
        dhdx(pos) = dhdx(pos) + c' * ( c' * xpos + d )^(p-1);
      end
    end
    dhdx(pos) = - 1/h * dhdx(pos);
  end

  h = -log( h^(1/p) );
% elseif( strcmpi( obstacle.type, 'sphere' ) ) % Sphere
elseif( obstacle.type_idx == 2 )
  xc = obstacle.center;
  r = obstacle.radius;
  
  dist = xpos - xc;
  h = r * r - dist' * dist;

  if( nargout == 2 )
    dhdx( pos ) = - 2 * dist;
  end
else
  error( 'Unknown obstacle type: %d.', obstacle.type );
end
