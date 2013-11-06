function aviobj = animate_trajectories( arg1_, arg2_, arg3_ )
% To encode the movie run the following command in brouwer:
% mencoder -ovc lavc -lavcopts vcodec=wmv2 <filename> -o <output>

if( nargin == 2 || isempty( arg3_ ) )
  warning( 'relax:animate_trajectories', 'Filename is empty, will display the animation only.' );
  filename = [];
elseif ( nargin == 3 )
  filename = arg3_;
else
  error( 'Needs 2 or 3 arguments: ( result, idx, <filename> )' );
end

if( ~isfield( arg1_, 'user' ) || ...
    ~isnumeric( arg2_ ) )
  error( 'Some of the arguments are wrong.' );
end

user = arg1_.user( arg2_ );
x = fwd_euler( user.x0, user.data.tau, user.data.u, user.data.d, user );
u = user.data.u;
d = user.data.d;
tau = user.data.tau;
  
t0 = user.t0;
tf = user.tf;
terminal_loc = user.terminal_loc;
px = user.idxs.px;
py = user.idxs.py;
pz = user.idxs.pz;
yaw = user.idxs.yaw;
pitch = user.idxs.pitch;
roll= user.idxs.roll;
obstacle = user.obstacle;

% Important parameters of the movie go here.
animsize = [ 600, 800 ];
fps = 60;
rot_offset = -60;
rot_speed = -15;

figanim = figure( 'position', [ 50, 50, animsize(1), animsize(2) ], 'Color', 'white' );

xlims = minmax( x(px,:) );
ylims = minmax( x(py,:) );
zlims = minmax( x(pz,:) );
lims_tolerance = 5;

dt = 1/fps;
if( ~isempty( filename ) )
  aviobj = avifile( filename, 'fps', fps );
end

for tj = t0:dt:tf
  [ idx, lambda ] = interp_time( tau, tj );
  clf;
  hold on;

  xinterp = ( 1 - lambda ) * x( :, idx ) + lambda * x( :, idx + 1 );
  xnew = [ x(:,1:idx), xinterp ];

  % Draw the path
  plot3( xnew( px, : ), xnew( py, : ), xnew( pz, : ), '-', ...
    'LineWidth', 2, ...
    'Color', 'b' );
  plot3( xinterp( px ), xinterp( py ), xinterp( pz ), 'ko', ...
    'MarkerFaceColor', [ 0 0 1 ] ,...
    'MarkerSize', 5 );

  % Draw the terminal location
  plot3( terminal_loc( px, 1 ), terminal_loc( py, 1 ), terminal_loc( pz, 1 ), 'ko', ...
    'MarkerFaceColor', [ 0 1 0 ] ,...
    'MarkerSize', 7 );

  % Draw vector depending on the current mode
  [ ~, didx ] = max( d(:,idx) );
  rotM = makehgtform( 'xrotate', xinterp( yaw ) ) * ...
    makehgtform( 'yrotate', xinterp( pitch ) ) * ...
    makehgtform( 'zrotate', xinterp( roll ) );
  e1 = -rotM(1:3,1:3) * [1;0;0];
  e2 = rotM(1:3,1:3) * [0;1;0];
  e3 = rotM(1:3,1:3) * [0;0;1];
  switch( didx )
    case{ 1 }
      quiver3( xinterp( px ), ...
        xinterp( py ), ...
        xinterp( pz ), ...
        0.25 * u(1,idx) * e3(1), ...
        0.25 * u(1,idx) * e3(2), ...
        0.25 * u(1,idx) * e3(3), ...
        0, ...
        'LineWidth', 1.5, ...
        'Color', [ 190, 0, 0 ] / 255 );
    case{ 2 }
      quiver3( xinterp( px ), ...
        xinterp( py ), ...
        xinterp( pz ), ...
        0.4 * u(2,idx) * e1(1), ...
        0.4 * u(2,idx) * e1(2), ...
        0.4 * u(2,idx) * e1(3), ...
        0, ...
        'LineWidth', 1.5, ...
        'Color', [ 190, 0, 190 ] / 255 );
      quiver3( xinterp( px ), ...
        xinterp( py ), ...
        xinterp( pz ), ...
        0.4 * u(2,idx) * e2(1), ...
        0.4 * u(2,idx) * e2(2), ...
        0.4 * u(2,idx) * e2(3), ...
        0, ...
        'LineWidth', 1.5, ...
        'Color', [ 190, 0, 190 ] / 255 );
    otherwise
      error( 'Unknown mode' );
  end

  % Draw the obstacles
  for k = 1:length( obstacle )
    switch( obstacle{k}.type )
      case{ 'sphere' }
        r = obstacle{k}.radius;
        xc = obstacle{k}.center;
      
        xlims(1) = min( xlims(1), xc(px) - r );
        xlims(2) = max( xlims(2), xc(px) + r );
        ylims(1) = min( ylims(1), xc(py) - r );
        ylims(2) = min( ylims(2), xc(py) + r );
        zlims(1) = min( zlims(1), xc(pz) - r );
        zlims(2) = min( zlims(2), xc(pz) + r );
      
        [ sx, sy, sz ] = sphere( 30 );
        h = surfl( xc( px ) + r * sx, xc( py ) + r * sy, xc( pz ) + r * sz );
        set( h, 'edgecolor','none' );
        colormap(gray);
        shading interp;
      otherwise
        warning( 'relax:animate_trajectories', 'Unknown obstacle type: %s.', obstacle{k}.type );
    end
  end

  % so that the plot looks reasonable...
  box on;
  % grid on;
  set( gca, 'XTick', -4:3, 'YTick', -2:5, 'ZTick', 0:12 );
  
  % irrational number avoids jumping at multiples of 90
  view( rot_speed * tj  + rot_offset + 1e-3 * sqrt(2), 25 );
  % axis vis3d
  axis( [ xlims(1) - lims_tolerance, xlims(2) + lims_tolerance, ...
          ylims(1) - lims_tolerance, ylims(2) + lims_tolerance, ...
          zlims(1), zlims(2) + lims_tolerance ] );
  axis equal;
  axis manual;

  if( ~isempty( filename ) )
    myframe = getframe( figanim );
    aviobj = addframe( aviobj, myframe );
  else
    pause( dt );
  end
end
if( ~isempty( filename ) )
  aviobj = close( aviobj );
end

end


function [ idx, lambda ] = interp_time( tau, t )
  assert( t >= tau(1) || t <= tau(end) );
  
  if( t == tau(end) )
    idx = length( tau ) - 1;
    lambda = 1;
    return;
  end
  
  idx = find( tau == t, 1, 'first' );
  if( ~isempty( idx ) )
    lambda = 0;
    return;
  end
  
  idx = find( tau < t, 1, 'last' );
  lambda = ( t - tau(idx) ) / ( tau(idx+1) - tau(idx) );
end
