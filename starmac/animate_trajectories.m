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
pz = user.idxs.pz;
theta = user.idxs.theta;
obstacle = user.obstacle;
L = user.starmac_params.length;

% Important parameters of the movie go here.
animsize = [ 1000, 400 ];
fps = 50;

figanim = figure( 'position', [ 50, 50, animsize(1), animsize(2) ], 'Color', 'white' );

xlims = minmax( x(px,:) );
zlims = minmax( x(pz,:) );
zlims(1) = min( -0.1, zlims(1) );
lims_tolerance = 0.5;

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
  
  plot( xnew( px, : ), xnew( pz, : ), 'LineWidth', 2 );

  line( xinterp( px ) - L * cos( xinterp( theta ) ) * [ .5; -.5 ], ...
    xinterp( pz ) + L * sin( xinterp( theta ) ) * [ .5; -.5 ], ...
    'Color', [0 0 0], 'LineWidth', 2 );
  
  [ ~, didx ] = max( d(:,idx) );
  switch( didx )
    case{ 1 }
      quiver( xinterp( px ) - L * cos( xinterp( theta ) ) * 0.5, ...
        xinterp( pz ) + L * sin( xinterp( theta ) ) * 0.5, ...
        0.25 * u(idx) * sin( xinterp( theta ) ), ...
        0.25 * u(idx) * cos( xinterp( theta ) ), ...
        'LineWidth', 1.5, ...
        'Color', [ 100, 100, 100 ] / 255 );
    case{ 2 }
      quiver( xinterp( px ), ...
        xinterp( pz ), ...
        0.25 * u(idx) * sin( xinterp( theta ) ), ...
        0.25 * u(idx) * cos( xinterp( theta ) ), ...
        'LineWidth', 1.5, ...
        'Color', [ 100, 100, 100 ] / 255 );
    case{ 3 }
      quiver( xinterp( px ) + L * cos( xinterp( theta ) ) * 0.5, ...
        xinterp( pz ) - L * sin( xinterp( theta ) ) * 0.5, ...
        0.25 * u(idx) * sin( xinterp( theta ) ), ...
        0.25 * u(idx) * cos( xinterp( theta ) ), ...
        'LineWidth', 1.5, ...
        'Color', [ 100, 100, 100 ] / 255 );
    otherwise
      error( 'Unknown mode' );
  end
  

  % Draw the terminal location
  plot( terminal_loc( px ), terminal_loc( pz ), 'ko', ...
    'MarkerFaceColor', [ 0 1 0 ] , 'MarkerSize', 7 );

  % Draw the obstacles
  for k = 1:length( obstacle )
    if( strcmpi( obstacle{k}.type, 'rectangle' ) )
      rectangle( 'Position', ...
        [ obstacle{k}.x-(obstacle{k}.width)/2, obstacle{k}.y-(obstacle{k}.height)/2, ...
          obstacle{k}.width, obstacle{k}.height ], ...
        'EdgeColor', [139 115 85]/255, ...
        'FaceColor', [139 115 85]/255, ...
        'LineStyle', '-', ...
        'LineWidth', 2 );
    else
      warning( 'mip:starmac:plot_trajectories', 'Unknown obstacle type: %s.', obstacle{k}.type );
    end
  end

  % so that the plot looks reasonable...
  xlabel( 'z', 'FontSize', 14 );
  ylabel( 'x', 'FontSize', 14 );
  % set( get( gca, 'XLabel' ), 'Position', [ 2 -1.05 ] )
  % set( get( gca, 'YLabel' ), 'Position', [ -1.05 2.5 ] )
  set( gca, 'XTick', -20:0.5:20, 'YTick', -1:0.5:10 );

  axis( [ xlims(1) - lims_tolerance, xlims(2) + lims_tolerance, ...
          zlims(1) - lims_tolerance, zlims(2) + lims_tolerance ] );
  axis equal;
  axis manual;

  % This fucking rectangle fucks up the axis, so we fix the axis limits before plotting it.
  rectangle( 'Position', [-2 -5 17 5], ...
    'EdgeColor', [139 115 85]/255, ...
    'FaceColor', [139 115 85]/255, ...
    'LineStyle', '-', ...
    'LineWidth', 2 );

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
