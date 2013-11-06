function plot_trajectories( arg1_, arg2_ )

if( nargin == 2 )
  if( ~isfield( arg1_, 'user' ) )
    error( 'Wrong argument.' );
  end
  user = arg1_.user( arg2_ );
  x = fwd_euler( user.x0, user.data.tau, user.data.u, user.data.d, user );
  u = user.data.u;
  d = user.data.d;
  tau = user.data.tau;
else
  error( 'Unknown number of arguments.' );
end

Ninputs = user.Ninputs;
Nmodes = user.Nmodes;
Nsamples = user.Nsamples;
u_min = user.u_min;
u_max = user.u_max;
t0 = user.t0;
tf = user.tf;
terminal_loc = user.terminal_loc;
obstacle = user.obstacle;
px = user.idxs.px;
py = user.idxs.py;
pz = user.idxs.pz;
color = { 'b', 'm', [0 1 0], 'y' };
scrsz = get( 0, 'ScreenSize' );

% Trajectory figure
figure('Position',[20 scrsz(4)-700 560 420]);
plot3( x( px, : ), x( py, : ), x( pz, : ), '-', 'linewidth', 2, 'Color', 'r' );
hold on;

% Plot the terminal location
plot3( terminal_loc( px, 1 ), terminal_loc( py, 1 ), terminal_loc( pz, 1 ), 'ko', ...
  'MarkerFaceColor', [ 0 1 0 ] ,...
  'MarkerSize', 7 );

axis([-6     4    -5     5     0    12]);
for k = 1:length(obstacle)
  if( strcmpi( obstacle{k}.type, 'rectangle' ) )
    warning( todo:unimplemented, 'Tell Humberto to program this routine' );
  elseif( strcmpi( obstacle{k}.type, 'polytope' ) )
    warning( todo:unimplemented, 'Tell Humberto to program this routine' );
    
  elseif( strcmpi( obstacle{k}.type, 'sphere' ) )
    r = obstacle{k}.radius;
    xc = obstacle{k}.center;
    
    [sx,sy,sz] = sphere(30);
    h = surfl( xc(px) + r * sx, xc(py) + r * sy, xc(pz) + r * sz);
    set(h, 'edgecolor','none');
    %     shading interp;
    colormap(gray);
    axis([-6    3.0000   -2.5   5  0  12.0000])
    %     shading interp;
  else
    warning( todo:unknown_type, 'Unknown obstacle type: %s.', obstacle{k}.type );
  end
end

axis( [-4, 3, -2, 5, 0, 12] );
axis equal;
grid off;
box on;
set( gca, 'XTick', [-4, 3], 'YTick', [-2, 5], 'ZTick', [0, 12] );
xlabel( 'x_1' );
ylabel( 'x_2' );
zlabel( 'x_3' );
set( get(gca,'XLabel'), 'Position', [-0.5, 5.2, -0.2] );
set( get(gca,'YLabel'), 'Position', [-4.2, 1.5, -0.2] );
set( get(gca,'ZLabel'), 'Position', [ 3.2, 5.2,    6] );
% plot3( linspace(3,0), zeros(100,1), zeros(100,1),'-', 'linewidth', 1, 'Color', 'k' )
% plot3( zeros(100,1), linspace(-3.5,0), zeros(100,1),'-', 'linewidth', 1, 'Color', 'k' )
view( -110, 15 );
% view(56, 30);
hold off;
% axis([-5.8957    3.0000   -2.4822    4.9961         0   12.0000]);
% axis equal;


time = linspace( t0, tf, Nsamples +  1 );

% Draw the continuous inputs
figure;
for cnt = 1:Ninputs
  amin = u_min( cnt ) * ones( size( u ) );
  amax = u_max( cnt ) * ones( size( u ) );
  alabel = [ 'u_' num2str(cnt) ];
  
  subplot( Ninputs, 1, cnt );
  hold on;
  grid on;
  plot( tau( 1:( end - 1 ) ), amin, 'r--' );
  plot( tau( 1:( end - 1 ) ), amax, 'r--' );
%   for i = 1:( length( tau ) - 1 )
%       if ( tau( i ) == tau ( i + 1 ) ) 
%           continue;
%       end
%       [ ~, didx ] = max( d( :, i ) );
%       plot( tau( i ), u( cnt, i ), 'MarkerSize', 5, 'Color', color{ didx } );
%   end
  
  plot( tau( 1:( end - 1 ) ), u( cnt, : ), '.-', 'MarkerSize', 5 );
  xlabel( 'Time' );
  ylabel( alabel );
  hold off;
end

% Draw the discrete inputs
figure;
[ tau_unique, I ] = unique( tau );
for cnt = 1:Nmodes
%   alabel = [ 'd_' num2str(cnt) ];
  if ( cnt == 1 )
      alabel = 'd_1';
  else
      alabel = 'd_2';
  end
  subplot( Nmodes, 1, cnt );
  hold on;
  grid on;
%   plot( tau( 1:( end - 1 ) ), d( cnt, : ), '.-', 'MarkerSize', 5 );
  plot( tau_unique( 1:(end - 1 ) ), d( cnt, I( 1:(end-1) ) ), 'r.', 'MarkerSize', 10, 'LineWidth', 1.5 );
  % xlabel( 'Time', 'FontSize', 21 );
  % ylabel( alabel, 'FontSize', 21 );
  xlabel( 'Time [s]' );
  ylabel( alabel );
  hold off;
end
