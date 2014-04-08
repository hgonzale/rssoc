function plot_trajectories( arg1_, arg2_ )

if( nargin == 1 )
  if( ~isfield( arg1_, 'user' ) )
    error( 'Wrong argument.' );
  end
  user = arg1_.user(end);
  x = fwd_euler( user.x0, user.data.tau, user.data.u, user.data.d, user );
  u = user.data.u;
  d = user.data.d;
  tau = user.data.tau;
elseif( nargin == 2 )
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
%instead of terminal loc, I have qxd, qyd, qzd, qpsid 
qxd = user.cost.terminal_cost.qxd;  
qyd = user.cost.terminal_cost.qyd; 
qzd = user.cost.terminal_cost.qzd;
qpsid = user.cost.terminal_cost.qpsid;
qx = user.idxs.qx;
qy = user.idxs.qy;
qz = user.idxs.qz;
%px = user.idxs.px;
%py = user.idxs.py;
%pz = user.idxs.pz;
qtheta = user.idxs.qtheta;
qphi = user.idxs.qphi;
qpsi = user.idxs.qpsi;
%vtheta = user.idxs.vtheta;
%vphi = user.idxs.vphi;
%vpsi = user.idxs.vpsi;

color = { 'b', 'm', [0 1 0], 'y' };
scrsz = get( 0, 'ScreenSize' );

% Trajectory figure
figure('Position',[20 scrsz(4)-700 560 420]);
plot3( x( qx, : ), x( qy, : ), x( qz, : ), '-', 'linewidth', 2, 'Color', 'r' );
hold on;

% Plot the terminal location (represented by qxd, qyd, qzd (and qpsid)) 
plot3( qxd, qyd, qzd, 'ko', ...
  'MarkerFaceColor', [ 0 1 0 ] ,...
  'MarkerSize', 7 );

%axis( [-4, 3, -2, 5, 0, 12] );
%axis equal;
grid on;
%box on;
%set( gca, 'XTick', [-4, 3], 'YTick', [-2, 5], 'ZTick', [0, 12] );
xlabel( 'x_1' );
ylabel( 'x_2' );
zlabel( 'x_3' );
set( get(gca,'XLabel'), 'Position', [-0.5, 5.2, -0.2] );
set( get(gca,'YLabel'), 'Position', [-4.2, 1.5, -0.2] );
set( get(gca,'ZLabel'), 'Position', [ 3.2, 5.2,    6] );
% plot3( linspace(3,0), zeros(100,1), zeros(100,1),'-', 'linewidth', 1, 'Color', 'k' )
% plot3( zeros(100,1), linspace(-3.5,0), zeros(100,1),'-', 'linewidth', 1, 'Color', 'k' )
%view( -110, 15 );
% view(56, 30);
hold off;
% axis([-5.8957    3.0000   -2.4822    4.9961         0   12.0000]);
% axis equal;

time = linspace( t0, tf, Nsamples +  1 );

% Figure of qtheta vs time 
%figure('Angle Position',[20 scrsz(4)-700 560 420]);
figure; 
plot( time, x( qtheta, : ), '-', 'linewidth', 2, 'Color', 'r' );
grid on;
xlabel( 'time' );
ylabel( 'qtheta' );


% Figure of qphi vs time 
%figure('Angle Position',[20 scrsz(4)-700 560 420]);
figure; 
plot( time, x( qphi, : ), '-', 'linewidth', 2, 'Color', 'r' );
grid on;
xlabel( 'time' );
ylabel( 'qphi' );

% Figure of qpsi vs time 
%figure('Angle Position',[20 scrsz(4)-700 560 420]);
figure; 
plot( time, x( qpsi, : ), '-', 'linewidth', 2, 'Color', 'r' );
grid on;
xlabel( 'time' );
ylabel( 'qpsi' );


% draw the continuous inputs
figure;
amin = u_min * ones( size( u ) );
amax = u_max * ones( size( u ) );

hold on;
grid on;
plot( tau( 1:( end - 1 ) ), amin, 'r--' );
plot( tau( 1:( end - 1 ) ), amax, 'r--' );
plot( tau( 1:( end - 1 ) ), u( 1, : ), '.-', 'markersize', 5 );
xlabel( 'time' );
ylabel( 'u' );
hold off;




% draw the discrete inputs
figure;
%[ tau_unique, I ] = unique( tau );
for cnt = 1:Nmodes
  alabel = [ 'd_' num2str(cnt) ];
  subplot( Nmodes, 1, cnt );
  hold on;
  grid on;
  plot( tau( 1:( end - 1 ) ), d( cnt, : ), '.-', 'markersize', 5 );
  %plot( tau_unique( 1:(end - 1 ) ), d( cnt, I( 1:(end-1) ) ), '.-', 'markersize', 5, 'linewidth', 1.5 );
  %xlabel( 'time', 'fontsize', 21 );
  %ylabel( alabel, 'fontsize', 21 );
  hold off;
end

end 