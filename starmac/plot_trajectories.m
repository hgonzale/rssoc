function plot_trajectories( arg1_, arg2_, arg3_ )

if ( nargin == 2 )
    if( ~isfield( arg1_, 'user' ) )
        error( 'Wrong argument.' );
    end
    user = arg1_.user( arg2_ );
    x = fwd_euler( user.x0, user.data.tau, user.data.u, user.data.d, user );
    u = user.data.u;
    d = user.data.d;
    tau = user.data.tau;
    draw_nominal_trajectory = 0;
elseif ( nargin == 3 )
    if( ~isfield( arg1_, 'user' ) )
        error( 'Wrong argument.' );
    end
    user = arg1_.user( arg2_ );
    x = fwd_euler( user.x0, user.data.tau, user.data.u, user.data.d, user );
    u = user.data.u;
    d = user.data.d;
    tau = user.data.tau;
    draw_nominal_trajectory = arg3_;
end

Nsamples = user.Nsamples;
Ninputs = user.Ninputs;
Nmodes = user.Nmodes;
u_min = user.u_min;
u_max = user.u_max;
% t0 = user.t0;
% tf = user.tf;
terminal_loc = user.terminal_loc;
px = user.idxs.px;
pz = user.idxs.pz;
theta = user.idxs.theta;
% dotx = user.idxs.dotx;
% dotz = user.idxs.dotz;
% dottheta = user.idxs.dottheta;
obstacle = user.obstacle;
L = user.starmac_params.length;
color = { 'b', [ 208 32 144 ]/255, 'r' };

% Trajectory figure
figure;
hold on;

% Parameters to draw the nominal trajectory
if ( draw_nominal_trajectory )
    zground = user.cost.instant_cost.zground;
    roll = user.cost.instant_cost.roll;
    t1 = user.cost.instant_cost.t1;
    t2 = user.cost.instant_cost.t2;
    t3 = user.cost.instant_cost.t3;
    Ntimes = 100;
    time = linspace(0,t1+t2+t3+1,Ntimes);
    
    for k = 1:Ntimes
        if( time(k) < t1 )
            xd(k) = time(k);
            zd(k) = zground;
            thetad(k) = roll;
        elseif( time(k) < t1 + t2 )
            xd(k) = t1 + cos( ( time(k) - t1 ) - pi/2 );
            zd(k) = zground + 2 + 2* sin( ( time(k) - t1 ) - pi/2 );
            thetad(k) = time(k) - t1 ;
        elseif( time(k) < t1 + t2 + t3 )
            xd(k) = time(k) - t2;
            zd(k) = zground;
            thetad(k) = roll;
        else
            xd(k) = t1 + t3;
            zd(k) = zground;
            thetad(k) = roll;
        end
    end
    plot( xd, zd, 'r--', 'LineWidth', 1.5);
end

% Draw the actual trajectory
figure;
hold on;
for i = 1:Nsamples
    helper = find( d( :, i ), 1, 'first' );
    plot( x( px, i:( i + 1 ) ), x( pz, i:( i + 1 ) ), '-', 'LineWidth', 4, 'Color', color{ helper } );
end

Del = 76;
for idx = 1:Del:( Nsamples + 1 )
    line( x( px, idx ) - L * cos( x( theta, idx ) ) * [ .5; -.5 ], ...
        x( pz, idx ) + L * sin( x( theta, idx ) ) * [ .5; -.5 ], ...
        'Color', [0 0 0], 'LineWidth', 2 );
    quiver( x( px, idx ), x( pz, idx ), ...
        0.25 * sin( x( theta, idx ) ), 0.25 * cos( x( theta, idx ) ), ...
        'LineWidth', 1.75, 'Color', [139 137 137]/255 );
end
idx = ( Nsamples + 1 );
line( x( px, idx ) - L * cos( x( theta, idx ) ) * [ .5; -.5 ], ...
    x( pz, idx ) + L * sin( x( theta, idx ) ) * [ .5; -.5 ], ...
    'Color', [0 0 0], 'LineWidth', 2 );
quiver( x( px, idx ), x( pz, idx ), ...
    0.25 * sin( x( theta, idx ) ), 0.25 * cos( x( theta, idx ) ), ...
    'LineWidth', 1.75, 'Color', [139 137 137]/255 );

% Draw the terminal location
plot( terminal_loc( px ), terminal_loc( pz ), 'ko', ...
    'MarkerFaceColor', [ 0 1 0 ] , 'MarkerSize', 7 );

% Draw the obstacles
for k = 1:length( obstacle )
    if( strcmpi( obstacle{k}.type, 'rectangle' ) )
        rectangle( 'Position', [ obstacle{k}.x-(obstacle{k}.width)/2 obstacle{k}.y-(obstacle{k}.height)/2 ...
            obstacle{k}.width obstacle{k}.height ], ...
            'EdgeColor', [139 115 85]/255,'FaceColor', [139 115 85]/255, 'LineStyle', '-', 'LineWidth', 2 );
    else
        warning( 'mip:starmac:plot_trajectories', 'Unknown obstacle type: %s.', obstacle{k}.type );
    end
end

% so that the plot looks reasonable...
xlabel( 'x_1', 'FontSize', 16 );
ylabel( 'x_2', 'FontSize', 16 );
set( get( gca, 'XLabel' ), 'Position', [ 3 -1.05 ] )
set( get( gca, 'YLabel' ), 'Position', [ -1.05 1 ] )
axis equal;
rectangle( 'Position', [-2 -1 17 .95], ...
    'EdgeColor', [139 115 85]/255, ...
    'FaceColor', [139 115 85]/255, ...
    'LineStyle', '-', ...
    'LineWidth', 2 );
axis( [ -1 7 -1 3 ] );
set( gca, 'XTick', [ -1 7 ], 'YTick', [ -1 3 ] );
hold off;


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
    plot( tau( 1:( end - 1 ) ), u( cnt, : ), '.-', 'MarkerSize', 5 );
    xlabel( 'Time' );
    ylabel( alabel );
    hold off;
end

% Draw the discrete inputs
figure;
for cnt = 1:Nmodes
    alabel = [ 'd_' num2str(cnt) ];
    subplot( Nmodes, 1, cnt );
    hold on;
    grid on;
    plot( tau( 1:( end - 1 ) ), d( cnt, : ), '.-', 'MarkerSize', 5 );
    xlabel( 'Time' );
    ylabel( alabel );
    hold off;
end
