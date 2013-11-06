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
px = user.idxs.px;
py = user.idxs.py;
color = { 'b', [ 208 32 144 ]/255, 'r' };
scrsz = get( 0, 'ScreenSize' );

% Trajectory figure
figure('Position',[20 scrsz(4)-700 560 420]);
hold on;
Del = 3;
for i = 1:Del:( Nsamples - Del + 1 )
    helper = find( d( :, i ), 1, 'first' );
    plot( tau( i:( i + Del ) ), x( px, i:( i + Del ) ), 'o-', 'linewidth', 2, 'Color', color{ helper } ); 
    plot( tau( i:( i + Del ) ), x( py, i:( i + Del ) ), 's-', 'linewidth', 2, 'Color', color{ helper } );
    ylabel('Tank Height', 'FontSize', 16 );  
    xlabel('Time [s]', 'FontSize', 16 );
end
axis( [ 0 10 2 3.51 ] );
hold off;

% figure('Position',[20 scrsz(4)-700 560 420]);
% hold on;
% for i = 1:Nsamples
%     helper = find( d( :, i ), 1, 'first' );
%     plot( tau( i:( i + 1 ) ), x( py, i:( i + 1 ) ), '-', 'linewidth', 2, 'Color', color{ helper } );
%     ylabel('x_2', 'FontSize', 16 );  
%     xlabel('Time [s]', 'FontSize', 16 );
% end
% axis( [ 0 10 2 3.51 ] );
% hold off;



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
  alabel = [ 'd_' num2str(cnt) ];
  subplot( Nmodes, 1, cnt );
  hold on;
  grid on;
%   plot( tau( 1:( end - 1 ) ), d( cnt, : ), '.-', 'MarkerSize', 5 );
  plot( tau_unique( 1:(end - 1 ) ), d( cnt, I( 1:(end-1) ) ), '.-', 'MarkerSize', 5, 'LineWidth', 1.5 );
  xlabel( 'Time', 'FontSize', 21 );
  ylabel( alabel, 'FontSize', 21 );
  hold off;
end

