function result = save_result( result, batch_name, curr_date )

if( nargin >= 3 )
  filename = [ 'results' filesep batch_name '_' curr_date '.mat' ];
else
  filename = [ 'results' filesep batch_name '_' datestr( now, 30 ) '.mat' ];
end

ls = dir( '.' );
found_results_dir = 0;
for k = 1:length( ls )
  if( ls(k).isdir && strcmp( ls(k).name, 'results' ) )
    found_results_dir = 1;
    break;
  end
end
if( ~found_results_dir )
  mkdir( 'results' )
end

save( filename, 'result' );
fprintf( 'Successfully saved file: %s\n', filename );
