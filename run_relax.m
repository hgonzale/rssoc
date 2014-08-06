function result = run_relax( config_file )

user = read_config( config_file );
beta = user.armijo.beta;

% set the Nsamples and  and construct the u0 and d0 that we initialize with
iter_Nsamples = user.Nsamples;
iter_tau = linspace( user.t0, user.tf, iter_Nsamples + 1 );
user = update_user( user );
len_cons = user.len_cons;

% check to make sure that each sample of the discrete mode sums to one
assert( max( abs( sum( user.d0, 1 ) - ones( 1, size( user.d0, 2 ) ) ) ) < 10 * eps );

s = state_encode( user, 0, user.u0, user.d0 );
[ ~ , iter_u, iter_d ] = state_decode( user, s );

for iter = 1:user.max_iter
  x = fwd_euler( user.x0, iter_tau, iter_u, iter_d, user, iter_Nsamples );
  if( len_cons > 0 )
    Psi = max_cons_fctn( user, iter_tau, iter_u, iter_d, x );
  else
    Psi = -Inf;
  end
  fprintf( 1, 'Iteration %d: Number of Samples = %d, Cost = %f, Psi = %f\n', ...
           iter-1, iter_Nsamples, obj_fctn( user, iter_tau, iter_u, iter_d, x ), Psi );

  user = update_user( user, iter_Nsamples, iter_tau, iter_u, iter_d );
  result.user( iter ) = user;

  fprintf( 1, '\tRunning optimality_fctn\n' );
  res_optfctn = optimality_fctn( user, [], [] );
  fprintf( 1, '\t\tTheta = %d\n', res_optfctn.value );
  
  result.optfctn( iter ) = res_optfctn;

  if ( res_optfctn.value > -user.numerical_tolerance )
    if( iter == 1 )
      % This is the code we should run whenever the optimality function is close to zero.
      [ iter_tau, iter_u, iter_d ] = resample_fctn( user, iter_tau, iter_u, iter_d );
      iter_Nsamples = size( iter_u, 2 );
      continue;
    else
      % We run this stopping condition instead, which is bad but practical.
      fprintf( 1, '\tSatisfied Numerical Tolerance\n' );
      fprintf( 1, '\tRunning pwm_fctn\n' );
      res_pwm = pwm_fctn( user, iter_tau, iter_u, iter_d, ...
        res_optfctn.value, res_optfctn.u_p, res_optfctn.d_p, res_insmode.k );

      result.res_pwm( iter ) = res_pwm;
      iter_u = res_pwm.u;
      iter_d = res_pwm.d;
      iter_tau = res_pwm.tau;
      iter_Nsamples = length( iter_tau ) - 1;
      break;
    end
  end

  fprintf( 1, '\tRunning armijo_fctn\n' );
  res_insmode = armijo_fctn( user, iter_tau, iter_u, iter_d, ...
                             res_optfctn.value, res_optfctn.u_p, res_optfctn.d_p );

  result.res_insmode( iter ) = res_insmode;

  if( mod( iter, user.pwm_iter ) == 0 )
    fprintf( 1, '\tRunning pwm_fctn\n' );
    res_pwm = pwm_fctn( user, iter_tau, iter_u, iter_d, ...
      res_optfctn.value, res_optfctn.u_p, res_optfctn.d_p, res_insmode.k );

    result.res_pwm( iter ) = res_pwm;

    if( res_pwm.k == Inf )
      [ iter_tau, iter_u, iter_d ] = resample_fctn( user, iter_tau, iter_u, iter_d );
      iter_Nsamples = size( iter_u, 2 );
    else
      iter_u = res_pwm.u;
      iter_d = res_pwm.d;
      iter_tau = res_pwm.tau;
      iter_Nsamples = length( iter_tau ) - 1;
    end
  else
    iter_u = iter_u + beta^res_insmode.k * ( res_optfctn.u_p - iter_u );
    iter_d = iter_d + beta^res_insmode.k * ( res_optfctn.d_p - iter_d );
    continue;
  end

end
user = update_user( user, iter_Nsamples, iter_tau, iter_u, iter_d );
result.user( iter + 1 ) = user;
if ( user.higher_order )
  [ x, cost ] = multistep_integ( user, iter_tau, iter_u, iter_d );
  fprintf( 1, 'Final cost = %d\n', obj_fctn( user, iter_tau, iter_u, iter_d, [ x; cost ] ));
else
  x = fwd_euler( user.x0, iter_tau, iter_u, iter_d, user );
  fprintf( 1, 'Final cost = %d\n', obj_fctn( user, iter_tau, iter_u, iter_d, x ));
end

save_result( result, user.batch_name );
end
