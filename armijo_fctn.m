function res = armijo_fctn( user, tau, u, d, theta, unew, dnew )

len_cons = user.len_cons;
alpha = user.armijo.alpha;
beta = user.armijo.beta;
kmax = user.armijo.kmax;
higher_order = user.higher_order;

tic();
if( len_cons > 0 )
    if ( higher_order )
        [ x, cost ] = multistep_integ( user, tau, u, d );
        Psi = max_cons_fctn( user, tau, u, d, [ x; cost ] );
    else
        Psi = max_cons_fctn( user, tau, u, d );
    end
    if( Psi <= 0 )
        if ( higher_order )
            J = obj_fctn( user, tau, u, d, [ x; cost ] );
        else
            J = obj_fctn( user, tau, u, d );
        end
        
        success = 0;
        for k = 0:kmax
            if ( higher_order )
                [ x, cost ] = multistep_integ( user, tau, u + beta^k * ( unew - u ),  d + beta^k * ( dnew - d ) );
                Jnew = obj_fctn( user, tau, u + beta^k * ( unew - u ), d + beta^k * ( dnew - d ), [ x; cost ] );
                Psinew =  max_cons_fctn( user, tau, u + beta^k * ( unew - u ), d + beta^k * ( dnew - d ), [ x; cost ] );
            else
                Jnew = obj_fctn( user, tau, u + beta^k * ( unew - u ), d + beta^k * ( dnew - d ) );
                Psinew = max_cons_fctn( user, tau, u + beta^k * ( unew - u ), d + beta^k * ( dnew - d ) );
            end
            
            if( Jnew - J <= alpha * beta^k * theta && ...
                    Psinew <= alpha * beta^k * theta )
                success = 1;
                break;
            end
        end
        if( ~success )
            warning( 'relax:armijo_fctn', 'Maximum number of interations reached' );
        end
    else % psi is positive
        success = 0;
        for k = 0:kmax
            if ( higher_order )
                [ x, cost ] = multistep_integ( user, tau, u + beta^k * ( unew - u ),  d + beta^k * ( dnew - d ) );
                Psinew =  max_cons_fctn( user, tau, u + beta^k * ( unew - u ), d + beta^k * ( dnew - d ), [ x; cost ] );
            else
                Psinew = max_cons_fctn( user, tau, u + beta^k * ( unew - u ), d + beta^k * ( dnew - d ) );
            end
            
            if( Psinew - Psi <= alpha * beta^k * theta )
                success = 1;
                break;
            end
        end
        if( ~success )
            warning( 'relax:armijo_fctn', 'Maximum number of interations reached' );
        end
    end
else % unconstrained problem.
    if ( higher_order )
        [ x, cost ] = multistep_integ( user, tau, u, d );
        J = obj_fctn( user, tau, u, d, [ x; cost ] );
    else
        J = obj_fctn( user, tau, u, d );
    end
    
    success = 0;
    for k = 0:kmax
        if ( higher_order )
            [ x, cost ] = multistep_integ( user, tau, u + beta^k * ( unew - u ),  d + beta^k * ( dnew - d ) );
            Jnew = obj_fctn( user, tau, u + beta^k * ( unew - u ), d + beta^k * ( dnew - d ), [ x; cost ] );
        else
            Jnew = obj_fctn( user, tau, u + beta^k * ( unew - u ), d + beta^k * ( dnew - d ) );
        end
        
        if( Jnew - J <= alpha * beta^k * theta )
            success = 1;
            break;
        end
    end
    if( ~success )
        warning( 'relax:armijo_fctn', 'Maximum number of interations reached' );
    end
end
tend = toc();

res.k = k;
res.time = tend;

