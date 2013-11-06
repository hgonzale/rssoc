function res = pwm_fctn( user, tau, u, d, theta, unew, dnew, mu )

Nsamples = user.Nsamples;
len_cons = user.len_cons;
alpha = user.armijo.alpha;
beta = user.armijo.beta;
alphabar = user.pwm.alphabar;
betabar = user.pwm.betabar;
eta = user.pwm.eta;
omega = user.pwm.omega;
higher_order = user.higher_order;

uprime = u + beta^mu * ( unew - u );
dprime = d + beta^mu * ( dnew - d );

N0 = log2( Nsamples );

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
        for k = 1:( N0 + eta )
            [ taurho, urho, drho ] = rho_fctn( user, k, tau, uprime, dprime );
            
            if ( higher_order )
                [ xrho, costrho ] = multistep_integ( user, taurho, urho, drho );
                Jnew = obj_fctn( user, taurho, urho, drho, [ xrho; costrho ] );
                Psinew = max_cons_fctn( user, taurho, urho, drho, [ xrho; costrho ] );
            else
                Jnew = obj_fctn( user, taurho, urho, drho );
                Psinew = max_cons_fctn( user, taurho, urho, drho );
            end
            
%             LHS = Jnew - J
%             RHS = ( alpha * beta^mu - alphabar * betabar^k ) * theta
%             LHS_post = alphabar * betabar^k
%             RHS_post = alpha * beta^mu
            if( Jnew - J <= ( alpha * beta^mu - alphabar * betabar^k ) * theta && ...
                    Psinew <= 0 && ...
                    alphabar * betabar^k <= ( 1 - omega ) * alpha * beta^mu ) % HG: we can do better than this logic
                %             if( Jnew - J <= ( alpha * beta^mu - alphabar * betabar^k ) * theta && ...
                %                     Psinew <= 0 ) % HG: we can do better than this logic
                success = 1;
                break;
            end
        end
        %         if ( ~success )
        %             k = N0 + 1;
        %         end
        if( ~success )
            k = Inf;
        end
    else % psi is positive
        success = 0;
        for k = 1:( N0 + eta )
            [ taurho, urho, drho ] = rho_fctn( user, k, tau, uprime, dprime );
            if ( higher_order )
                [ xrho, costrho ] = multistep_integ( user, taurho, urho, drho );
                Psinew = max_cons_fctn( user, taurho, urho, drho, [ xrho; costrho ] );
            else
                Psinew = max_cons_fctn( user, taurho, urho, drho );
            end
            
            if( Psinew - Psi <= ( alpha * beta^mu - alphabar * betabar^k ) * theta && ...
                    alphabar * betabar^k <= ( 1 - omega ) * alpha * beta^mu ) % HG: we can do better than this logic
                success = 1;
                break;
            end
        end
        if( ~success )
            k = Inf;
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
    for k = 1:( N0 + eta )
        [ taurho, urho, drho ] = rho_fctn( user, k, tau, uprime, dprime );
        if ( higher_order )
            [ xrho, costrho ] = multistep_integ( user, taurho, urho, drho );
            Jnew = obj_fctn( user, taurho, urho, drho, [ xrho; costrho ] );
        else
            Jnew = obj_fctn( user, taurho, urho, drho );
        end
        
        %     Jnew - J
        %      ( alpha * beta^mu - alphabar * betabar^k ) * theta
        %      alphabar * betabar^k
        %      ( 1 - omega ) * alpha * beta^mu
        if( Jnew - J <= ( alpha * beta^mu - alphabar * betabar^k ) * theta && ...
                alphabar * betabar^k <= ( 1 - omega ) * alpha * beta^mu ) % HG: we can do better than this logic
            success = 1;
            break;
        end
    end
    if( ~success )
        k = Inf;
    end
end
tend = toc();

res.u = urho;
res.d = drho;
res.tau = taurho;
res.k = k;
res.time = tend;
