function res = optimality_fctn( user, u_p, d_p )

Nmodes = user.Nmodes;
Nsamples = user.Nsamples;
Ninputs = user.Ninputs;
len_cineq = user.len_cineq;

if( nargin == 1 )
    u_p = randn( Ninputs, Nsamples );
    d_p = rand( Nmodes, Nsamples );
else
    if( isempty( u_p ) )
        u_p = zeros( Ninputs, Nsamples );
    end
    if( isempty( d_p ) )
        d_p = zeros( Nmodes, Nsamples );
    end
end

% precomputing the STM, dfdu, and f
tic();
user = state_variation_precompute( user );
tend = toc();

if( strcmp( user.optfctn.solver, 'tomlab' ) == 1 )
    Name = [ user.batch_name '_optimality_fctn' ];
    
    switch( user.optfctn.tomlab_solver )
        case{ 'filterSQP','snopt' }
            user.qp_cons = 1; % Turning this flag off to make it clear that we are solving it with quadratic constraints
            
            s_0 = state_encode( user, 0, u_p, d_p );
            s_L = state_encode( user, -Inf, user.u_min, zeros( Nmodes, 1 ) );
            s_U = state_encode( user, Inf,  user.u_max, ones( Nmodes, 1 ) );
            c_L = -Inf * ones( len_cineq, 1 );
            c_U = zeros( len_cineq, 1 );
            fLowBnd = [];
            [ A, b_L, b_U ] = dinput_encode_lin_cons( user ); % ensures that the sum of the discrete inputs at each sample sum to one
            ConsPattern = [];
            Prob = conAssign( @optfctn_obj_f, @optfctn_obj_g, [], [], ...
                s_L, s_U, Name, s_0, ...
                [], fLowBnd, ...
                A, b_L, b_U, ...
                @optfctn_cons_f, @optfctn_cons_g, [], ConsPattern, c_L, c_U );
        case{ 'cplex', 'minlpBB', 'qpSolve', 'lssol', 'knitro', 'qpopt', 'sqopt' }
            if ( user.qp_cons )
                s_0 = state_encode( user, 0, u_p, d_p );
                s_L = state_encode( user, -Inf, user.u_min, zeros( Nmodes, 1 ) );
                s_U = state_encode( user, Inf,  user.u_max, ones( Nmodes, 1 ) );
                [ A, b_L, b_U ] = dinput_encode_lin_cons( user ); % ensures that the sum of the discrete inputs at each sample sum to one
                [ F, c, qc ] = encode_QP_cons( user ); % generates the cost and matrices for the QPQC formulation
                Prob = miqqAssign( F, c, A, b_L, b_U, s_L, s_U, s_0, qc,...
                    [], [], [], [], Name);
            else
                % the optimization is taking place here over (up - u) and
                % (dp - d) therefore all the upper and lower bounds are
                % different here.
                s_0 = state_encode( user, 0, u_p, d_p );
                dummyD = ones( Nmodes, Nsamples );
                s_L = state_encode( user, -Inf, repmat( user.u_min, 1, Nsamples ) - user.data.u, 0 * dummyD - user.data.d );
                s_U = state_encode( user, 0, repmat( user.u_max, 1, Nsamples ) - user.data.u, dummyD - user.data.d );
                [ F, c, A, b_L, b_U ] = encode_LP_cons( user ); % generates the cost and matrics QP formulation
                Prob = qpAssign( F, c, A, b_L, b_U, s_L, s_U, s_0, Name );
            end
    end
    
    
    switch( user.optfctn.tomlab_solver )
        case 'knitro'
            Prob.KNITRO.options.BLASOPTION = 1;
            Prob.KNITRO.options.BLASOPTIONLIB = 'libblas64_mp.so';
    end
    Prob.user = user;
    Prob.Solver.Alg = 1; % Use Newton's method with BFGS Hessian updates
    %     if( checkDerivs( Prob, s_0, 0 ) ~= 0 )
    %         error( 'The derivatives are WRONG!.' );
    %     end
    
    tic();
    tlres = tomRun( user.optfctn.tomlab_solver, Prob );
    tend = toc() + tend;
    
    if( tlres.ExitFlag ~= 0 )
        warning( 'relax:optimality_fctn', strcat( 'Solver %s (%s) could not find an optimal solution.\n', ...
            '\tExitFlag = %d.\n', ...
            '\tExitText = "%s".' ), ...
            user.optfctn.solver, user.optfctn.tomlab_solver, tlres.ExitFlag, tlres.ExitText );
    end
    
    if ( user.qp_cons )
        [ res.zeta, res.u_p, res.d_p ] = state_decode( user, tlres.x_k );
    else
        [ res.zeta, up_u, dp_d ] = state_decode( user, tlres.x_k );
        res.u_p = up_u + user.data.u;
        res.d_p = dp_d + user.data.d;
    end
    
    res.value = tlres.f_k;
    res.s_0 = s_0; % keep in mind if qp_cons == 0 then this is dp_d
    res.s_end = tlres.x_k'; % keep in mind if qp_cons == 0 then this is dp_d
    res.inform = tlres.ExitFlag;
    res.time = tend;
    res.user = user;
    res.tlres = tlres;
    if ( res.value > user.numerical_tolerance )
        warning('what the fuck');
    end
elseif( strcmp( user.optfctn.solver, 'cplex' ) == 1 )
    % the optimization is taking place here over (up - u) and
    % (dp - d) therefore all the upper and lower bounds are
    % different here.
    s_0 = state_encode( user, 0, u_p, d_p );
    dummyD = ones( Nmodes, Nsamples );
    s_L = state_encode( user, -Inf, repmat( user.u_min, 1, Nsamples ) - user.data.u, 0 * dummyD - user.data.d );
    s_U = state_encode( user, 0, repmat( user.u_max, 1, Nsamples ) - user.data.u, dummyD - user.data.d );
    [ F, c, A, b_L, b_U ] = encode_LP_cons( user ); % generates the cost and matrics QP formulation
    
    % cplex cannot handle upper and lower bound constraints simulateneously
    % so we have reformat things to make it more amenable
    Aineq = [ -A; A ];
    bineq = [ -b_L; b_U ];
    
    %     call cplex to solve:
    %     min     0.5*x'*F*x+f*x or c*x
    %     st.     Aineq*x     <= bineq
    %             s_L <= x <= s_U
    [ res.s_end, res.value, res.inform, tlres ] = cplexqp( F, c, Aineq, bineq, [], [], s_L, s_U );
    
    % time to decode the output
    [ res.zeta, up_u, dp_d ] = state_decode( user, res.s_end );
    res.u_p = up_u + user.data.u;
    res.d_p = dp_d + user.data.d;
    res.s_0 = s_0; % keep in mind this is dp_d
    res.time = tend + tlres.time;
    res.user = user;
    res.tlres = tlres;
    if ( res.value > user.numerical_tolerance )
        warning('what the fuck');
    end
elseif( strcmp( user.optfctn.solver, 'quadprog' ) == 1 )
    % the optimization is taking place here over (up - u) and
    % (dp - d) therefore all the upper and lower bounds are
    % different here.
    s_0 = state_encode( user, 0, u_p, d_p );
    dummyD = ones( Nmodes, Nsamples );
    s_L = state_encode( user, -Inf, repmat( user.u_min, 1, Nsamples ) - user.data.u, 0 * dummyD - user.data.d );
    s_U = state_encode( user, 0, repmat( user.u_max, 1, Nsamples ) - user.data.u, dummyD - user.data.d );
    [ F, c, A, b_L, b_U ] = encode_LP_cons( user ); % generates the cost and matrics QP formulation
    
    % quadprog cannot handle upper and lower bound constraints simulateneously
    % so we have reformat things to make it more amenable
    Aineq = [ -A; A ];
    bineq = [ -b_L; b_U ];
    
    %     call quadprog to solve:
    %     min     0.5*x'*F*x+f*x or c*x
    %     st.     Aineq*x     <= bineq
    %             s_L <= x <= s_U
    tic;
    [ res.s_end, res.value, res.inform, ~ ] = cplexqp( F, c, Aineq, bineq, [], [], s_L, s_U );
    res.time = tend + toc();
    
    % time to decode the output
    [ res.zeta, up_u, dp_d ] = state_decode( user, res.s_end );
    res.u_p = up_u + user.data.u;
    res.d_p = dp_d + user.data.d;
    res.s_0 = s_0; % keep in mind this is dp_d
    res.user = user;
    res.tlres = [];
    if ( res.value > user.numerical_tolerance )
        warning('what the fuck');
    end
elseif( strcmp( user.optfctn.solver, 'gurobi' ) == 1 )
    % the optimization is taking place here over (up - u) and
    % (dp - d) therefore all the upper and lower bounds are
    % different here.
    s_0 = state_encode( user, 0, u_p, d_p );
    dummyD = ones( Nmodes, Nsamples );
    s_L = state_encode( user, -Inf, repmat( user.u_min, 1, Nsamples ) - user.data.u, 0 * dummyD - user.data.d );
    s_U = state_encode( user, 0, repmat( user.u_max, 1, Nsamples ) - user.data.u, dummyD - user.data.d );
    [ F, c, A, b_L, b_U ] = encode_LP_cons( user ); % generates the cost and matrics QP formulation
    
    % gurobi cannot handle upper and lower bound constraints simulateneously
    % so we have reformat things to make it more amenable
    Aineq = [ -A; A ];
    bineq = [ -b_L; b_U ];
    
    % setup problem for gurobi to solve
    Prob.A = sparse( Aineq );
    Prob.Q = sparse( 0.5 * F );
    Prob.obj = c;
    Prob.rhs = bineq;
    Prob.sense = '<';
    Prob.lb = s_L;
    Prob.ub = s_U;
    
    %     call gurobi to solve:
    %     min     x'*Q*x+obj*x
    %     st.     A*x     <= rhs
    %             lb <= x <= ub
    stupidparams.outputflag = 0;
    tlres = gurobi( Prob, stupidparams );
    
    % time to decode the output
    res.s_end = tlres.x;
    res.value = tlres.objval;
    res.inform = tlres.status;
    
    
    [ res.zeta, up_u, dp_d ] = state_decode( user, res.s_end );
    res.u_p = up_u + user.data.u;
    res.d_p = dp_d + user.data.d;
    res.s_0 = s_0; % keep in mind this is dp_d
    res.time = tend + tlres.runtime;
    res.user = user;
    res.tlres = tlres;
    if ( res.value > user.numerical_tolerance )
        warning('what the fuck');
    end
elseif( strcmp( user.optfctn.solver, 'mosek' ) == 1 )
    % the optimization is taking place here over (up - u) and
    % (dp - d) therefore all the upper and lower bounds are
    % different here.
    s_0 = state_encode( user, 0, u_p, d_p );
    dummyD = ones( Nmodes, Nsamples );
    s_L = state_encode( user, -Inf, repmat( user.u_min, 1, Nsamples ) - user.data.u, 0 * dummyD - user.data.d );
    s_U = state_encode( user, 0, repmat( user.u_max, 1, Nsamples ) - user.data.u, dummyD - user.data.d );
    [ F, c, A, b_L, b_U ] = encode_LP_cons( user ); % generates the cost and matrics QP formulation
    
    %     call mosek to solve:
    %     min     1/2*x'*Q*x+c*x
    %     st.     b_L <=  A*x <= b_U
    %             lb <= x <= ub
    tic;
    tlres = mskqpopt( F, c, A, b_L, b_U, s_L, s_U, [], 'minimize echo(0)' );
    res.time = tend + toc();
    
    % time to decode the output
    res.s_end = tlres.sol.itr.xx;
    res.value = tlres.sol.itr.pobjval;
    res.inform = tlres.sol.itr.solsta;
    
    
    [ res.zeta, up_u, dp_d ] = state_decode( user, res.s_end );
    res.u_p = up_u + user.data.u;
    res.d_p = dp_d + user.data.d;
    res.s_0 = s_0; % keep in mind this is dp_d
    res.user = user;
    res.tlres = tlres;
    if ( res.value > user.numerical_tolerance )
        warning('what the fuck');
    end
else
    error( 'Unsupported OCP solver: %s', user.optfctn.solver );
end


%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

function out = optfctn_obj_f( s, Prob )

[ zeta, ~, ~ ] = state_decode( Prob.user, s );

out = zeta; % RV: This is a dummy variable, do not confuse with \zeta(\xi,\eta) in the paper.


%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

function out = optfctn_obj_g( ~, Prob )

len_s = Prob.user.len_s;
zidx = Prob.user.idxs.zeta;

out = zeros( 1, len_s );
out( 1, zidx ) = 1; % RV: This is a dummy variable, do not confuse with \zeta(\xi,\eta) in the paper.

