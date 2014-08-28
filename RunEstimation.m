function thetaHat = RunEstimation(dataR, dataS, theta_0, n, spec, opt, shareHat, identifiable, paramsid)
%% Run Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bounds on all parameters
ub  = ones( n.theta, 1 ) * 30;
lb  = ones( n.theta, 1 ) * -30;

% Define objective function
obj = @(x) LogLike( x, dataR );
constraints = @(x) ShareConstraints(x, dataS, identifiable, n, shareHat);

HPattern = zeros(n.theta);
index = [ (1:n.conGroup+1)'
    paramsid.beta_1(identifiable.beta_1)
    paramsid.beta_2(identifiable.beta_2)
    paramsid.S(identifiable.S)];
HPattern(index,:) = 1;
HPattern(:,index) = 1;

for m = 1:n.market
    index_m = paramsid.delta(identifiable.delta(:,m),m);
    HPattern(index_m, index_m) = 1;
end

JPattern = HPattern;
JPattern(index, :) = [];

HPattern = sparse(HPattern);
JPattern = sparse(JPattern);

% index = 1+n.conGroup+(n.maxChoice-1)*(n.conChar + n.prodChar);

% Run Optimization
tic

optOption   = optimset( 'MaxIter', opt.maxIter, ...
    'MaxFunEvals', opt.maxFunEvals, ...
    'Display', opt.display, ...
    'TolFun', opt.tolFun, ...
    'TolCon', opt.tolCon, ...
    'TolX', opt.tolX, ...
    'Algorithm', opt.algorithm, ...
    'SubproblemAlgorithm', 'cg', ...
    'GradObj', opt.gradObj, ...
    'GradConstr', opt.gradConstr );

optOption2   = optimset( 'MaxIter', opt.maxIter, ...
    'MaxFunEvals', opt.maxFunEvals, ...
    'Display', opt.display, ...
    'TolFun', opt.tolFun, ...
    'TolCon', opt.tolCon, ...
    'TolX', opt.tolX, ...
    'Algorithm', opt.algorithm, ...
    'JacobPattern', JPattern, ...
    'HessPattern', HPattern, ...
    'SubproblemAlgorithm', 'cg', ...
    'GradObj', opt.gradObj, ...
    'GradConstr', opt.gradConstr );

if spec.solver == 1    
    if spec.constraint == 1
        [ thetaHat, MLE.value ] = ...
            knitromatlab( obj, theta_0, [], [], [], [], lb, ub, ...
            constraints, [], optOption2 );
        
    elseif spec.constraint == 2
        [ thetaHat, MLE.value ] = ...
            knitromatlab( obj, theta_0, [], [], [], [], lb, ub, ...
            [], [], optOption );
    end
    
elseif spec.solver == 2    
    optOption   = optimset( 'MaxIter', opt.maxIter, ...
        'MaxFunEvals', opt.maxFunEvals, ...
        'Display', opt.display, ...
        'TolFun', opt.tolFun, ...
        'TolX', opt.tolX );
    
    [ thetaHat, MLE.value ] = fminsearch( obj, theta_0, optOption );
    
elseif spec.solver == 3
    
    optOption   = optimset( 'MaxIter', opt.maxIter, ...
        'MaxFunEvals', opt.maxFunEvals, ...
        'Display', opt.display, ...
        'TolFun', opt.tolFun, ...
        'TolCon', opt.tolCon, ...
        'TolX', opt.tolX, ...
        'Algorithm', opt.algorithm, ...
        'GradObj', opt.gradObj, ...
        'GradConstr', opt.gradConstr, ...
        'DerivativeCheck', 'off' );
    
    if spec.constraint == 1
        
        [ thetaHat, MLE.value ] = ...
            fmincon( obj, theta_0, [], [], [], [], lb, ub, ...
            constraints, optOption);
        
    elseif spec.constraint == 2
        
        [ thetaHat, MLE.value ] = ...
            fmincon( obj, theta_0, [], [], [], [], lb, ub, ...
            [], optOption );
    end
    
elseif spec.solver == 4
    funcs.objective = obj;
    funcs.gradient = @(x) LogLikeGrad(x, dataR);
    
    options.lb = lb;
    options.ub = ub;
    
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.max_iter = opt.maxIter;
    options.ipopt.print_level = 5;
    options.ipopt.print_frequency_iter = 10;
    options.ipopt.jacobian_regularization_value = 1e-3;
    options.ipopt.limited_memory_max_history = 30;
    
    thetaHat = ipopt(theta_0, funcs, options);
    MLE.value = obj(thetaHat);
end

% % Estimation time in seconds
% MLE.time    = toc;
% 
% %% Compute MLE Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Define negative log likelihood function
% nLoglike    = @( x, data, cens, freq ) obj( x );
% 
% % Compute MLE covariance matrix by finite difference
% fprintf('\n Computing MLE Covariance Matrix...');
% MLE.covOpt  = statset( 'GradObj', opt.gradObj );
% %MLE.cov     = mlecov( thetaHat, data, 'nloglf', nLoglike, ...
% %                      'options', MLE.covOpt );
% 
% % BHHH
% 
% d_nLogLike = [];
% for k = 1:numel(dataR)
%     [ probChosen, d_probChosen ]  = ProbitProb( thetaHat(dataR{k}.pick), dataR{k});
%     
%     d = zeros(n.theta, dataR{k}.n.con);
%     d(dataR{k}.pick, :) = d_probChosen;
%     
%     d_nLogLike  = [d_nLogLike, bsxfun( @rdivide, d, probChosen')];
% end
% 
% if spec.constraint == 2
%     MLE.cov = inv(d_nLogLike*d_nLogLike');
% else
%     deltaid         = sort(paramsid.delta(identifiable.delta));
%     betaid          = 1:n.theta;
%     betaid(deltaid) = [];
%     
%     [~, ~, ~, d_s]  = ShareConstraints( thetaHat, dataS, identifiable, n, shareHat);
%     
%     d_s_delta       = d_s(deltaid,:);
%     d_s_beta        = d_s(betaid,:);
%     d_f_delta       = d_nLogLike(deltaid,:);
%     d_f_beta        = d_nLogLike(betaid,:);
%     
%     d_delta_beta    = -d_s_beta/d_s_delta;
%     d_f             = d_f_beta + d_delta_beta*d_f_delta;
%     
%     d_theta_beta    = zeros(n.theta - n.delta, n.theta);
%     d_theta_beta(:, deltaid) = d_delta_beta;
%     d_theta_beta(:, betaid)  = eye(n.theta - n.delta);
%     
%     MLE.cov = d_theta_beta'/(d_f*d_f')*d_theta_beta;
% end
% 
% % Standard errors
% MLE.se      = sqrt( diag( MLE.cov ) );
% 
% % t-statistics
% MLE.tStat   = thetaHat ./ MLE.se;
% 
% % p > |t|
% MLE.p       = 2 * ( 1 - normcdf( abs( MLE.tStat ) ) );
% 
% % 95% Confidence Interval
% MLE.ci_lb   = thetaHat - 1.96 * MLE.se;
% MLE.ci_ub   = thetaHat + 1.96 * MLE.se;
% 
% % MLE value
% MLE.value   = -MLE.value;
% 
% % AIC
% MLE.AIC     = 2 * ( n.theta - MLE.value );
% 
% % BIC
% MLE.BIC     = -2 * MLE.value + n.theta * log( n.con );
