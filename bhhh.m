function [ MLE ] = bhhh( thetaHat, dataR, n, spec)

d_nLogLike = [];
for k = 1:numel(dataR)
    [ probChosen, d_probChosen ]  = ProbitProb( thetaHat(dataR{k}.pick), dataR{k});
    
    d = zeros(n.theta, dataR{k}.n.con);
    d(dataR{k}.pick, :) = d_probChosen;
    
    d_nLogLike  = [d_nLogLike, bsxfun( @rdivide, d, probChosen')];
end

if spec.constraint == 2
    MLE.cov = inv(d_nLogLike*d_nLogLike');
else
    deltaid         = sort(paramsid.delta(identifiable.delta));
    betaid          = 1:n.theta;
    betaid(deltaid) = [];
    
    [~, ~, ~, d_s]  = ShareConstraints( thetaHat, dataS, identifiable, n, shareHat);
    
    d_s_delta       = d_s(deltaid,:);
    d_s_beta        = d_s(betaid,:);
    d_f_delta       = d_nLogLike(deltaid,:);
    d_f_beta        = d_nLogLike(betaid,:);
    
    d_delta_beta    = -d_s_beta/d_s_delta;
    d_f             = d_f_beta + d_delta_beta*d_f_delta;
    
    d_theta_beta    = zeros(n.theta - n.delta, n.theta);
    d_theta_beta(:, deltaid) = d_delta_beta;
    d_theta_beta(:, betaid)  = eye(n.theta - n.delta);
    
    MLE.cov = d_theta_beta'/(d_f*d_f')*d_theta_beta;
end

% Standard errors
MLE.se      = sqrt( diag( MLE.cov ) );

% t-statistics
MLE.tStat   = thetaHat ./ MLE.se;

% p > |t|
MLE.p       = 2 * ( 1 - normcdf( abs( MLE.tStat ) ) );

% 95% Confidence Interval
MLE.ci_lb   = thetaHat - 1.96 * MLE.se;
MLE.ci_ub   = thetaHat + 1.96 * MLE.se;

% MLE value
MLE.value   = -MLE.value;

% AIC
MLE.AIC     = 2 * ( n.theta - MLE.value );

% BIC
MLE.BIC     = -2 * MLE.value + n.theta * log( n.con );
end

