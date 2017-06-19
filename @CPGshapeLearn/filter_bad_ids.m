function [good_ids] = filter_bad_ids(obj,results)

% NOTE: this function is for symmetric Matsuoka CPG with 2 neurons!

% this function takes only CPG's with not NaN period
% also, this function filtered CPG's with abnormally big Fourier coeff

%% % % % taking only CPGs with oscilations: 
% 1) take according to Jonathan's autocorrelation method (for verifaction
%       purposes)
periods = vertcat(results(:).periods);
ids_Jon = (~(isnan(periods)))';
% 2) taking the period from the nonlinear least squar (LSQ) fit method
ids_LSQ = false(1,length(periods));
for i=1:length(periods)
    temp = results(i).periods_LSQ;
    ids_LSQ(1,i) = ~(isnan(temp(1,1)));
end
ids_good_period = ids_Jon & ids_LSQ;

%% % % % taking only CPGs with reasonable coeffients:
ids_good_bias = true(1,length(periods));
ids_good_Sine = true(1,length(periods));
ids_good_cos = true(1,length(periods));

coef_tol = 10; % the thrushold for home much big is too big

for i=1:length(periods)
    % check bias coef:
    temp = results(i).bias_coef;
    if abs(temp) > coef_tol
        ids_good_bias(1,i) = false;
    end
    
    % check Sine coef:
    temp = results(i).sine_coef;
    for j=1:size(temp,2)
        if temp(1,j) > coef_tol
            ids_good_Sine(1,i) = false;
            break;
        end
    end
    
    % check cosSine coef:
    temp = results(i).cos_coef;
    for j=1:size(temp,2)
        if temp(1,j) > coef_tol
            ids_good_cos(1,i) = false;
            break;
        end
    end
end
ids_good_coef = ids_good_bias & ids_good_Sine & ids_good_cos;

%% % % % The "good" ids:
good_ids = ids_good_period & ids_good_coef;


end

