function wantedParam_vec = get_param(obj,wantedParam_name)
% this function takes the results data structure and extract the wanted
% parameter (as a row vector)

% inputs: 
%   *) 'wantedParam_name' - the name of the wanted parameter as a string

% outputs:
%   *) 'wantedParam_vec' - a row vector with the wanted parameters values.

results = obj.sim_results;
ids = obj.ids;

switch wantedParam_name
    case 'tau'
        wantedParam_vec = vertcat(results(ids).tau);
    case 'b'
        wantedParam_vec = vertcat(results(ids).b);
    case 'a'
        wantedParam_vec = vertcat(results(ids).a);
    case 'c'
        wantedParam_vec = vertcat(results(ids).c);
    case {'freq','f'}
        periods_LSQ = vertcat(results(ids).periods_LSQ);
        % because the other harmonics are just multiplication
        wantedParam_vec = 1./periods_LSQ(:,1);
    case 'periods'
        periods_LSQ = vertcat(results(ids).periods_LSQ);
        wantedParam_vec = periods_LSQ(:,1);
    case {'A0','bias_coef'}
        wantedParam_vec = vertcat(results(ids).bias_coef);
    case {'A1','sine_coef_1'}
        sine_coef = vertcat(results(ids).sine_coef);
        wantedParam_vec = sine_coef(:,1);
    case {'A2','sine_coef_2'}
        sine_coef = vertcat(results(ids).sine_coef);
        wantedParam_vec = sine_coef(:,2);
    case {'A3','sine_coef_3'}
        sine_coef = vertcat(results(ids).sine_coef);
        wantedParam_vec = sine_coef(:,3);
    case {'B1','cos_coef_1'}
        cos_coef = vertcat(results(ids).cos_coef);
        wantedParam_vec = cos_coef(:,1);
    case {'B2','cos_coef_2'}
        cos_coef = vertcat(results(ids).cos_coef);
        wantedParam_vec = cos_coef(:,2);
    case {'B3','cos_coef_3'}
        cos_coef = vertcat(results(ids).cos_coef);
        wantedParam_vec = cos_coef(:,3);
    otherwise
        error('no parameter with that name...');
end

wantedParam_vec = wantedParam_vec';

%TODO:
% think how and for who to generate the desired parameters.


end

