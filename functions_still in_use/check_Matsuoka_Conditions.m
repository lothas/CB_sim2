function [ cond0, cond1, cond2, cond2UC, cond2_Sim, cond2UC_Sim] = ...
                            check_Matsuoka_Conditions(T,tau,a,b,c,activeNeurons_fromSim)
    
% Input:
% T, tau - time constants of the oscillator
% a - weigths matrix
% b - the inhibition coefficient
% 'activeNeurons_fromSim' - '1'= neuron is active in simulation
%                           '0'= neuron is inactive in simulation
% showFlag - if '1' then display the calc

% Outputs:
% cond0 - (((T-tau)^2)>=(4*tau*T*b))
% cond1 - necessary condition for active neuron according to one of
%           Matsuokas Lemmas.
% cond2 - linerazied system not stable (periodic solution) according to the
%           necessary condition and eigValues in RHP.
% cond2UC - linerazied system not stable (periodic solution) according to the
%           necessary condition and eigValues in the unit circle (discrete time conversion).
% cond2_Sim - linerazied system not stable (periodic solution) according to the
%           the neurons that were found in the Sim and eigValues in RHP.
% cond2UC_Sim - linerazied system not stable (periodic solution) according to the
%           the neurons that were found in the Sim and eigValues in the unit circle (discrete time conversion).

a_prime = a/(1+b);

cond0 = (((T-tau)^2)>=(4*b*tau*T));

cond1 = zeros(1,length(c));

%% testing with lemma 3
%         for i=1:length(c)  
%             cond1(i) = (c(i) > a_prime(i,:)*c(:)); 
%         end

%% testing with lemma 4
for i=1:length(c) % testing with lemma 4
    u_k_i = (c - a_prime(:,(1:end ~= i))*c(1:end ~= i));
    u_k_ii = max(0,u_k_i);
    v_j_i = (c - a_prime(:,(1:end ~= i))*u_k_ii(1:end ~= i));
    v_j_ii = max(0,v_j_i);
    cond1(i) = (c(i) > a_prime(i,:)*v_j_ii(:)); 
end

%% checking linearized matrix from the necessary condition:
ActiveneuronsIndex = find(cond1);
howManyActiveNeurons = length(ActiveneuronsIndex);

if ( howManyActiveNeurons > 0) % at least one neuron is active 
    As = a(ActiveneuronsIndex,ActiveneuronsIndex); % make As from only the active neurons
    I = eye(length(As));

    LinSys = [ (-As-I)/tau , (-b/tau)*I; % TODO: this is not for any genral B vector!
                I/T        , -I/T          ;];

    eigLinSys = eig(LinSys);
    eigLinSys_DT = eig(expm(LinSys)*0.01);

    if (isempty(find(real(eigLinSys) > 0, 1 )))
        % check the continuous time stability condition
        cond2 = 0;
    else
        cond2 = 1;
    end

    for j=1:length(eigLinSys_DT)
        if norm([real(eigLinSys_DT(j)),imag(eigLinSys_DT(j))]) > 1
            % check the sampled time stability condition
            cond2UC = 1;
            break; % no need to check again
        else
            cond2UC = 0;
        end
    end
else
    cond2 = 0; % canot check linearized system
    cond2UC = 0;
end

%% checking linearized matrix from the Sim results:
ActiveneuronsIndex_fromSim = find(activeNeurons_fromSim);
howManyActiveNeurons_fromSim = length(ActiveneuronsIndex_fromSim);

if ( howManyActiveNeurons_fromSim > 0) % at least one neuron is active 
    As1 = a(ActiveneuronsIndex_fromSim,ActiveneuronsIndex_fromSim); % make As from only the active neurons
    I1 = eye(length(As1));

    LinSys_Sim = [ (-As1-I1)/tau , (-b/tau)*I1; % TODO: this is not for any genral B vector!
                I1/T        , -I1/T          ;];

    eigLinSys_Sim = eig(LinSys_Sim);
    eigLinSys_DT_Sim = eig(expm(LinSys_Sim)*0.01);

    if (isempty(find(real(eigLinSys_Sim) > 0, 1 )))
        % check the continuous time stability condition
        cond2_Sim = 0;
    else
        cond2_Sim = 1;
    end

    for j=1:length(eigLinSys_DT_Sim)
        if norm([real(eigLinSys_DT_Sim(j)),imag(eigLinSys_DT_Sim(j))]) > 1
            % check the sampled time stability condition
            cond2UC_Sim = 1;
            break; % no need to check again
        else
            cond2UC_Sim = 0;
        end
    end
else
    cond2_Sim = 0; % canot check linearized system
    cond2UC_Sim = 0;
end


end

