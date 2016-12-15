function [cond2] = check_cond2_forAll(T,tau,a,b,c,cond1, showFlag )
    
% Input:
% T, tau - time constants of the oscillator
% a - weigths matrix
% b - the inhibition coefficient
% showFlag - if '1' then display the calc

% Outputs:
% cond2 - linerazied system not stable (periodic solution) according to the
%           necessary condition and eigValues in RHP.


a_prime = a/(1+b);

v = de2bi((0:1:15)');

v(:,find(cond1)) = 1;

v1 = unique(v,'rows');

for i=1:(2^length(find(~cond1)))
    
    neurons2Check = v1(i,:);
    neurons2Check_index = find(neurons2Check);
    As = a(neurons2Check_index,neurons2Check_index); % make As from only the active neurons
    I = eye(length(As));

    LinSys = [ (-As-I)/tau , (-b/tau)*I; % TODO: this is not for any genral B vector!
                I/T        , -I/T          ;];

    eigLinSys = eig(LinSys);
    if (isempty(find(real(eigLinSys) > 0, 1 )))
        % check the continuous time stability condition
        cond2 = 0;
        break;
    else
        cond2 = 1;
    end
end



end


