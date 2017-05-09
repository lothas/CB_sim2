
clear all; close all; clc

syms x1 x2 x3 x4 u1 u2 u3 u4 real

x = [x1;x2;x3;x4];
u = [u1;u2;u3;u4];
y = x;%max([zeros(4,1),x],[],2);
%%
% fileName = 'MatsRandomRes_test.mat';
fileName = 'MatsRandomRes_02_02_2017.mat';
% fileName = 'MatsRandomRes_2-4_02_2017.mat';MatsRandomRes_02_02_2017
load(fileName)

periods = horzcat(results(:).periods);
periods = periods(1,:);

osc = ~isnan(periods);
good_ids = find(osc);
count = false(1,length(periods));
%%
for i=1:length(periods)
    
    tau = results(i).Tr;
    T = results(i).Ta;
    b = results(i).b;
    s = results(i).c;
    W_ij = results(i).W;

    sami = W_ij * y;
    f_i = (1/tau) * (-x + sami + s - b*u);
    f_i_prime = (1/T) * (-u + y);

    f = [f_i;f_i_prime];

    J = double(jacobian(f', [x',u']));
    
    % place '1' in every element when the Jacobian is negative
    J_bol = (J<0);
    
    % remove the elements on the diagonal, they are not interesting for
    % checking monotonic matrix
    J_bol_no_diagonal = J_bol-diag(diag(J_bol));
    
    % check if all of the elements that are not on the diagonal are
    % non-negative:
    if sum(sum(J_bol_no_diagonal)) == 0
        count(1,i) = true;
    end
end

disp(['the amount of ocsilattory CPGs is: ',num2str(sum(osc))]);
disp(['the amount of monotonic CPGs is: ',num2str(sum(count))]);
disp(['the amount of CPGs which fullfil both: ',num2str(sum(count & osc))]);
% jacobian([x*y*z, y^2, x + z], [x, y, z])