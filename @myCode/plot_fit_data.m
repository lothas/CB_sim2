function plot_fit_data(obj,freq_net_est,problemType)
% plot the fitting data. for example, in the symmetric 2 neurons case we
% can plot on a 3D graph the change of 'a','tau' and the frequency on the
% z-axis and than we can compare it for example to the Matsuoka estimation formula.

% NOTE: pay attention to changes in the seq order, what row of 'sampl' is
%       what CPG parameter!

targ_test = obj.targ_test;
sampl_test = obj.sampl_test;

freq_from_code = targ_test;

switch problemType
    case '2N_symm'
        freq_Matsuoka_est = zeros(1,size(targ_test,2));
        % seqOrder should be: 1)tau 2)b 3)a 4)s
        tau = sampl_test(1,:);
        b = sampl_test(2,:);
        a = sampl_test(3,:);
        T = 5.*tau; % note: make sure that it is indeed '5' in the sim
        for j=1:size(targ_test,2);
            freq_Matsuoka_est(1,j) = obj.MatsuokaEstimation(tau(1,j),T(1,j),b(1,j),a(1,j));
        end
        
        if mean(tau)==tau % than plot 2D graph
            figure;
            scatter(a,freq_from_code,'o'); hold on
            scatter(a,freq_Matsuoka_est,'x');
            scatter(a,freq_net_est,'d');
            xlabel('a');    ylabel('freq [Hz]');   grid on;
            legend('freq from code','Matsuoka est','MoE or NN estimation');
            title('frequency over "a"'); hold off;
        else
            figure;
            scatter3(a,tau,freq_from_code,'o'); hold on
            scatter3(a,tau,freq_Matsuoka_est,'x');
            scatter3(a,tau,freq_net_est,'d');
            xlabel('a');    ylabel('tau');    zlabel('freq [Hz]'); grid on;
            legend('freq from code','Matsuoka est','MoE or NN estimation');
            title('frequency over "a" and "tau"'); hold off;
        end
        
    case '2N_general'
            % inputs seqOrder should be: 1)tau 2)b 3)w12 4)w21
            w12 = sampl_test(3,:);
            w21 = sampl_test(4,:);
            
            figure;
            scatter3(w12,w21,freq_from_code,'o'); hold on
            scatter3(w12,w21,freq_net_est,'x');
            xlabel('w_{21}');    ylabel('w_{21}');    zlabel('freq [Hz]'); grid on;
            legend('freq from code','MoE or NN estimation');
            title('frequency over the CPG weights'); hold off;
        
        %TODO: add more plotting options!
end


end

