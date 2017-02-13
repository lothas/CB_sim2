function plot_fit_data(obj,method,problemType)
% plot the fitting data. for example, in the symmetric 2 neurons case we
% can plot on a 3D graph the change of 'a','tau' and the frequency on the
% z-axis and than we can compare it for example to the Matsuoka estimation formula.

% NOTE: pay attention to changes in the seq order, what row of 'sampl' is
%       what CPG parameter!

targ_test = obj.targ_test;
sampl_test = obj.sampl_test;

freq_from_code = targ_test;

switch method
    case {'our_MoE','my_MoE'}
        freq_net_est = obj.my_MoE_out.out_from_test;
        gateOut = obj.my_MoE_out.gateNet(sampl_test);
        [g_max,g_max_ind] = max(gateOut,[],1);
    case {'paper_MoE'}
        freq_net_est = obj.paper_MoE_out.out_from_test;
        gateOut = obj.paper_MoE_out.gateOut_from_test;
        [g_max,g_max_ind] = max(gateOut,[],1);
    case {'NN'}
        freq_net_est = obj.NN.out_from_test;
    otherwise
        error('invalid method');
end

if obj.expertCount == 2
    colors = [1,0,0; 0,0,1]; 
elseif obj.expertCount > 2
    colors = rand(obj.expertCount,3); % if more than 2 experts
end
            
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
            
            % plotting a 2 experts output (freq over 'a') with color coding
            if exist('g_max','var')==1 % only plot this for MoE methods
                figure;
                scatter(a,freq_from_code,'go'); hold on
                scatter(a,freq_Matsuoka_est,'bo');
                for j=1:obj.expertCount
                    for i=1:size(a,2)
                        if g_max_ind(1,i) == j
                            if g_max(1,i) > 0.5 % only samples with more than 50% get fillet dot
                                scatter(a(1,i),freq_net_est(1,i),'ko','MarkerFaceColor',colors(j,:));
                            else
                                scatter(a(1,i),freq_net_est(1,i),'ko');
                            end
                        end
                    end
                end
                xlabel('a');    ylabel('freq [Hz]'); grid on;
                legendNames = {'real freq','Matsuoka est','our MoE'}; %'Matsuoka estimation'
                legend(legendNames);
                hold off
            end
        else
            figure;
            scatter3(a,tau,freq_from_code,'o'); hold on
            scatter3(a,tau,freq_Matsuoka_est,'x');
            scatter3(a,tau,freq_net_est,'d');
            xlabel('a');    ylabel('tau');    zlabel('freq [Hz]'); grid on;
            legend('freq from code','Matsuoka est','MoE or NN estimation');
            title('frequency over "a" and "tau"'); hold off;
        end
        
        [MSE_matsuoka,~] = obj.NN_perf_calc(freq_from_code,freq_Matsuoka_est,0,0,'test' )
        disp(['Matsuoka estimation perf = ',num2str(MSE_matsuoka)]);
        
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
        
                        % plotting a 2 experts output (freq over 'a') with color coding
            if exist('g_max','var')==1 % only plot this for MoE methods
                figure;
                scatter3(w12,w21,freq_from_code,'go'); hold on
                for j=1:obj.expertCount
                    for i=1:size(w12,2)
                        if g_max_ind(1,i) == j
                            if g_max(1,i) > 0.5 % only samples with more than 50% get fillet dot
                                scatter3(w12(1,i),w21(1,i),freq_net_est(1,i),'ko','MarkerFaceColor',colors(j,:));
                            else
                                scatter3(w12(1,i),w21(1,i),freq_net_est(1,i),'ko');
                            end
                        end
                    end
                end
                xlabel('w12');  ylabel('w21');    zlabel('freq [Hz]'); grid on;
                legendNames = {'real freq','our MoE'}; %'Matsuoka estimation'
                legend(legendNames);
                hold off
            end
        %TODO: add more plotting options!
end


end

