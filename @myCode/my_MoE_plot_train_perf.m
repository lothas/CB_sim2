function my_MoE_plot_train_perf(obj)
% plot only graphs which are relevant to training sessions:

legendNames = cell(1,obj.expertCount);
for j=1:obj.expertCount
    legendNames{1,j} = ['#',num2str(j),' expert'];
end

iterNum = 1:obj.numOfIteretions;

if obj.my_MoE_out.competetiveFlag==1 || obj.my_MoE_out.competetiveFlag==2
    figure;
    subplot(2,2,1) % cluster size over iteration number
    plot(iterNum,obj.my_MoE_out.expertsTrainData.expert_i_GroupSize,'-o'); hold on;
    title('cluster size over #iteration');
    xlabel('#iteretion');   ylabel('group size [#points]');
    legend(legendNames);

    subplot(2,2,2) % expert perf over #interation
    plot(iterNum,obj.my_MoE_out.expertsTrainData.Experts_perf_mat,'-o'); hold on;
    title('expert MSE over #interation');
    xlabel('#iteretion');   ylabel('performance [MSE]');
    legend(legendNames);
    
    subplot(2,2,3) % expert perf over #interation
    plot(iterNum,double(obj.my_MoE_out.expertsTrainData.emptyGroupIndecator),'-o'); hold on;
    title('indication on empty clusters: "1" means empty');
    xlabel('#iteretion');   ylabel('"1"=empty   "0"-not empty');

end

figure; 
subplot(2,1,1) % total MSE error over #iteration
plot(iterNum,obj.my_MoE_out.Moe_perf_over_iter,'b-o');
xlabel('#iteration'); ylabel('MoE MSE error');
title('total MSE error over #iteration');

subplot(2,1,2) % gateNet perf over #interation
plot(iterNum,obj.my_MoE_out.gateTraniData.gateNN_perf_vec,'-o'); hold on;
title('gateNet perf (MSE) over #interation');
xlabel('#iteretion');   ylabel('performance [crossentropy]');

% plot the weights of the experts:
for k=1:floor(obj.expertCount/9)
    lastExpertIND = (k-1)*9;
    figure; hold on;
    for j=1:9
        currentExptIND = lastExpertIND+j;
        subplot(3,3,j);
        title_temp = (['expert #',num2str(currentExptIND)]);
        obj.NN_weights_matrix_visualize(obj.my_MoE_out.expertsNN{1,currentExptIND},title_temp);
    end
end

figure;
plotregression(obj.targ_train,obj.my_MoE_out.out_from_train,'train',...
    obj.targ_valid,obj.my_MoE_out.out_from_valid,'validation',...
    obj.targ_test,obj.my_MoE_out.out_from_test,'test');

end

