function [X_data,runTime] = plot_gen_runTime(obj,gen_num)
% plot the runtime of each generation in MOGA
% 
% Inputs:
% *) 'gen_num' - the generation to focus

% generation numbers:
X_data = 1:gen_num;
runTime = zeros(numel(obj.data_names),length(X_data));



for i=1:numel(obj.data_names)
       runTime(i,:) =  obj.data{1,i}.GA.totGenTime(1,X_data);
end

if nargout == 0
    figure; hold on;
    ax = gca;
    for i=1:numel(obj.data_names)
       plot(ax,X_data,runTime(i,:)); 
    end
    title('generation runTime');
    xlabel('generation num');
    ylabel('runTime [sec]');
    legend(obj.Legends);
    grid minor;
    hold off;
end


end

