function [currentRes] = plotCPG(results,num)
% this funtion takes the results and a CPG num and plot the CPG
% output

% Inputs: *) 'results' - structure contain the samples
%         *) 'num'- the relevant CPG num
%
% Outputs: *) 'currentRes' - the structure containing the relevant
%                           information

% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 30; % 15
MML.nNeurons = 4;

nSamples = length(results);

switch num
    case {'random','rand','Random'}
        examplNum = randsample(1:nSamples,1);
    otherwise
        examplNum = num;
end

currentRes = results(examplNum);
[~, ~, signal] = MML.runSim(currentRes.seq);

figure;
plot(signal.T,signal.X); grid minor;
xlabel('time[sec]');    ylabel('X_i');
title('X_i over time');
set(gca,'FontSize',15);

figure;
ax1 = axes('Position',[0 0 1 1],'Visible','off');
ax2 = axes('Position',[.3 .1 .6 .8]);
plot(ax2,signal.T,signal.signal(1,:)); hold on;
plot(ax2,signal.T,signal.signal(2,:)); grid minor;
xlabel({'time[sec]','seq=',num2str(currentRes.seq(1:4))});    
ylabel('Torque_i');
title('CPG output over time');
legend('y_2 - y_1','y_3-y_4');
set(ax2,'FontSize',15);

seq = currentRes.seq;

descr = {sprintf('CPG #%d',examplNum);
    'seq: ';
    sprintf('tau = %0.3f    b = %0.3f',seq(1,1),seq(1,2));
    sprintf('c_{i} = %0.3f, %0.3f, %0.3f, %0.3f',seq(1,3:6));
    sprintf(['W_{ij} = %0.3f, %0.3f, %0.3f, %0.3f \n \t \t'...
    '%0.3f, %0.3f, %0.3f, %0.3f \n \t \t'...
    '%0.3f, %0.3f, %0.3f, %0.3f'],seq(1,7:18));
    ' ';
    sprintf(['IC = %0.3f, %0.3f, %0.3f, %0.3f \n \t \t',...
    '%0.3f, %0.3f, %0.3f, %0.3f'],currentRes.x0);
    ' ';
    'identified period:';
    sprintf('AutoCorr = %0.3f,    %0.3f',currentRes.periods(:));
    sprintf('NL-LS = %0.3f,    %0.3f',currentRes.periods_LSQ(:))};
axes(ax1) % sets ax1 to current axes
text(.025,0.6,descr)
% sets ax2 back to current axis so I can zoom
axes(ax2)
end

