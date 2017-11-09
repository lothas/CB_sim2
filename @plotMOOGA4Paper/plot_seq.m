function [ax] = plot_seq(obj,seq,Title)
%PLOT_SEQ plots the CPG output for a given sequence

[out, ~, signal] = obj.MML.runSim(seq');

titleAdd = sprintf('    periods: ');
hold on;
for i=1:size(signal.signal,1)
    plot(signal.T,signal.signal(i,:));
    titleAdd = [titleAdd,sprintf('    %.3f',out.periods(i))];
end
xlabel('time[sec]');    ylabel('CPG output');
title( {Title,titleAdd});

ax = gca;


end