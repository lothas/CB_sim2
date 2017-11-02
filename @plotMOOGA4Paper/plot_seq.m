function [ax] = plot_seq(obj,seq,Title)
%PLOT_SEQ plots the CPG output for a given sequence

[out, ~, signal] = obj.MML.runSim(seq');
plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');
xlabel('time[sec]');    ylabel('CPG output');
title( {Title,...
    ['    periods: ',...
    num2str(out.periods(1)),'    ',...
    num2str(out.periods(2))]});

ax = gca;


end