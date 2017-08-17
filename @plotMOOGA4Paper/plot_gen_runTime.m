function plot_gen_runTime(obj,gen_num)
% plot the runtime of each generation in MOGA
% 
% Inputs:
% *) 'gen_num' - the generation to focus

% generation numbers:
X_data = 1:gen_num;

figure; hold on;
ax = gca;

for i=1:4
   Y =  obj.data{1,i}.GA.totGenTime(1,X_data);
   plot(ax,X_data,Y); 
    
end

title('generation runTime');
xlabel('generation num');
ylabel('runTime [sec]');
legend(obj.legends);
grid minor;
hold off;

end

