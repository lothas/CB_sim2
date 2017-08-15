function plot_gen_runTime(data,gen_num,names)
% plot the runtime of each generation in MOGA
% 
% Inputs:
% *) 'data' - cell array contain MOGa results
% *) 'gen_num' - the generation to focus
% *) 'names' - casaes names for legend

% generation numbers:
X_data = 1:gen_num;

figure; hold on;
ax = gca;

for i=1:4
   Y =  data{1,i}.GA.totGenTime(1,X_data);
   plot(ax,X_data,Y); 
    
end

title('generation runTime');
xlabel('generation num');
ylabel('runTime [sec]');
legend(names);
grid minor;
hold off;

end

