function view_change_in_param_as_video(data,...
    paramNames,gen_num,...
    fitnessOrder,seqOrder_extend,titleAdd)
% takes two parameters and present them on a 2D graph that each frame
% represents a generation:

filename = 'testAnimated1.gif';

h = figure;

ax = cell(1,4);
for i=1:4
    ax{1,i} = subplot(2,2,i);
    set(ax{1,i},'NextPlot','replace');
end

for g=1:gen_num
    
    for i=1:4

        Y = cell(1,2);

        for j=1:2
            switch paramNames{1,j}
                case {'VelFit','NrgEffFit',...
                     'VelRangeFit #1','VelRangeFit #2','VelRangeFit #3',...
                     'VelRangeFit #4','VelRangeFit #5','VelRangeFit #6',...
                     'VelRangeFit #7','VelRangeFit #8','EigenFit'}
                    ID = strcmp(paramNames{1,j},fitnessOrder);
                    Y{1,j} = squeeze(data{1,i}.GA.Fit(:,ID,g));
                    Axis = [0,0.5,0,0.5];


                case {'tau','b','c_1','c_2','c_3','c_4',...
                      'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
                      'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}',...
                      'ks_\tau','ks_c1','ks_c2','ks_c3','ks_c4'}
                     ID = strcmp(paramNames{1,j},seqOrder_extend);
                     Ytemp = squeeze(data{1,i}.GA.Seqs(:,ID,g));
                     Y{1,j} = (Ytemp-min(min(Ytemp)))/...
                         (max(max(Ytemp))-min(min(Ytemp))); % norm 'Y'
                     Axis = [0,1,0,1];
                otherwise
                    error('invalid Param...');
            end
        end

        scatter(ax{1,i},Y{1,1},Y{1,2});
        xlabel(ax{1,i},paramNames{1,1});
        ylabel(ax{1,i},paramNames{1,2});
        title(ax{1,i},{['gen #',num2str(g)],['case: ',titleAdd{1,i}]});
        axis(ax{1,i},Axis);
        
        drawnow
        
        
        
    end

    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if g == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
      
    pause(0.1);
end


end

