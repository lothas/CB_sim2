function DrawKey(GA, Keys, Segments)
%DRAWKEY Summary of this function goes here
%   Detailed explanation goes here

colors = {[236, 93, 93]; [80, 217, 66]; [95, 180, 246];
          [105, 219, 212]; [236, 177, 93]; [224, 236, 93];
          [124, 117, 210]; [210, 117, 172]};
fontsize = 12;
      
if nargin<3
    Keys = GA.Gen.Keys;
    Segments = GA.Gen.Segments;
end

seg_width = 1/sum(Segments);

figure('units', 'normalized', 'Position', [0.2, 0.3, 0.44, 0.55])
axis equal
axis off
hold on
i = 1;
x = 0;
x0 = x;
for s = 1:size(Segments,2)
    if Segments(s) == 0
        continue
    end
    
    for l = 1:Segments(s)
        rectangle('Position', [x, 0, seg_width, seg_width], ...
            'FaceColor', colors{s}/255.0, 'EdgeColor', [0 0 0])
        
        % Add gene number label
        text(x+seg_width/2, seg_width/2, int2str(i), ...
            'HorizontalAlignment', 'center', 'FontSize', fontsize); 
        
        i = i + 1;
        x = x + seg_width;
    end
    
    % Add label
    text((x+x0)/2, 1.7*seg_width, Keys{1,s}, 'FontSize', fontsize, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold'); 
%     text((x+x0)/2, 1.5*seg_width, Keys{1,s}, 'Rotation', 90, ...
%         'FontWeight', 'bold', 'FontSize', fontsize); 
    
    x0 = x;
end

end

