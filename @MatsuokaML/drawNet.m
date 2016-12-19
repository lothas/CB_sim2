function drawNet( ~, W )
%DRAWNET Summary of this function goes here
%   Detailed explanation goes here
% Implemented only for 4 neurons

nNeurons = 4;
radius = 0.3;
centers = [-1,1;1,1;1,-1;-1,-1];
conColor = [0 0 0];
conrad = 0.2*radius;

maxW = max(abs(W));
if maxW > 10
    % Normalize weights between -10 and 10
    W = W/maxW*10;
end

figure

% Draw weights
Wmat = zeros(nNeurons, nNeurons);
k = 1;
for i = 1:nNeurons
    for j = 1:nNeurons
        if i == j
            continue
        end
        Wmat(i,j) = W(k);
        drawCon(j,i,1,W(k));
        k = k+1;
    end
end

% Draw neurons
for i = 1:size(centers,1)
    pos = [centers(i,:) - [radius, radius], ...
        2*radius, 2*radius];
    rectangle('Position',pos, 'Curvature',[1,1], ...
        'FaceColor',[0.8 0.9 0.7], 'LineWidth', 2)
end
% Add numbers
text(centers(:,1), centers(:,2), {'1','2','3','4'}, ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 24, 'FontWeight', 'bold')

axis_lim = 1.5;
axis([-axis_lim, axis_lim, -axis_lim, axis_lim])
axis equal
axis off
    function drawCon(i,j,dir,weight)
        if weight == 0
            return
        end
        
        p1 = centers(i,:);
        p2 = centers(j,:);
        v = p2-p1;
        vp = [-v(2), v(1)];
        d = norm(v);
        if dir == 1
            p11 = p1+0.75*radius*v/d+0.3*radius*vp/d;
            p21 = p2-1.25*radius*v/d+0.3*radius*vp/d;
        else
            p11 = p1+1.25*radius*v/d-0.3*radius*vp/d;
            p21 = p2-0.75*radius*v/d-0.3*radius*vp/d;
        end
            
        line([p11(1) p21(1)],[p11(2) p21(2)],...
            'LineWidth', abs(weight), 'Color', conColor)
        
        % Add inhibitory/exitatory connection
        if dir == 1
            conpos = [p21 - [conrad, conrad], 2*[conrad, conrad]];
        else
            conpos = [p11 - [conrad, conrad], 2*[conrad, conrad]];
        end
        if weight>0
            % Draw filled circle
            rectangle('Position', conpos, 'Curvature', [1,1], ...
                'EdgeColor', conColor, 'FaceColor', conColor, ...
                'LineWidth', 2)
        else
            % Draw empty circle
            rectangle('Position', conpos, 'Curvature', [1,1], ...
                'EdgeColor', conColor, 'FaceColor', [1 1 1], ...
                'LineWidth', 2)
        end
    end
end


