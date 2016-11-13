function dist = GetDistance( GA, Gen, type )
%GETDISTANCE Summary of this function goes here
%   Detailed explanation goes here
    if nargin<3
        type = 'seq';
    end
    
    if strcmp(type, 'seq')
        data = GA.Seqs(:,:,Gen);
    else
        data = GA.Fit(:,:,Gen);
    end
    
    Nn = GA.Population;
    dist = zeros(Nn,1);
    for i = 1:Nn
        % Find closest neighbor
        dist_mat = data - repmat(data(i,:),Nn,1);
        dist_vec = sqrt(diag(dist_mat*dist_mat'));
        
        % Remove self
        dist_vec(i) = [];
        dist(i) = min(dist_vec);
    end
end

