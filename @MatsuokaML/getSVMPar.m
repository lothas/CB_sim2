function [c, Worig, W, Tr, Ta, n_iter] = ...
    getSVMPar(obj, c, Worig, W, Tr, Ta)
%GETSVMPAR Use SVM to make sure that the selected parameters will
%result in the desired period
    in = [c;
         W(~logical(eye(size(W))));
         Tr];
     
    % Normalize X
    in = bsxfun(@rdivide, bsxfun(@minus, in, ...
        obj.norm_par(:,1)), obj.norm_par(:,2));

    count = 0;
    while n_iter<200
        if predict(obj.SVM, in')
            % Param. should converge to the desired period range
            break
        else
            % Param. won't converge to the desired period range
            % Try new parameters
            [c, Worig, W, Tr, Ta] = obj.getRandPar();

            in = [c;
                 W(~logical(eye(size(W))));
                 Tr];
             
            % Normalize X
            in = bsxfun(@rdivide, bsxfun(@minus, in, ...
                       obj.norm_par(:,1)), obj.norm_par(:,2));

            n_iter = n_iter + 1;
        end
    end
%     disp(['Number of iterations: ',int2str(n_iter)]);
end

