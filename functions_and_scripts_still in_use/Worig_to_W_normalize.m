function W_norm = Worig_to_W_normalize(Worig,c )
% this function takes W and convert it to W_hat

% INPUTS: *) 'Worig' - a 12 X N matrix with all of the weights.
%         *) 'c' - a 4 X N matrix with all of the tonic inputs.
%         
% OUTPUTS: *) 'W_norm' - a 12 X N matrix contain all of the normalize weights

N = size(c,2);
W_norm = zeros(size(Worig));

for i=1:N
   a = Worig(:,i);
   c_i = c(:,i);
   Aorig = [0    ,a(1) ,a(2) ,a(3);
        a(4) ,0    ,a(5) ,a(6);
        a(7) ,a(8) ,0    ,a(9);
        a(10),a(11),a(12),0   ];

   A = (diag(1./c_i)*(diag(c_i)*Aorig)')';

   W_norm(:,i) = [A(1,2);A(1,3);A(1,4);...
                    A(2,1);A(2,3);A(2,4);...
                    A(3,1);A(3,2);A(3,4);...
                    A(4,1);A(4,2);A(4,3)];
end


% % % % ANOTHER METHOD:
% % we can get the normalized coupling weights straight from the sim results.
% N = length(results);
% W_norm = zeros(12,N);
% for i=1:N
%     A = results(i).W;
%     W_norm(:,i) = [A(1,2);A(1,3);A(1,4);...
%                     A(2,1);A(2,3);A(2,4);...
%                     A(3,1);A(3,2);A(3,4);...
%                     A(4,1);A(4,2);A(4,3)];
% end

end

