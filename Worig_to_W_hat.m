function [ W ] = Worig_to_W_hat( Worig,c )
% this function takes W and convert it to W_hat

   for i=1:size(Worig,2)
       a = Worig(:,i);
       c_i = c(:,i);
       Aorig = [0    ,a(1) ,a(2) ,a(3);
            a(4) ,0    ,a(5) ,a(6);
            a(7) ,a(8) ,0    ,a(9);
            a(10),a(11),a(12),0   ];

       A = (diag(1./c_i)*(diag(c_i)*Aorig)')';
       
       W(:,i) = [A(1,2);A(1,3);A(1,4);...
              A(2,1);A(2,3);A(2,4);...
              A(3,1);A(3,2);A(3,4);...
              A(4,1);A(4,2);A(4,3)];
   end
end

