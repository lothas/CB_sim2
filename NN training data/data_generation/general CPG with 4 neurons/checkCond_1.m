function cond1 = checkCond_1(seq,seqOrder)
% this function check the condition:
%       cond1(i) = (c(i) > a_prime(i,:)*c(:));

% checking which neuron are expected to be active.

%the 'seqOrder' is usually: {'tau','b','c_1','c_2','c_3','c_4',...
%     'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
%     'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};
%

seq = seq(1,1:18);

cond1 = false(1,4);

tau = seq(strcmp('tau',seqOrder));
T = 5*tau;
% WARNING: make sure that T is indeed 5*tau

b = seq(strcmp('b',seqOrder));
if b<0
   disp('warning b<0 ');
end

% only for 4 neurons CPG:
c_1 = seq(strcmp('c_1',seqOrder));
c_2 = seq(strcmp('c_2',seqOrder));
c_3 = seq(strcmp('c_3',seqOrder));
c_4 = seq(strcmp('c_4',seqOrder));

c = [c_1,c_2,c_3,c_4];
if find(c<0,1)
    disp('warning one of c<0 ');
end
   
w_12 = seq(strcmp('w_{12}',seqOrder));
w_13 = seq(strcmp('w_{13}',seqOrder));
w_14 = seq(strcmp('w_{14}',seqOrder));
w_21 = seq(strcmp('w_{21}',seqOrder));
w_23 = seq(strcmp('w_{23}',seqOrder));
w_24 = seq(strcmp('w_{24}',seqOrder));
w_31 = seq(strcmp('w_{31}',seqOrder));
w_32 = seq(strcmp('w_{32}',seqOrder));
w_34 = seq(strcmp('w_{34}',seqOrder));
w_41 = seq(strcmp('w_{41}',seqOrder));
w_42 = seq(strcmp('w_{42}',seqOrder));
w_43 = seq(strcmp('w_{43}',seqOrder));

W = [0   ,w_12,w_13 ,w_14;
    w_21 ,0   ,w_23 ,w_24;
    w_31 ,w_32 ,0   ,w_34;
    w_41 ,w_42,w_43 ,0  ];

% % % normalize:
W = (diag(1./c)*(diag(c)*W)')';

a_prime = W/(1+b);

% % % testing with lemma 3:
for i=1:length(c)  
    cond1(i) = (c(i) > a_prime(i,:)*c(:)); 
end

% % % % testing with lemma 4:
% for i=1:length(c) % testing with lemma 4
%     u_k_i = (c - a_prime(:,(1:end ~= i))*c(1:end ~= i));
%     u_k_ii = max(0,u_k_i);
%     v_j_i = (c - a_prime(:,(1:end ~= i))*u_k_ii(1:end ~= i));
%     v_j_ii = max(0,v_j_i);
%     cond1(i) = (c(i) > a_prime(i,:)*v_j_ii(:)); 
% end

end