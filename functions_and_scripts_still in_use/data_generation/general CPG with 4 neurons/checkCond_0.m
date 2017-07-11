function cond0 = checkCond_0(seq,seqOrder)
% this function check the condition:
%       cond0 = (((T-tau)^2)>=(4*b*tau*T));
%the 'seqOrder' is usually: {'tau','b','c_1','c_2','c_3','c_4',...
%     'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
%     'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};
%

tau = seq(strcmp('tau',seqOrder));
T = 5*tau;
% WARNING: make sure that T is indeed 5*tau

b = seq(strcmp('b',seqOrder));

cond0 = (((T-tau)^2)>=(4*b*tau*T));

end
