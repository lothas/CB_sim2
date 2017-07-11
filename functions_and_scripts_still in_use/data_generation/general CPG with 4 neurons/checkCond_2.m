function cond2 = checkCond_2(seq,seqOrder,cond1)
% this function check the condition:
%       checking the eigenvalues of the linearize matrix


%the 'seqOrder' is usually: {'tau','b','c_1','c_2','c_3','c_4',...
%     'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
%     'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};
%

seq = seq(1,1:18);

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
% cancel if not normalizing the weights.
W = (diag(1./c)*(diag(c)*W)')';

% % % % % check the second condition
ActiveneuronsIndex = find(cond1); % find which neurons are expected to be active
howManyActiveNeurons = length(ActiveneuronsIndex);

if ( howManyActiveNeurons > 0) % at least one neuron is active 
    As = W(ActiveneuronsIndex,ActiveneuronsIndex); % make As from only the active neurons
    I = eye(length(As));

    LinSys = [ (-As-I)/tau , (-b/tau)*I; % TODO: this is not for any genral B vector!
                I/T        , -I/T          ;];

    eigLinSys = eig(LinSys);

    if (isempty(find(real(eigLinSys) > 0, 1 )))
        % check the continuous time stability condition
        cond2 = false;
    else
        cond2 = true;
    end

else
    cond2 = false; % canot check linearized system
end


end