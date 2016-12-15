function condMat = checkCondMat(nSims,results,periods,neuronActive)

% Inputs:
% nSims - number of simulations.
% results- contain the genome sequence of the CPG parameters 
% 'activeNeurons_fromSim' - '1'= neuron is active in simulation
%                           '0'= neuron is inactive in simulation
%
%
% Outputs:
% condMat - logic variable of all of the conditions and simulation data.
%

% % DEFINE AS 'logic': much better:) 
A_positive = false(nSims,1);    % all of the weigths are positive
has_period = false(nSims,1);    % sequnce has a period
cond0 = false(nSims,1);         %(((T-tau)^2)>=(4*tau*T))
cond1 = false(nSims,4);         % activ neurons
cond2 = false(nSims,1);         % LinSys stability,  eig in RHP
cond2UC = false(nSims,1);       % LinSys stability, eig out of unit circle
cond2_Sim = false(nSims,1);     % LinSys stability (active neurons from the Sim),  eig in RHP
cond2UC_Sim = false(nSims,1);   % LinSys stability (active neurons from the Sim), eig out of unit circle

check_period =  isnan(periods(1,:));

for i=1:nSims
   % check if this sequence has a period
   if check_period(1,i)
       has_period(i,1) = 0;
   else
       has_period(i,1) = 1;
   end
   
   % define the CPG parameters
   tau = results(i).seq(1);
   T = tau*5;
   b = results(i).seq(2);
   if b<0
       disp(['warning b<0 at i = ',num2str(i)]);
   end
   c = (results(i).seq(3:6))';
   if find(c<0,1)
        disp(['warning one of c<0 at i = ',num2str(i)]);
   end
   
   a = results(i).seq(7:18);
   Aorig = [0    ,a(1) ,a(2) ,a(3);
        a(4) ,0    ,a(5) ,a(6);
        a(7) ,a(8) ,0    ,a(9);
        a(10),a(11),a(12),0   ];
   A = (diag(1./c)*(diag(c)*Aorig)')';
   % TODO: use getGenes for better performance!!
    
   checkPositiveA = find((A < 0),1,'first');

   % check if all of the weigths are positive
   if (~isempty(checkPositiveA))
       A_positive(i,1) = 0;
   else
       A_positive(i,1) = 1;
   end
   
   % check the conditions:
   [ cond0(i,1), cond1(i,:), ~, cond2UC(i,1), cond2_Sim(i,1), cond2UC_Sim(i,1)] = ...
                            check_Matsuoka_Conditions(T,tau,A,b,c,neuronActive(i,:));
   
   % checking linearization for all possible stationary solutions: 
   [cond2(i,1)] = check_cond2_forAll(T,tau,A,b,c,cond1(i,:), 0 );
end

condMat = [ A_positive ,cond0, cond1, cond2, cond2UC, cond2_Sim, cond2UC_Sim, has_period];

end

