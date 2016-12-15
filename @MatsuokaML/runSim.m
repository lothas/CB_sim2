function [out, sim, signal] = runSim(obj, sequence)
%SIM Runs a simulation of a CPG with Matsuoka neurons using the given
%parameters

% if nargin<3
%     beta = obj.Sim.Con.beta;
% end

%     cond1 = (obj.Sim.Con.tau_ratio + 1)*obj.Sim.Con.beta > max(obj.Sim.Con.W(:))
%     if ~cond1
%         obj.Sim.Con.beta = 2*max(obj.Sim.Con.W(:))/(obj.Sim.Con.tau_ratio + 1);
%     end

    % Instantiate copy of obj.Sim.
    sim.Mod.SetKeys = {};
    sim.Env.SetKeys = {};
    sim.Con = copy(obj.Sim.Con);
            
    % Setup tau_r, tau_a, c, W and feedback gains using genome
    % and beta
    sim = obj.Gen.Decode(sim, sequence);
%     sim.Con.beta = beta;
    
    sim.Con.s_in = 0;
    sim.Con = sim.Con.Adaptation();
    
    [x0, X, T] = runMatsuokaSim();
    [y, periods, signals, pos_work, neg_work, neuronActive, ...
        neuronOsc, perError1, perOK1, perError2, perOK2] = ...
        obj.processResults(X, T);
%     periods
%     [simFreq, amp] = obj.processResultsFFT(X, T, 0);
    
    % Plot results
    if obj.doPlot
        delta_y = max(X(:))-min(X(:));
        if ishandle(2)
            close(2)
        end
        figure(2)
        set(gcf,'units','normalized','position',[0.1 0.1 0.8 0.8])
        subplot(2,3,1)
        hold on
        for n = 1:obj.nNeurons
            plot(T,X(:,2*n-1)+(2.5-n)*delta_y);
        end
        subplot(2,3,4)
        hold on
        for n = 1:obj.nNeurons
            plot(T,X(:,2*n)+(2.5-n)*delta_y);
        end
        subplot(2,3,[2 3 5 6])
        hold on
        
        for n = 1:obj.nNeurons/2
            disp(['Output from pair ',int2str(n),':']);
            disp({'Period', 'Max y', 'Min y', 'Work +', 'Work -';
                periods(n), max(signals(:,n)), min(signals(:,n)), ...
                pos_work(n), neg_work(n)})

            plot(T, signals(:,n));

            title(['Period: ',num2str(max(periods),3),' s'])
        end
    end

    % Prepare output
    out.x0 = x0;
    out.periods = periods;
%     out.period_Rea = 1/simFreq;
    out.perError1 = perError1;
    out.perOK1 = perOK1;
    % ^ Returns 1 if the calculated period was verified to be correct
    out.perError2 = perError2;
    out.perOK2 = perOK2;
%     out.amp = amp;
    out.pos_work = pos_work;
    out.neg_work = neg_work;
    
    out.neuronActive = neuronActive;
    out.neuronOsc = neuronOsc;
    
    clear signal
    signal.T = T;
    signal.X = X;
    signal.y = y;
    signal.signal = signals;
    
    function Xt = derivative(~, X)
        Xt = sim.Con.Derivative(0, X);
%         Xt = zeros(size(X));
%         
%         y = max(X(1:2:end),0);
%         Wy = W*y;
%         
%         for i = 1:obj.nNeurons
%             Xt(2*i-1:2*i) = ...
%                 [(-X(2*i-1) - Wy(i,:) + c(i) - b*X(2*i))/Tr;
%                  (-X(2*i) + y(i))/Ta];
%         end
    end

    function [x0, X, T] = runMatsuokaSim()
        % Setup initial conditions
        x0 = zeros(2*obj.nNeurons,1);
%         x0(1:2:end) = (1-2*rand(obj.nNeurons,1)).*c/(1-b);
        x0(1:obj.nNeurons) = randn(obj.nNeurons,1)/6.*sim.Con.Amp0;

        tend = min(3*obj.tEnd, obj.tEnd/sim.Con.tau);
        tSpan = 0:obj.tStep:tend; % Time span
        
        % Set sim options
        options = odeset('AbsTol',obj.absTol,'RelTol',obj.relTol);

        % Run simulation
        [T,X] = ode45(@derivative, tSpan, x0, options);
    end
end
