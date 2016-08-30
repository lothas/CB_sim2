function [out, signal] = sim(obj, b, c, W, Tr, Ta)
%SIM Runs a simulation of a CPG with Matsuoka neurons using the given
%parameters

    [x0, X, T] = runSim();
    [y, periods, signals, pos_work, neg_work] = obj.processResults(X, T);

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
    out.pos_work = pos_work;
    out.neg_work = neg_work;
    
    clear signal
    signal.T = T;
    signal.X = X;
    signal.y = y;
    signal.signal = signals;
    
    function Xt = derivative(~,X)
        Xt = zeros(size(X));
        
        y = max(X(1:2:end),0);
        Wy = W*y;
        
        for i = 1:obj.nNeurons
            Xt(2*i-1:2*i) = ...
                [(-X(2*i-1) - Wy(i,:) + c(i) - b*X(2*i))/Tr;
                 (-X(2*i) + y(i))/Ta];
        end
    end

    function [x0, X, T] = runSim()
        % Setup initial conditions
        x0 = zeros(2*obj.nNeurons,1);
%         x0(1:2:end) = (1-2*rand(obj.nNeurons,1)).*c/(1-b);
        x0(1:2:end) = randn(obj.nNeurons,1)/6.*c;

        tSpan = 0:obj.tStep:obj.tEnd; % Time span
        
        % Set sim options
        options = odeset('AbsTol',obj.absTol,'RelTol',obj.relTol);

        % Run simulation
        [T,X] = ode45(@derivative, tSpan, x0, options);
    end
end
