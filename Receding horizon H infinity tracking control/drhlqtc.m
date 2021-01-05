function [x,y,u] = drhlqtc(x0,A,B,C,Q,R,Qf,N,yr,is)
    % DRHTC Discrete-time RHTC for Unconstrained Systems
    % convert reference vector to row vector
    [s1,s2] = size(yr); 
    if (s2 == 1)
        yr = yr';
    end
    % check future reference length for simulation
    if (length(yr) < (is + N))
        disp('The future reference is of too small length');
        return
    end
    % Riccati solution K(N-1) for RHTC
    [K_N_1 K_history_vec] = drde2(A,B,C,Q,R,Qf,N);
    % initialization of history variables
    xi = x0; % state of plant x(i)
    ui_history = []; % history of rhc u(i)
    xi_history = x0; % history of state of plant x(i)
    yi_history = C*x0; % history of output of plant y(i)
    for i=1:is
        % time-varying feed feedward gain g for RHTC
        [g_1 g_history] = dvde(A,B,C,Q,R,Qf,N,K_history_vec,yr,i);
        % receding horizon tracking controller u(i)
        ui = -inv(R+B'*K_N_1*B)*B'*(K_N_1*A*xi+g_1);
        % plant is controlled by rhtc u(i) at time [i,i+1]
        xi = A*xi + B*ui;
        yi = C*xi;
        % history u(i), x(i), y(i)
        ui_history = [ui_history ui];
        xi_history = [xi_history xi];
        yi_history = [yi_history yi];
    end
    x = xi_history; % state trajectory x
    y = yi_history; % output trajectory y
    u = ui_history; % RHTC history vector u
end