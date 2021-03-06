function [x, y, u, w] = drhtc_hinf (x0, A, B, Bw, C, Q, R, Rw, Qf, gamma_2, N, yr, is)
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
    [K_N_1 K_history_vec] = drde2_hinf(A,B,Bw,C,Q,R,Rw,Qf,gamma_2,N);
    % initialization of history variables
    xi = x0; % state of plant x(i)
    ui_history = []; % history of rhc u(i)
    xi_history = x0; % history of state of plant x(i)
    yi_history = C*x0; % history of output of plant y(i)
    wi_history = [];
    system_order = size(A,1);
    for i=1:is
        % time-varying feed feedward gain g for RHTC
        [g_1 g_history] = dvde_hinf(A ,B ,Bw ,C ,Q ,R ,Rw ,Qf ,gamma_2, N, K_history_vec, yr, i);
        % receding horizon tracking controller u(i)
        Lambda = eye(system_order) + K_N_1 * ( B*inv(R)*B' - 1/gamma_2 * Bw*inv(Rw)*Bw');
        ui = -inv(R) * B' * inv( Lambda )*( K_N_1 * A * xi + g_1 ) ;
        wi = 1/gamma_2 * inv(Rw) * Bw' * inv ( Lambda ) * ( K_N_1 * A * xi + g_1 ) ;
        % plant is controlled by rhtc u(i) at time [i,i+1]
        xi = A*xi + B*ui + Bw*wi;
        yi = C*xi;
        ui_history = [ui_history ui];
        xi_history = [xi_history xi];
        yi_history = [yi_history yi];
        wi_history = [wi_history wi];
    end
    x = xi_history; % state trajectory x
    y = yi_history; % output trajectory y
    u = ui_history; % RHTC history vector u
    w = wi_history; % RHTC history vector u
end