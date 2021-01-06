function [g_1,g_history ] = dvde(A,B,C,Q,R,Qf,N,K_history_vec,yr,i)
    %DVDE Discrete-time Vector Differential Equation Solver
    % system dimension
    n = size(A,1);
    [s1,s2] = size(K_history_vec);
    % boundary condition
    g_N = -C'*Qf*yr(:,i+N); 
    g_history = g_N; 
    g_j = g_N;
    % solve Vector Difference Equation
    for j=(N-1):-1:1
        K_j = vec2mtx(K_history_vec(:,(n*n*j-15):n*n*j));
        g_j = A'*inv(eye(n)+K_j*B*inv(R)*B')*g_j - C'*Q*yr(:,i+j);
        g_history = [g_history g_j];
    end
    % time-varying feed-forward gain g(1) for RHTC
    [m,n] = size(g_history); 
    g_1 = g_history(:,n);
end