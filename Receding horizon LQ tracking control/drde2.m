function [K_N_1,K_history_vec] = drde2(A,B,C,Q,R,Qf,N)
    % DRDE2 Discrete-time Riccati Difference Equation Solver
    % This RDE appears in LQ Tracking Problem

    % system dimension
    n = size(A,1);
    % boundary condition
    K_0 = C'*Qf*C; K_0_vec = mtx2vec(K_0); K_history_vec = K_0_vec;
    K_i = K_0;
    % solve Riccati Differential Equation 2
    for i=1:N-1
        K_i = A'*K_i*inv(eye(n)+B*inv(R)*B'*K_i)*A + C'*Q*C;
        K_i_vec = mtx2vec(K_i);
        K_history_vec = [K_history_vec K_i_vec];
    end
    % constant feedback gain K(N-1) for RHTC
    [s1,s2] = size(K_history_vec); 
    K_N_1 = vec2mtx(K_history_vec(:,s2));
end