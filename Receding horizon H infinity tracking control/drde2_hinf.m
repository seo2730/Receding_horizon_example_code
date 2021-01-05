function [K_N_1,K_history_vec] = drde2_hinf (A, B, Bw, C, Q, R, Rw, Qf, gamma_2, N )
    % system dimension
    n = size(A,1);
    % boundary condition
    K_0 = C'*Qf*C;
    K_0_vec = mtx2vec(K_0); 
    K_history_vec = K_0_vec; 
    K_i = K_0;
    % solve Riccati Differential Equation
    for i=1:N-1
        if min (real (eig (Rw - 1/gamma_2*Bw'*K_i*Bw)))< 0
            error ('error');
        end
        Lambda = eye(n) + K_i * ( B*inv(R)*B' - 1 / gamma_2 * Bw * Rw^(-1) * Bw' ) ;
        K_i = A' * inv(Lambda) * K_i *A + C'*Q*C;
        K_i_vec = mtx2vec(K_i);
        K_history_vec = [K_history_vec K_i_vec];
    end
    % constant feedback gain K(N-1) for RHTC
    [s1,s2] = size(K_history_vec); 
    K_N_1 = vec2mtx(K_history_vec(:,(s2-n*n+1):s2));
end