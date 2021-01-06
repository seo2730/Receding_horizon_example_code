function [g_1,g_history] = dvde_hinf(A, B, Bw, C, Q, R, Rw, Qf, gamma_2, N, K_history_vec, yr, i)
    n = size(A,1);
    g_N = -C'*Qf*yr(:,i+N); 
    g_history = g_N; 
    g_j = g_N;
    for j = (N-1) : -1 : 1
        K_j = vec2mtx(K_history_vec(:,(n*n*j-15):n*n*j));
        Lambda = eye(n) + K_j * (B*inv(R)*B' - 1/gamma_2 * Bw * inv(Rw) * Bw');
        %g_j = A'*inv(eye(n)+K_j*B*inv(R)*B')*g_j - C'*Q*yr(:,i+j);
        
        con = Rw - 1/gamma_2*Bw'*K_j*Bw;
        
        con_eig = eig(con);
        for i=1:rank(con)
           if con_eig(i) > 0
              flag = 1;           
           else
               flag = 0;
           end
        end
        
        if flag == 1
            g_j = -A'*inv(Lambda)*g_j - C'*Q*yr(:,i+j);
            g_history = [g_history g_j];
        end
    end
    % time-varying feed-forward gain g(1) for RHTC
    [m,n] = size(g_history); 
    g_1 = g_history(:,n);
end