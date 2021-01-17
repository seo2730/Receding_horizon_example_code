# Receding horizon H infinity tracking control

## Function drhtc_hinf.m

![image](https://user-images.githubusercontent.com/42115807/104828974-36b22680-58b2-11eb-8880-2638d405a0ee.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104829023-ea1b1b00-58b2-11eb-8a63-0c93d86fc5ad.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104829027-fa32fa80-58b2-11eb-80ef-c83fe88685ed.png)<br>

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

## Function drde2_hinf.m

![image](https://user-images.githubusercontent.com/42115807/104829041-146cd880-58b3-11eb-88fa-25bfd26aac3e.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104829050-19318c80-58b3-11eb-94c7-08581e1e3c24.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104829027-fa32fa80-58b2-11eb-80ef-c83fe88685ed.png)<br>

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
            Lambda = eye(n) + K_i * ( B*inv(R)*B' - 1 / gamma_2 * Bw * inv(Rw) * Bw' ) ;
            K_i = A' * inv(Lambda) * K_i *A + C'*Q*C;
            K_i_vec = mtx2vec(K_i);
            K_history_vec = [K_history_vec K_i_vec];
        end
        % constant feedback gain K(N-1) for RHTC
        [s1,s2] = size(K_history_vec); 
        K_N_1 = vec2mtx(K_history_vec(:,(s2-n*n+1):s2));
    end

## Function dvde_hinf.m
![image](https://user-images.githubusercontent.com/42115807/104829067-4847fe00-58b3-11eb-948b-77e0f19dc15f.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104829027-fa32fa80-58b2-11eb-80ef-c83fe88685ed.png)<br>
위 식들이 성립할려면 아래 식이 성립되어야한다.<br>
![image](https://user-images.githubusercontent.com/42115807/104829083-6d3c7100-58b3-11eb-973a-0320ff8a0b9b.png)<br>


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

# Result
![image](https://user-images.githubusercontent.com/42115807/104829127-21d69280-58b4-11eb-84d1-ff7f856a0d6e.png)
