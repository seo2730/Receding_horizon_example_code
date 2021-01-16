# Recedgin horizon LQ tracking control(Free terminal cost) 

## Function drhlqtc.m

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

## Function drde2.m

![image](https://user-images.githubusercontent.com/42115807/104813997-dc807980-584f-11eb-9799-a6a251458cf0.png)<br>

#### Forward computation

![image](https://user-images.githubusercontent.com/42115807/104814077-50bb1d00-5850-11eb-8f97-4cc357d67abb.png)<br>


    function [K_N_1,K_history_vec] = drde2(A,B,C,Q,R,Qf,N)
        % DRDE2 Discrete-time Riccati Difference Equation Solver
        % This RDE appears in LQ Tracking Problem

        % system dimension
        n = size(A,1);
        % boundary condition
        K_0 = C'*Qf*C;
        K_0_vec = mtx2vec(K_0);
        K_history_vec = K_0_vec;
        K_i = K_0;
        % solve Riccati Differential Equation 2
        for i=1:N-1
            K_i = A'*K_i*inv(eye(n)+B*inv(R)*B'*K_i)*A + C'*Q*C;
            K_i_vec = mtx2vec(K_i);
            K_history_vec = [K_history_vec K_i_vec];
        end
        % constant feedback gain K(N-1) for RHTC
        [s1,s2] = size(K_history_vec);
        K_N_1 = vec2mtx(K_history_vec(:,(s2-n*n+1):s2));
    end

## Function dvde.m

    function [g_1,g_history ] = dvde(A,B,C,Q,R,Qf,N,K_history_vec,yr,i)
        %DVDE Discrete-time Vector Differential Equation Solver
        % system dimension
        n = size(A,1);
        [s1,s2] = size(K_history_vec);
        % boundary condition
        g_N = -C'*Qf*yr(:,i+N); g_history = g_N; g_j = g_N;
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

## Function mtx2vec.m

    function vec = mtx2vec(mtx)
        vec = reshape(mtx', 1, []);
        %vec = reshape(mtx, 1, [])';
    end

## Function vec2mtx.m

    function mtx = vec2mtx(vec)
        [s1,s2] = size(vec);
        n = sqrt(s2);
        mtx = reshape(vec',n,[]);
        %mtx = reshape(vec,(s1),[])';
    end
