# Receding horizon LQ tracking control(Free terminal cost) 

## Function drhlqtc.m

![image](https://user-images.githubusercontent.com/42115807/104814439-3124f400-5852-11eb-943e-16ea962ce2d7.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104814404-05a20980-5852-11eb-9f0b-e477d01ee916.png)<br>
-> 위 식은 Matrix inverse 공식 이용하면 된다.<br>
-> yr는 reference value이다.

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

![image](https://user-images.githubusercontent.com/42115807/104814182-d2ab4600-5850-11eb-8453-f58a5fc784c3.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104814172-bc9d8580-5850-11eb-9fcc-390a1c7b524f.png)<br>

#### Forward computation

![image](https://user-images.githubusercontent.com/42115807/104814077-50bb1d00-5850-11eb-8f97-4cc357d67abb.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104814222-100fd380-5851-11eb-8c6b-36019ea67334.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104814232-21f17680-5851-11eb-86fc-c4ee970f9fb4.png)<br>

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

![image](https://user-images.githubusercontent.com/42115807/104814263-4c433400-5851-11eb-9780-dee1b854d1e0.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104814296-7d236900-5851-11eb-9515-aeec34cf3307.png)<br>

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
Matrix에서 Vector로 변경해주는 함수이다.<br>

    function vec = mtx2vec(mtx)
        vec = reshape(mtx', 1, []);
        %vec = reshape(mtx, 1, [])';
    end

## Function vec2mtx.m
Vector에서 Matrix로 변경해주는 함수이다.<br>

    function mtx = vec2mtx(vec)
        [s1,s2] = size(vec);
        n = sqrt(s2);
        mtx = reshape(vec',n,[]);
        %mtx = reshape(vec,(s1),[])';
    end
    
## Result
![image](https://user-images.githubusercontent.com/42115807/104829106-de7c2400-58b3-11eb-9323-cd954507cbcc.png)

