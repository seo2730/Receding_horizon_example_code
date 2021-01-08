function [B_bar , C_bar , G_bar] = MakeBigMatrices2(A,B,C,G,Q,R,horizon)
    C_bar = C; 
    G_bar = zeros(2,1);
    Q_stack = Q ; R_stack = R ; A_i = A; A_i_2 = A;
    % parameter setting
    N_input = size(B,2); N_output = size(C,1);
    N_system_noise = size(G,2); N_order = size(A,1);
    
    % Error check
    IsInput = 1;
    % initialization
    if B == 0
        IsInput = 0;
        B_bar = 0;
    else
        B_bar = zeros(N_output,N_input); 
    end
    
    % main procedure
    for i = 1 : horizon
        if i==1
            if IsInput == 1
                B_bar = [zeros(N_output, N_input)];
            end
            G_bar = [zeros(N_output,N_system_noise)];
            C_bar = C;
            
        elseif i==2
            B_bar = [B_bar zeros(N_output, N_input);  C*B  zeros(N_output, N_input)];
            G_bar = [G_bar zeros(N_output, N_system_noise); C*G zeros(N_output, N_system_noise)];
            C_bar = [C_bar; C * A_i];
            
            B_i = C*B;
            G_i = C*G;
        else
            B_bar = [B_bar zeros(N_output*(i-1), N_input);  C*A_i_2*B B_i zeros(N_output, N_input)];
            G_bar = [G_bar zeros(N_output*(i-1), N_system_noise); C*A_i_2*G G_i zeros(N_output, N_system_noise)];
            C_bar = [C_bar; C * A_i];
            
            B_i = [C*A_i_2*B B_i];
            G_i = [C*A_i_2*G G_i];
            
            A_i_2 = A_i_2 * A;
        end
            Q_stack = daug(Q_stack,Q); 
            R_stack = daug(R_stack,R);
            A_i = A_i * A;
    end
    %Xi = G_bar * Q_stack *G_bar' + R_stack;
end