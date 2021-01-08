function [B_bar , C_bar , G_bar , Xi] = MakeBigMatrices(A,B,C,G,Q,R,horizon)
    % Error check
    A_inv = inv(A); IsInput = 1;
    % initialization
    if B == 0
        IsInput = 0;
        B_bar=0;
    else
        B_bar = C * A_inv * B;
    end
    
    C_bar = C * A_inv; 
    G_bar = -C * A_inv * G;
    Q_stack = Q ; R_stack = R ; A_inv_i = A_inv;
    % parameter setting
    N_input = size(B,2); N_output = size(C,1);
    N_system_noise = size(G,2); N_order = size(A,1);
    % main procedure
    for i = 2 : horizon
        A_inv_i = A_inv_i * A_inv ;
        if IsInput == 1
            B_bar = [B_bar -C_bar*A_inv*B; zeros(N_output, N_input*(i-1)) -C*A_inv*B ];
        end
            G_bar = [G_bar -C_bar*A_inv*G; zeros(N_output, N_system_noise*(i-1)) -C*A_inv*G ] ;
            C_bar = [C * A_inv_i ;C_bar];
            Q_stack = daug(Q_stack,Q); 
            R_stack = daug(R_stack,R);
    end
    Xi = G_bar * Q_stack *G_bar' + R_stack ;
end