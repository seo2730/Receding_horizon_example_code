# Receding horizon filter(General)

## MakeBigMatrices2.m
![image](https://user-images.githubusercontent.com/42115807/104829730-dd4df580-58b9-11eb-8a3c-8afe8e27cd2d.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104829736-e63ec700-58b9-11eb-82c2-9e1bfe9a2b63.png)<br>
![image](https://user-images.githubusercontent.com/42115807/104829751-138b7500-58ba-11eb-90c2-e02632a4180d.png)<br>

    for i = 1 : horizon
        if i==1
            if IsInput == 1
                B_bar = [zeros(N_output, N_input)];
            end
            G_bar = [zeros(N_output,N_system_noise)];
            C_bar = C;
            
        elseif i==2
            if IsInput == 1
                B_bar = [B_bar zeros(N_output, N_input);  C*B  zeros(N_output, N_input)];
                B_i = C*B;
            end
            G_bar = [G_bar zeros(N_output, N_system_noise); C*G zeros(N_output, N_system_noise)];
            C_bar = [C_bar; C * A_i];
            B_H = [A_i*G B_H];
            
            G_i = C*G;
            
            Q_stack = daug(Q_stack,Q); 
            R_stack = daug(R_stack,R);
        else
            if IsInput == 1
                B_bar = [B_bar zeros(N_output*(i-1), N_input);  C*A_i_2*B B_i zeros(N_output, N_input)];
                B_i = [C*A_i_2*B B_i];
            end
            G_bar = [G_bar zeros(N_output*(i-1), N_system_noise); C*A_i_2*G G_i zeros(N_output, N_system_noise)];
            C_bar = [C_bar; C * A_i];
            B_H = [A_i*G B_H];
            
            G_i = [C*A_i_2*G G_i];
            
            Q_stack = daug(Q_stack,Q); 
            R_stack = daug(R_stack,R);
            
            A_i_2 = A_i_2 * A;
        end
        
            A_i = A_i * A;      
    end

![image](https://user-images.githubusercontent.com/42115807/104829818-edb2a000-58ba-11eb-8409-c400c5a08f70.png)<br>

    W11 = C_bar'*inv(R_stack)*C_bar;
    W12 = C_bar'*inv(R_stack)*G_bar;
    W22 = G_bar'*inv(R_stack)*G_bar + inv(Q_stack);
    
    W = [W11  W12;
         W12' W22];
         
![image](https://user-images.githubusercontent.com/42115807/104829842-4b46ec80-58bb-11eb-98e0-223d66b741b9.png)<br>
반복문에 위 식을 넣고 마지막에 아래 코드를 추가한다.<br> 
     
     B_H = [A_i B_H];

## Robustness_of_MVF_filter.m
    [B_bar, C_bar, G_bar, B_H, W22, R_N] = MakeBigMatrices2(A,B,C,G,Q,R,N_horizon);
위 함수로 각 변수들을 구한 뒤<br>
![image](https://user-images.githubusercontent.com/42115807/104829866-a8db3900-58bb-11eb-9a8d-a5889e1c6129.png)<br>
위 식을 아래 코드 표현 가능하다.<br>
    
    CG = [C_bar'; G_bar'];
    H = B_H*inv(W22)*CG*inv(R_N);

## Result
#### Reference와 FIR 필터 IIR 필터(Kalman Filter) 비교
![image](https://user-images.githubusercontent.com/42115807/104829510-f2c22000-58b7-11eb-8f5f-3dc8e973a455.png)<br>

#### Reference와 FIR 필터 IIR 필터(Kalman Filter) error 비교
![image](https://user-images.githubusercontent.com/42115807/104829529-269d4580-58b8-11eb-9828-526f99b64cf1.png)
