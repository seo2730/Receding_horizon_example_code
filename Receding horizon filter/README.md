# Make Big Matrices
### MVFIR Filterss with nosingular A(Batch form : 행렬 크기가 horizon 갯수에 따라 다름)
![image](https://user-images.githubusercontent.com/42115807/103977753-3bddea00-51bd-11eb-8ca4-75b4b67e1d6c.png)<br>
![image](https://user-images.githubusercontent.com/42115807/103977762-426c6180-51bd-11eb-95a9-4191ec842d1e.png)<br>
![image](https://user-images.githubusercontent.com/42115807/103977771-49936f80-51bd-11eb-8122-a595894bfe57.png)<br>
![image](https://user-images.githubusercontent.com/42115807/103977780-4ef0ba00-51bd-11eb-942c-eef883b9fd59.png)<br>
위 식들에 대한 코드가<br>

1.B_bar 초기값<br>

    if B == 0
        IsInput = 0;
        B_bar=0;
    else
        B_bar = C * A_inv * B;
    end

2.C_bar 초기값<br>

    C_bar = C * A_inv;

3.G_bar 초기값<br>

    G_bar = -C * A_inv * G;

4.B_bar, C_bar, G_bar 행렬 업데이트

    for i = 2 : horizon
        A_inv_i = A_inv_i * A_inv ;
        if IsInput == 1
            B_bar = [B_bar -C_bar*A_inv*B; zeros(N_output, N_input*(i-1)) -C*A_inv*B ];
        end
            G_bar = [G_bar -C_bar*A_inv*G; zeros(N_output, N_system_noise*(i-1)) -C*A_inv*G ] ;
            C_bar = [C * A_inv_i ;C_bar];
            Q_stack = daug(Q_stack,Q); R_stack = daug(R_stack,R);
    end

![image](https://user-images.githubusercontent.com/42115807/103978196-464cb380-51be-11eb-952c-ace72d0537a9.png)<br>
위 식의 코드는<br>

     Xi = G_bar * Q_stack *G_bar' + R_stack ;
     
<br>
<br>

# Robustness of MVF filters

### System and parameters

![image](https://user-images.githubusercontent.com/42115807/103979047-3930c400-51c0-11eb-8af7-a48c907107ac.png)<br>
![image](https://user-images.githubusercontent.com/42115807/103979094-4e0d5780-51c0-11eb-913e-c95d0566c6f2.png)<br>
![image](https://user-images.githubusercontent.com/42115807/103980467-4a2f0480-51c3-11eb-935f-a59a22aa1904.png)<br>

    A = [0.9305      0  0.1107; 
         0.0077 0.9802 -0.0173; 
         0.0142      0  0.8953];
    B = [0.0217  0.2510; 
         0.0192 -0.0051; 
         0.0247  0.0030];
    C = [1 0 0; 
         0 1 0]; 
    G = [1 1 1]';  
    Q = 0.02; 
    R = 0.02*eye(2);
    D_1 = 0; D_2 = 0;
    N_sample = 250; 
    N_order = size(A,1); 
    N_horizon = 10;

<br>
### Making big matrices for FIR filters

    [B_bar, C_bar, G_bar, Xi] = MakeBigMatrices (A,B,C,G,Q,R,10);
    
![image](https://user-images.githubusercontent.com/42115807/103978498-ec002280-51be-11eb-9315-a933f14b74ad.png)<br>
    
    H = inv(C_bar'*inv(Xi) * C_bar) * C_bar'*inv(Xi);
 
<br>
### Parameter initialize

    intial_state= [0 0 0 ]'; 
    x = intial_state;
    IIR_x_hat = [0 0 0]'; 
    FIR_x_hat = IIR_x_hat;
    P = 0.0*eye(N_order); 
    real_state = zeros(N_order , N_sample);
    estimated_state = zeros(N_order, N_sample); 
    real_state(:,1)=x;
    estimated_state(:,1) = IIR_x_hat;
    FIR_estimated_state = zeros(N_order, N_sample);
    measurements = zeros(2*N_sample);
    delta_A = 0.1*[1 0 0;0 1 0;0 0 0.1 ];
    delta_C = 0.1*[0.1 0 0;0 0.1 0 ];
  
<br>
### Main procedure

    for i = 1 : N_sample-1
        if( i > 50 && i < 101 )
            x = (A+delta_A)*x + G*randn(1)*0.02;
            y = (C+delta_C)*x + randn(2,1)*0.04;
        else
            x = A*x + G*randn(1)*0.02;
            y = C*x + randn(2,1)*0.04;
        end
        real_state(:,i+1)=x;
        % IIR_filter : one step predicted estimate
        IIR_x_hat = A * IIR_x_hat + A * P *C' * inv( R + C * P * C') * ( y - C*IIR_x_hat);
        P = A * inv( eye(N_order) + P * C' * inv(R) * C) * P * A' + G * Q * G';
        estimated_state(:,i+1)=IIR_x_hat;
        measurements(2*i-1:2*i) = y;
        % FIR filter
        if i>10
            FIR_x_hat = H * (measurements(2*i-19:2*i))';
        end
            FIR_estimated_state(:,i+1) = FIR_x_hat;
    end
    
**코드 설명 :**<br>

    if( i > 50 && i < 101 )
        x = (A+delta_A)*x + G*randn(1)*0.02;
        y = (C+delta_C)*x + randn(2,1)*0.04;
    else
        x = A*x + G*randn(1)*0.02;
        y = C*x + randn(2,1)*0.04;
    end

delta 조건 때문에 이런 식으로 코드를 짜보았다. 뒤에 random 함수는 noise를 추가하였다.<br>

    % IIR_filter : one step predicted estimate
    IIR_x_hat = A * IIR_x_hat + A * P *C' * inv( R + C * P * C') * ( y - C*IIR_x_hat);
    P = A * inv( eye(N_order) + P * C' * inv(R) * C) * P * A' + G * Q * G';
    estimated_state(:,i+1)=IIR_x_hat;
    measurements(2*i-1:2*i) = y;
    
IIR filter는 Kalma filter와 같다.<br>

    % FIR filter
    if i>10
        FIR_x_hat = H * (measurements(2*i-19:2*i))';
    end
        FIR_estimated_state(:,i+1) = FIR_x_hat;

FIR filter는 horizon 크기가 10이므로 10 이상이 되어야 추정이 가능하다.<br>
![image](https://user-images.githubusercontent.com/42115807/103981006-45b71b80-51c4-11eb-88f1-60635733339e.png)<br>
위 사진 코드가<br>

    FIR_x_hat = H * (measurements(2*i-19:2*i))'

# Result
![image](https://user-images.githubusercontent.com/42115807/103981197-b3634780-51c4-11eb-98ee-fa06c87cccdb.png)<br>
![image](https://user-images.githubusercontent.com/42115807/103981211-bb22ec00-51c4-11eb-916a-2f51d6de8e79.png)<br>
![image](https://user-images.githubusercontent.com/42115807/103981242-c8d87180-51c4-11eb-8ec4-e7da14e27992.png)<br>
![image](https://user-images.githubusercontent.com/42115807/103981259-d130ac80-51c4-11eb-8929-d2a22699c4ad.png)<br>
Blue line : Kalman Filter, Red line : FIR filter<br>
<br>
결과가 다른 이유는 random 함수 때문이다. 대체적으로 FIR가 robust한 것을 볼 수 있다

