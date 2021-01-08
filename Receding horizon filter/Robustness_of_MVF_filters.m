% Systems and paramteres %
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

% Making big matrices for FIR filters
[B_bar, C_bar, G_bar, Xi] = MakeBigMatrices (A,B,C,G,Q,R,10);
H = inv(C_bar'*inv(Xi) * C_bar) * C_bar'*inv(Xi);

% Parameter initialize
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

% main procedure
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

% Plot
plot( 1:N_sample , real_state(2,:)-estimated_state(2,:), 1:N_sample , real_state(2,:)-FIR_estimated_state(2,:))