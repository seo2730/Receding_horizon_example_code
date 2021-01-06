clear;
% state-space model of original system and its discretization
M = 1; bb = 0.01; Ac = [0 1 0 0;0 -bb/M 0 0;0 0 0 1;0 0 0 -bb/M];
Bc = [0 0; 1/M 0; 0 0; 0 1/M]; Cc = [1 0 0 0; 0 0 1 0]; Dc = [0 0; 0 0]; Ts = 0.055;
[A, B, C, D] = c2dm(Ac,Bc,Cc,Dc,Ts,'zoh'); 
Bw = [0.016 0.002; 0.01 0.009; 0.008 0; 0 0.0005];

% design parameter in cost function
% - state weighting matrix: Q
% - control weighting matrix: R
% - disturbance weighting matrix: Rw
% - terminal weighting matrix: Qf
% - prediction horizon: N
system_order = size ( A , 1); input_order = size ( B , 2 );
dist_order = size ( Bw , 2 );
Q = eye(system_order); R = eye(input_order); Rw = 1.2*eye(dist_order); N = 5;

% arbitrary future references available over [i,i+N]
y1r_1 = 1:-0.01:0.01; y1r_2 = zeros(1,100); y1r_3 = 0:0.01:0.99;
y1r_4 = ones(1,100); y2r_1 = ones(1,100); y2r_2 = 1:-0.01:0.01;
y2r_3 = zeros(1,100); y2r_4 = 0:0.01:0.99; y1r = [y1r_1 y1r_2 y1r_3 y1r_4 ones(1,100)]; y2r = [y2r_1 y2r_2 y2r_3 y2r_4 ones(1,100)]; yr = [y1r;y2r];

% simulation step
is = 440;
% Discrete-time RHTC for Unconstrained Systems
% initial state
x0 = [ 1; 0; 1; 0]; Qf = 100*eye(2); Q = eye(2);
gamma_2 = 10; 
[x_hinf,y_hinf,u_hinf,w] = drhtc_hinf (x0, A, B, Bw, C, Q, R, Rw, Qf, gamma_2, N, yr, is);
[x,y,u] = drhlqtc(x0,A,B,C,Q,R,Qf,N,yr,is);

for k = 0 : is
    y_rhtc_ws(k+1,:) = [k*0.01 y(:,k+1)'];
end

for k = 0 : is
    r_rhtc_ws(k+1,:) = [k*0.01 yr(:,k+1)'];
end

for k = 0 : is
    y_rhtc_ws_hinf(k+1,:) = [k*0.01 y_hinf(:,k+1)'];
end

plot(r_rhtc_ws(:,2) , r_rhtc_ws(:,3) , y_rhtc_ws(:,2) ,y_rhtc_ws(:,3) ,'-.' , y_rhtc_ws_hinf(:,2) , y_rhtc_ws_hinf(:,3));
legend('Reference','LQ RHTC','H infinity RHTC') ;


