M = 1; bb = 0.01; Ac = [0 1 0 0;0 -bb/M 0 0;0 0 0 1;0 0 0 -bb/M];
Bc = [0 0 ;1/M 0 ; 0 0;0 1/M]; Cc = [1 0 0 0 ;0 0 1 0]; Dc = [0 0;0 0]; 
Ts = 0.055; [A B C D] = c2dm(Ac,Bc,Cc,Dc,Ts,'zoh');
% design parameter in cost function
% - state weighting matrix: Q
% - control weighting matrix: R
% - initial state: x0
Q = eye(2); R = eye(2); Qf = eye(2); x0 = [1;0;1;0]; N=10;
% arbitrary future references available over [i,i+N]
y1r_1 = 1:-0.01:0.01; y1r_2 = zeros(1,100); y1r_3 = 0:0.01:0.99;
y1r_4 = ones(1,100); y2r_1 = ones(1,100); y2r_2 = 1:-0.01:0.01;
y2r_3 = zeros(1,100); y2r_4 = 0:0.01:0.99; y1r = [y1r_1 y1r_2 y1r_3 y1r_4 ones(1,100)]; 
y2r = [y2r_1 y2r_2 y2r_3 y2r_4 ones(1,100)]; yr = [y1r;y2r];
% simulation step
is = 440;
% Discrete-time LQTC for Unconstrained Systems
[x,y,u] = drhlqtc(x0,A,B,C,Q,R,Qf,N,yr,is);