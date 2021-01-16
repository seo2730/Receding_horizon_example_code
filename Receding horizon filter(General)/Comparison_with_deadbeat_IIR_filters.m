% Systems and paramteres %
A = [0 1 0 1; 
     0 0 1 0; 
     0 0 0 0; 
     0 0 0 1];

B = [0 0;
     1 0;
     0 1;
     0 0];

G = [1 0 0 0]';

C = [1 0 0 0; 
     0 1 0 0];

N = [0 0 1 0; 
     0 0 0 1];

W = [-1 0 0 1;
      0 0 0 1;
      0 0 1 0];

Ao = [0 0 0;
      1 0 0;
      0 0 0];

L = [0 -1;
     1  0;
     0  0];
 
M = [0 0 1;
     0 1 0];
 
Q = 0.02^2; 
R = 0.04*eye(2);
N_sample = 250; 
N_order = size(A,1); 
N_horizon = 10;

% Making big matrices for FIR filters
[B_bar, C_bar, G_bar] = MakeBigMatrices2(A,B,C,G,Q,R,10);
% H = inv(C_bar'*inv(Xi) * C_bar) * C_bar'*inv(Xi);