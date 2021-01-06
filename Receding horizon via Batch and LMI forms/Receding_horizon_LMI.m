A = [ 0.6831 0.0353 ; 0.0928 0.6124 ] ; B = [0.6085;0.0158];
Q =eye(2); R =3;
SystemDim = size(A,1); InputDim = size(B,2);
% Getting started LMISYS description

setlmis([]);
% Defining LMI variables
Y=lmivar(2,[InputDim SystemDim]); X=lmivar(1,[SystemDim 1]);

% Defining LMI
lmiterm([-2 1 1 X],1,1); % LMI (1,1) : X
lmiterm([-2 2 1 Y],B,1); % LMI (2,1) : B*Y
lmiterm([-2 2 1 X],A,1); % LMI (2,1) : A*X
lmiterm([-2 2 2 X],1,1); % LMI (2,2) : X
lmiterm([-2 3 1 X],sqrt(Q),1); % LMI (3,1) : Q^(1/2)*X
lmiterm([-2 3 3 0],1); % LMI (3,3) : 1
lmiterm([-2 4 1 Y],sqrt(R),1); % LMI (4,1) : R^(1/2)*Y
lmiterm([-2 4 4 0],1); % LMI (4,4) : 1

LMISYS = getlmis; 

[tmin,xfeas]=feasp(LMISYS);

X=dec2mat(LMISYS,xfeas,X);

H = Y*inv(X);
Qf = inv(X);