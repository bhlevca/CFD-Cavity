function [P] = SOLVEP(imax,jmax,dt,dx,dy,RHP,UBAR,VBAR,CMP)

% Description:
%
% This function solves the pressure equation.  The coefficient matrix is
% given by CMP and the RHS of the equation is stored in the matrix RHP.
% This matrix is transformed into a vector prior to solving the system.


%% Define parameters
L = (imax-1)*(jmax-1);

%% Initialize Storage
F  = zeros(L,1);
P  = zeros(imax+1,jmax+1);

%% Transpose RHS into vector
K = 1;
for j = 2:jmax
    for i = 2:imax
        F(K) = RHP(i,j);
        K=K+1;
    end
end

%% Impose pressure at corner = 0
F(1) = 0;

%% Solve system
x3= CMP\F;

%% Put solution into matrix form
Q = 1;
for j = 2:jmax
    for i = 2:imax
        P(i,j) = x3(Q);
        Q=Q+1;
    end
end

%% Take care of boundary conditions
P(2:imax,1) = P(2:imax,2);          % Left
P(2:imax,jmax+1) = P(2:imax,jmax);  % Right
P(1,2:jmax) = P(2,2:jmax);          % Bottom
P(imax+1,2:jmax) = P(imax,2:jmax);  % Top

return
end
