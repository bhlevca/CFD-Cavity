function [RHP] = RHSP(imax,jmax,dt,dx,dy,UBAR,VBAR)

% Description:
%
% This function computes the RHS of the pressure Poisson equation.  The RHS
% is simple to calculate where RHS = div(v_hat)/dt.  Although the RHS
% should be in vector form to solve the linear system, we will save it as a
% matrix here for ease of computing and then transfer it into a vector
% before solving the linear system.

%% Initialize Storage
RHP = zeros(imax+1,jmax+1);

%% Compute RHS using velocity components on staggered mesh
for i = 2:imax
    for j = 2:jmax
        RHP(i,j) = (UBAR(i,j)-UBAR(i-1,j))/dt/dx + ...
                   (VBAR(i,j)-VBAR(i,j-1))/dt/dy;
    end
end

%% Change sign
RHP = -RHP;

return
end