function [U,V] = REALUNV(imax,jmax,dt,dx,dy,UBAR,VBAR,P,VI)

% Description:
%
% This function performs the velocity corrector step using the fractional
% time-step velocity and the pressure. The corrector step is defined as
% given by
%
%    v^(n+1) = -dt*grad(P) + v_hat
%
% where the n+1 denotes the updated time step.

%% Initialize Storage
U = zeros(imax,jmax+1);
V = zeros(imax+1,jmax);

%% Compute U velocity
for i = 2:imax-1
    for j = 2:jmax
        U(i,j) = UBAR(i,j) - (P(i+1,j)-P(i,j))*dt/dx;
    end
end

%% Apply U BC's
U(1,2:jmax) = 0;                            % Left
U(imax,2:jmax) = 0;                         % Right (moving)
U(2:imax-1,1) = -U(2:imax-1,2);             % Bottom
U(2:imax-1,jmax+1) = -U(2:imax-1,jmax);     % Top

%% Compute V velocity
for i = 2:imax
    for j = 2:jmax-1
        V(i,j) = VBAR(i,j) - (P(i,j+1)-P(i,j))*dt/dy;
    end
end

%% Apply V BC's
V(2:imax,1) = 0;                                % Bottom
V(2:imax,jmax) = 0;                             % Top
V(1,2:jmax-1) = -V(2,2:jmax-1);                 % Left
V(imax+1,2:jmax-1) = 2*VI - V(imax,2:jmax-1);   % Right (moving)

return
end


