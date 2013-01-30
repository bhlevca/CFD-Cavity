function [UBAR,OLDL1] = SOLVEUBAR(imax,jmax,dt,dx,dy,Re,U,V,OLDL1)

% Description:
%
% This function solves for the new fractional step x velocity.  The
% advective terms are interpolated onto nodes prior to differencing.
% Boundary conditions are imposed directly for any points on a boundary and
% using linear interpolation for those points staggered across a boundary.

%% Initialize Storage
UBAR = zeros(imax,jmax+1);
RL   = zeros(imax,jmax+1);

%% Compute simplification parameter
alpha = -dt/Re/(dx^2*dy^2);

%% Compute explicit terms
for i = 2:imax-1
    for j = 2:jmax
        A1 = .5*(U(i+1,j)+U(i,j)); % TIME ADVANCE U(I+1,J)
        A2 = .5*(U(i-1,j)+U(i,j)); % TIME ADVANCE U(I,J)
        A3 = .5*(U(i,j+1)+U(i,j)); % COEFFS FOR D(UV)/DY
        A4 = .5*(U(i,j-1)+U(i,j)); % COEFFS FOR D(UV)/DY
        A5 = .5*(V(i+1,j)+V(i,j)); % COEFFS FOR D(UV)/DY
        A6 = .5*(V(i,j-1)+V(i+1,j-1)); % COEFFS FOR D(UV)/DY
        RL(i,j) = -((A1^2-A2^2)/dx + (A3*A5-A4*A6)/dy);
    end
end

%% Compute fractional step velocity
for i = 2:imax-1
    for j = 2:jmax
        UBAR(i,j) = U(i,j) + dt*RL(i,j) - (dt/2)*OLDL1(i,j) - ...
            alpha*(dy^2*(U(i+1,j)-2*U(i,j)+U(i-1,j)) + ...
                   dx^2*(U(i,j+1)-2*U(i,j)+U(i,j-1)));
    end
end

%% Impose BC's
UBAR(1   ,2:jmax) = 0;  % Left
UBAR(imax,2:jmax) = 0;  % Right
UBAR(1:imax,     1) = -UBAR(1:imax,2);      % Bottom
UBAR(1:imax,jmax+1) = -UBAR(1:imax,jmax);   % Top

%% UPDATE OLDL1
OLDL1(2:imax-1,2:jmax) = RL(2:imax-1,2:jmax);

return
end
