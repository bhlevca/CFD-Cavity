function [VBAR,OLDL2] = SOLVEVBAR(imax,jmax,dt,dx,dy,Re,U,V,OLDL2,VI)

% Description:
%
% This function solves for the new fractional step y velocity.  The
% advective terms are interpolated onto nodes prior to differencing.
% Boundary conditions are imposed directly for any points on a boundary and
% using linear interpolation for those points staggered across a boundary.

%% Initialize Storage
VBAR = zeros(imax+1,jmax);
RL   = zeros(imax+1,jmax);

%% Compute simplifying factor
alpha = -dt/Re/(dx^2*dy^2);

%% Compute explicit terms
for i = 2:imax
    for j = 2:jmax-1
        B1 = .5*(U(i,j+1)+U(i,j)); % COEFFS FOR D(UV)/DX
        B2 = .5*(V(i+1,j)+V(i,j)); % COEFFS FOR D(UV)/DX
        B3 = .5*(U(i-1,j+1)+U(i-1,j)); % COEFFS FOR D(UV)/DX
        B4 = .5*(V(i-1,j)+V(i,j)); % COEFFS FOR D(UV)/DX
        B5 = .5*(V(i,j+1)+V(i,j)); % TIME ADVANCE V(I,J+1)
        B6 = .5*(V(i,j-1)+V(i,j)); % TIME ADVANCE V(I,J)        
        RL(i,j) = -((B1*B2-B3*B4)/dx + (B5^2-B6^2)/dy);
    end
end

%% Compute fractional step velocity
for i = 2:imax
    for j = 2:jmax-1
        VBAR(i,j) = V(i,j)+ dt*RL(i,j) - (dt/2)*OLDL2(i,j) - ...
            alpha*(dy^2*(V(i+1,j)-2*V(i,j)+V(i-1,j)) + ...
                   dx^2*(V(i,j+1)-2*V(i,j)+V(i,j-1)));
    end
end

%% Impose BC's
VBAR(2:imax,   1) = 0;      % Bottom
VBAR(2:imax,jmax) = 0;      % Top
VBAR(1,     1:jmax) = -VBAR(2,1:jmax);          % Left
VBAR(imax+1,1:jmax) = 2*VI - VBAR(imax,1:jmax); % Right (moving wall)

%% UPDATE OLDL2
OLDL2(2:imax,2:jmax-1) = RL(2:imax,2:jmax-1);

return
end


