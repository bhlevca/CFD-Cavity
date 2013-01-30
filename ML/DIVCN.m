function [DIV,CN] = DIVCN(imax,jmax,dt,dx,dy,U,V)

% Description:
%
% This function computes the divergence of the flowfield on the pressure 
% nodes of the mesh using the staggered velocities.  By definition the 
% divergence should be zero (machine precision) for all points.  This is 
% imposed by solving the Poisson style system for the pseudo-pressure.  %
% Also, the Courant number is computed.

%% Initialize Storage
DIV = zeros(imax+1,jmax+1);

%% Compute divergence at all nodes
for i = 2:imax
    for j =2:jmax
        DIV(i,j) = (1/dx)*(U(i,j)-U(i-1,j))+...
                   (1/dy)*(V(i,j)-V(i,j-1));
    end
end

%% Compute Courant number
CN = 0;
for i = 2:imax
    for j = 2:jmax
        Us = .5*(U(i,j)+U(i-1,j));
        Vs = .5*(V(i,j)+V(i,j-1));
        CR1 = abs(Us/dx)+abs(Vs/dy);
        if CR1 > CN
            CN = CR1;
        end
    end
end
CN = CN*dt;

return
end