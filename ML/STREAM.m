function [psi] = STREAM(imax,jmax,dx,dy,U,V)

% Description:
%
% This function solves for the stream function using the V velocity. The
% streamfunction is defined as follows:
%
%    u =  d(psi)/dy
%    v = -d(psi)/dx
%
% The values are computed on the non-staggered mesh.
% The velocity is the velocity on the node.


%% Initialize Variables
psi = zeros(imax,jmax);

%% Compute stream funciton
for i =2:imax-1
    psi(i,2:jmax-1) = psi(i-1,2:jmax-1) - dx*V(i,2:jmax-1);
end

return
end