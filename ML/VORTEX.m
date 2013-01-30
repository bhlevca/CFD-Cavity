function [VOR] = VORTEX(imax,jmax,dx,dy,U,V)

% Description:
%
% This function computes the vorticity on the grid nodes using the velocity
% componenents on the staggered mesh.

%% Initialize Storage
VOR = zeros(imax,jmax);

%% Compute Vorticity
for i=1:imax
    for j=1:jmax
        VOR(i,j)=((V(i+1,j)-V(i,j))/dx)-...
                 ((U(i,j+1)-U(i,j))/dy);
    end
end

%% Deal with corners by taking average of adjacent two points
VOR(1,1) = .5*(VOR(1,2)+VOR(2,1));
VOR(1,jmax) = .5*(VOR(2,jmax)+VOR(1,jmax-1));
VOR(imax,1) = .5*(VOR(imax-1,1)+VOR(imax,2));
VOR(imax,jmax) = .5*(VOR(imax-1,jmax)+VOR(imax,jmax-1));

return
end

