function [CM] = coeff_matrix(imax,jmax,dx,dy)

% Desciption:
% 
% This function defines the coefficient matrix for the pressure solver.
% The pressure equation is a  2D poisson equation with Neumman BC's
% on all boundaries.  To prevent the matrix from being singular an
% arbitrary point is set to zero relative pressure.  Thus, the solution at
% any point is w.r.t. this value.  Use of a sparse matrix for the
% coefficients speeds up the solution by over a factor of 5.
%
% The solution uses a typical 5 point stencil.  Edges and corners must be
% handled appropriately for Neumann BC's.

%% Define parameters
L = (imax-1)*(jmax-1); % Length of the matrix
jump = imax-1;

%% Initialize Storage
CM = sparse(L,L);

%% Fill in sparse matrix
ind = 0;
for j = 2:jmax      % Loop over all interior pressure points
    for i = 2:imax
        ind = ind+1;
        
        if(i == 2 && j == 2)                    % Bottom-left corner
            CM(ind,ind)   = -1/dx^2 - 1/dy^2;
            CM(ind,ind+1) = 1/dx^2;
            CM(ind,ind+jump) = 1/dy^2;
        elseif(i == 2 && j == jmax)             % Top-left corner
            CM(ind,ind)   = -1/dx^2 - 1/dy^2;
            CM(ind,ind+1) = 1/dx^2;
            CM(ind,ind-jump) = 1/dy^2;
        elseif(i == imax && j == 2)             % Bottom-right corner
            CM(ind,ind)   = -1/dx^2 - 1/dy^2;
            CM(ind,ind-1) = 1/dx^2;
            CM(ind,ind+jump) = 1/dy^2;
        elseif(i == imax && j == jmax)          % Top-right corner
            CM(ind,ind)   = -1/dx^2 - 1/dy^2;
            CM(ind,ind-1) = 1/dx^2;
            CM(ind,ind-jump) = 1/dy^2;
        elseif(i == 2)                          % Left edge
            CM(ind,ind)   = -1/dx^2 - 2/dy^2;
            CM(ind,ind+1) = 1/dx^2;
            CM(ind,ind+jump) = 1/dy^2;
            CM(ind,ind-jump) = 1/dy^2;
        elseif(i == jmax)                       % Right edge
            CM(ind,ind)   = -1/dx^2 - 2/dy^2;
            CM(ind,ind-1) = 1/dx^2;
            CM(ind,ind+jump) = 1/dy^2;
            CM(ind,ind-jump) = 1/dy^2;
        elseif(j == 2)                          % Bottom edge
            CM(ind,ind)   = -2/dx^2 - 1/dy^2;
            CM(ind,ind+1) = 1/dx^2;
            CM(ind,ind-1) = 1/dx^2;
            CM(ind,ind+jump) = 1/dy^2;
        elseif(j == jmax)                       % Top edge
            CM(ind,ind)   = -2/dx^2 - 1/dy^2;
            CM(ind,ind+1) = 1/dx^2;
            CM(ind,ind-1) = 1/dx^2;
            CM(ind,ind-jump) = 1/dy^2;
        else                                    % All others
            CM(ind,ind)   = -2/dx^2 - 2/dy^2;
            CM(ind,ind+1) = 1/dx^2;
            CM(ind,ind-1) = 1/dx^2;
            CM(ind,ind+jump) = 1/dy^2;
            CM(ind,ind-jump) = 1/dy^2;
        end
    end
end

%% Reverse Sign
CM = -CM;

%% Set corner element to be 0
CM(1,:) = 0;
CM(1,1) = CM(2,2);

return
end

        