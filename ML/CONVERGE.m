function [cmax] = CONVERGE(U,OLDU)

% Description:
%
% This function determines when the system has converged to a steady
% solution. Currently, the code is set up to compute the infinity norm of
% the difference between the solution at the current time step and the
% previous time step.

%% Compute infinity norm (treat solution as a vector, not a matrix)
cmax = max(max(abs(U-OLDU)));

return
end

