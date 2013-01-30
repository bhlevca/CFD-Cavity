function [] = RESULTS_PLOT(x,y,UP,VP,psi,VOR,P)

% Description:
%
% This function plots results for the driven cavity flow.

% Plot velocity vectors
figure; 
quiver(x,y,UP',VP','k');
title('Velocity Vectors')
xlabel('x')
ylabel('y')
axis equal
axis([0 1 0 1])

% Plot stream function contours
figure
contour(x,y,psi',15,'k')
axis equal
axis ([0 1 0 1]);
xlabel('x')
ylabel('y')
title('Stream-function contours')

% Plot vorticity contours
figure
contour(x,y,VOR',100,'k')
axis equal
axis([0 1 0 1]);
xlabel('x')
ylabel('y')
title('Vorticity contours')

% Plot pressure contours
x_p = zeros(length(x)-1,1);
for i = 1:length(x)-1
    x_p(i) = (x(i) + x(i+1))/2;
end

y_p = zeros(length(y)-1,1);
for i = 1:length(y)-1
    y_p(i) = (y(i) + y(i+1))/2;
end

figure
contourf(x_p,y_p,P(2:end-1,2:end-1)',20)
axis equal
axis([0 1 0 1]);
xlabel('x')
ylabel('y')
title('Pressure contours')


return
end