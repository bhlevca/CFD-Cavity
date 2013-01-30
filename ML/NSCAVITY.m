% This routine calculates Navier-Stokes equation for lid-driven cavity
% flow in  a square cavity. Vorticity and stream function plots are output 
% in subprogram, RESULTS_PLOT.m

%======================================================================
%                        Main function
%======================================================================
clc
clear all
close all
format long

% dt/dx/dy : Time/Length Distance
% UBAR/VBAR : Velocity Components Before Correction
% U/V : Velocity Components after Correction 
% CMP : Coefficient Matrix for Pressure solver
% RHP : Right hand side matrix for the Pressure solver
% P : Pressure
% VI : Lid moving velocity
% CN : Courant Number
% d : Diffusive Term
% psi : Stream function psi(x,y)
% cmax : Convergence maximum

%% Define inputs
imax = 50;      % Grid size in x direction
jmax = 50;      % Grid size in y direction
ite  = 20000;   % Max number of iteration
err  = 1e-7;    % Error
Re   = 100;       % Reynolds Number
VI   = 1;       % Wall Velocity V

%% Compute some factors
dx = 1/(imax-1);
dy = 1/(jmax-1);
x  = 0:dx:1;      % x direction dimension
y  = 0:dy:1;      % y direction dimension
dt = .025*Re*dx^2;

%% --------- Initialize Storage -----------------

% Pressure ----------------------
P   = zeros(imax+1,jmax+1);
RHP = zeros(imax+1,jmax+1);
DIV = zeros(imax+1,jmax+1);

% Vertical Velocity -------------
VBAR  = zeros(imax+1,jmax);
V     = zeros(imax+1,jmax);
OLDL2 = zeros(imax+1,jmax);

% Lid velocity on the right side
VBAR(imax+1,:) = 2*VI - VBAR(imax,:);

% Horizontal Velocity -----------
UBAR  = zeros(imax,jmax+1);
U     = zeros(imax,jmax+1);
OLDL1 = zeros(imax,jmax+1);

% Velocity on colocated mesh ---
UP = zeros(imax,jmax);
VP = zeros(imax,jmax);


%% ---------- Compute coefficient Matrix -----------------
[CMP] = coeff_matrix(imax,jmax,dx,dy);

%% ---------- Begin time integration -------------------
TOTAL_t = 0; % Total time spent
for i = 1:ite
    TOTAL_t = TOTAL_t + dt;
    [UBAR,OLDL1] = SOLVEUBAR(imax,jmax,dt,dx,dy,Re,U,V,OLDL1);
    [VBAR,OLDL2] = SOLVEVBAR(imax,jmax,dt,dx,dy,Re,U,V,OLDL2,VI);
    [RHP] = RHSP(imax,jmax,dt,dx,dy,UBAR,VBAR);
    [P] = SOLVEP(imax,jmax,dt,dx,dy,RHP,UBAR,VBAR,CMP);
    OLDU = U;
    [U,V] = REALUNV(imax,jmax,dt,dx,dy,UBAR,VBAR,P,VI);
    [DIV,CN] = DIVCN(imax,jmax,dt,dx,dy,U,V);
    if i > 10
        [cmax] = CONVERGE(U,OLDU);
        disp(['It = ',int2str(i),...
              '; CN = ',num2str(CN),...
              '; CMAX = ',num2str(cmax)])
        if cmax < err
            break;
        end
    end
end
Num_Iteration = i;

disp(' ')
disp(['Total Iterations = ',int2str(i),...
      '; total time = ',num2str(TOTAL_t),...
      '; CN = ',num2str(CN)])


%% Compute velocity on nodes
for j = 1:jmax
    for i = 1:imax
        UP(i,j) = 0.5*(U(i,j+1) + U(i,j));
        VP(i,j) = 0.5*(V(i+1,j) + V(i,j));
    end
end

%% Compute stream function and vorticity
[psi] = STREAM(imax,jmax,dx,dy,UP,VP);
[VOR] = VORTEX(imax,jmax,dx,dy,U,V);

%% Plot results
RESULTS_PLOT(x,y,UP,VP,psi,VOR,P)
