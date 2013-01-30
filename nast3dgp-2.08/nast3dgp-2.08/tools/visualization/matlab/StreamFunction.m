function [psi,om]=StreamFunctionVorticity(x,y,flag,u,v)
%
% [psi,om]=StreamFunctionVorticity2D(x,y,u,v) returns
% for a 2D velocity field (u,v) on a cartesian grid
% the streamfunction (psi) and the vorticity (om)
% 
% The velocity values are given on a staggered mesh including ghostvalues, 
% based on the particular boundary conditions
%
% The grid is given by the x and y coordinates of the gridlines, which bound the
% mass control volumes. If [m,n]==size(u)==size(v) then
% 
%  x=x[1...n-1]  ,  y=y[1...m-1]
% 
% In addition the flag is used to 'flag-out' certain cells
% flag(i,j)!=0 denotes obstacle cells, where no fluid is present
% 
% we assume psi==0 On any boundary  
%

[m,n]=size(u) ;

%%%%%%%%%%%%
% Omega
%
om =zeros(m-1,n-1) ;

dudy=om ;
dvdx=om ;

for j=1:n-1
 dudy(:,j)=0.5 * ( (u(1:m-1,j+1)-u(1:m-1,j))  + (u(2:m,j+1)-u(2:m,j)) ) / (y(j+1)-y(j)) ;
end

for i=1:m-1
 dvdx(i,:)=0.5 * ( (v(i+1,1:n-1)-v(i,1:n-1))  + (v(i+1,2:n)-v(i,2:n)) ) / (x(i+1)-x(i)) ;
end
 
om=dudy-dvdx ;

%%%%%%%%%%%%
% Psi
%
psi=zeros(m,n-1) ;

for j=2:n-1
  f=flag(:,j-1) ;
  du=u(:,j)   ;
  dy=y(j)-y(j-1)     ;
  psi(:,j)=(f==0).*(psi(:,j-1)+du*dy) ;
end

return
