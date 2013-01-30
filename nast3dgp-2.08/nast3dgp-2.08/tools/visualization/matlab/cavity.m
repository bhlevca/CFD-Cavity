%
% small example how to visualize with matlab
%

% read velocity components
u=ReadNaSt3D('t.u') ;
v=ReadNaSt3D('t.v') ;
p=ReadNaSt3D('t.p') ;
f=ReadNaSt3D('t.fl') ;

% extract 2D slices
%u=reshape(U.a(:,:,4),[U.n(1),U.n(2)]) ;
%v=reshape(V.a(:,:,4),[U.n(1),U.n(2)]) ;
%f=reshape(F.a(:,:,4),[U.n(1),U.n(2)]) ;

% get gridlines
%x=U.x ;
%y=V.y ;

% streamfunction
[psi,om]=StreamFunction(u.x,u.y,f.a(:,:,5),u.a(:,:,5),v.a(:,:,5)) ;

% visualize result
%[xg,yg]=ndgrid(x,y(1:U.n(2)-1)) ;
%contour(xg,yg,psi) ;
%contour(xg,yg,om) ;


% plot
figure (1); 
%pcolor (psi(2:m,1:m-1)); 
%l=[ -0.09 -0.08 -0.075 -0.0725 -0.07 -0.0675 -0.065 -0.0625 -0.06
%-0.0575 -0.055  -0.525 -0.05 -0.04 -0.03 -0.02 -0.01 -0.005 -0.0025
%-0.00125 -0.000675 -0.0003375 -0.00016875 -8.4375e-5 -4.21875e-5 0
%4.21875e-5 ];
l=[ -1e-10 -1e-7 -1e-5 -1e-4 -0.01 -0.03 -0.05 -0.07 -0.09 -0.1 -0.11 -0.115 -0.1175 1e-8 1e-7 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.3e-3 3e-3];
[c,h] = contour (psi', l); colorbar;
%[c,h] = contour ((psi(2:m,1:m-1))'); colorbar;
%clabel(c,h);
maxindex = 50;
axis ([1 maxindex-1 1 maxindex-1]);
title ('Stream Function');

figure (2); 
%l = [ -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 -0.5 -0.25 -0.125 0 0.125 0.25 0.5
%1 2 3 4 5 6 7 8 9 10];
l = [ -3.0 -2.0 -1.0 -0.5 0 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0];
[c,h] = contour (om', l); colorbar;
%clabel(c,h);
axis ([1 maxindex-1 1 maxindex-1]);
title ('Vorticity');
%omega ((maxindex - 1)/2, (maxindex - 1)/2)

figure (3); 
s=4;
%quiver (v,u,s); 
quiver (u.a(:,:,5)',v.a(:,:,5)',s);
axis ([1 maxindex 1 maxindex]);
title ('Velocity Field');

figure (4); 
%l=-0.25:1/25:1.25;
lc=[ -0.2 -0.1 -0.05 -0.04 -0.03 -0.02 -0.01 -0.005 0 0.005 0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 ];
lg=min(min(p.a(:,:,5))):(min(min(p.a(:,:,5))) - max(max(p.a(:,:,5))))/5:max(max(p.a(:,:,5)));
l=[lc lg];
[c,h] = contour (p.a(:,:,5)', l); colorbar;
%clabel(c,h);
axis ([1 maxindex 1 maxindex]);
title ('Pressure');

figure (5); 
surf (p.a(:,:,5)');shading interp;
axis ([1 maxindex 1 maxindex min(min(p.a(:,:,5)))*0.95 max(max(p.a(:,:,5)))*1.05]);
title ('Pressure');
