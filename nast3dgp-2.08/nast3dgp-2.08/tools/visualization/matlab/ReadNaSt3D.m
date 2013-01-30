function U = ReadExplorer(fname) 
%  U = ReadExplorer(fname)  
%  reads rectilinear data from a file in MatLab/Explorer format 
%  U is a struct, which contains:
%
%  U.n(3)       ; dimensions
%  U.x(1..n(1)) ; x coordinates of grid lines 
%  U.y(1..n(2)) ; y coordinates of grid lines 
%  U.z(1..n(3)) ; z coordinates of grid lines
%  U.a(1..n(1) , 1..n(2) , 1..n(3)) ; data
%

[fid,msg]=fopen(fname,'r') ; 

if fid < 0 
  msg 
  fname
end 

% read dimensions

U.n=fread(fid,3,'int') ;

% read Grid

U.x=fread(fid,U.n(1),'float') ;
U.y=fread(fid,U.n(2),'float') ;
U.z=fread(fid,U.n(3),'float') ;

% read Data 

U.a=zeros(U.n(1),U.n(2),U.n(3)) ;

for k=1:U.n(3)
 s=fread(fid,[U.n(1),U.n(2)],'float') ;
 U.a(:,:,k)=s ;
end 

fclose(fid) ;

