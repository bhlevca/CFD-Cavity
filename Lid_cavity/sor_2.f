C****************************************************************************C
C                                                                            C
C     This program solve the Navier-Stokes (NS) equations using the          C
C     Successive Over Relaxation (SOR) method. The solution is spatially     C
C     second order accurate. This code is written by Prof. Ercan Erturk.     C
C     Visit http://www.cavityflow.com                                        C
C                                                                            C
C********************************************C*******************************C
C                                            C
C     s(i,j) ==> streamfunction variable     C
C     v(i,j) ==> vorticity variable          C
C     x(i)   ==> x-coordinate                C
C     y(j)   ==> y-coordinate                C
C     dh     ==> grid spacing                C
C     Re     ==> Reynolds Number             C
C                                            C
C********************************************C

      program main

      implicit double precision (a-h,o-z)

      parameter(N=128)

      common / flow variables /
     > s(0:N,0:N),v(0:N,0:N),
     > s_old(0:N,0:N),v_old(0:N,0:N)
      common / geometry /
     > x(0:N),y(0:N)


       Re=1000.d0

       dh=1.0d0/dble(N)

      do 1 k=0,N
       x(k)=dble(k)/dble(N)
       y(k)=dble(k)/dble(N)
 1    continue

C     Relaxation parameter
       beta=0.6d0

C     Initial guess
c     Note that when homogeneous initial guess is used, in the first iteration 
c     residual_3 is indeterminate. In order to avoid this, instead of using 
c     exactly zero initial values, at interior points use very very small 
c     numbers that could be considered as zero.
      do 2 k=0,N
      do 22 j=0,N
       s(k,j)=1.d-32
       v(k,j)=1.d-32
 22    continue
       s(0,k)=0.d0
       v(0,k)=0.d0
       s(N,k)=0.d0
       v(N,k)=0.d0
       s(k,0)=0.d0
       v(k,0)=0.d0
       s(k,N)=0.d0
       v(k,N)=0.d0
 2    continue

C     Record the CPU time at start
       call cpu_time(time_start)

C     Start the iterations
      do 999 iteration=1,1000000

C     Update old variables
      do 3 i=1,N-1
      do 3 j=1,N-1
       s_old(i,j)=s(i,j)
       v_old(i,j)=v(i,j)
 3    continue

C     Update Streamfunction variable
      do 4 i=1,N-1
      do 4 j=1,N-1
       s(i,j)=beta*(   (
     >        s(i-1,j)+s(i+1,j)+s(i,j-1)+s(i,j+1)
     >       +dh**2.*v(i,j)
     >                 )/4.d0
     >             )+(1.d0-beta)*s(i,j)
 4    continue
      
C     Calculate Vorticity at the wall
C     NOTE:For these boundary conditions please refer to:
C          T. Stortkuhl, C. Zenger, S. Zimmer, "An Asymptotic Solution for
C          the Singularity at the Angular Point of the Lid Driven Cavity",
C          International Journal of Numerical Methods for Heat & Fluid Flow
C          1994, Vol 4, pp 47--59
      do 5 k=1,N-1
       v(k,0)=(
     >       -(s(k-1,1)+s(k,1)+s(k+1,1))/(3.d0*dh**2.)
     >       -(0.5d0*v(k-1,0)+0.5d0*v(k+1,0)
     >       +0.25d0*v(k-1,1)+v(k,1)+0.25d0*v(k+1,1))/(9.d0)
     >        )*9.d0/2.d0
       v(k,N)=(-1.d0/dh
     >       -(s(k-1,N-1)+s(k,N-1)+s(k+1,N-1))/(3.d0*dh**2.)
     >       -(0.5d0*v(k-1,N)+0.5d0*v(k+1,N)
     >       +0.25d0*v(k-1,N-1)+v(k,N-1)+0.25d0*v(k+1,N-1))/(9.d0)
     >        )*9.d0/2.d0
       v(0,k)=(
     >       -(s(1,k-1)+s(1,k)+s(1,k+1))/(3.d0*dh**2.)
     >       -(0.5d0*v(0,k-1)+0.5d0*v(0,k+1)
     >       +0.25d0*v(1,k-1)+v(1,k)+0.25d0*v(1,k+1))/(9.d0)
     >        )*9.d0/2.d0
       v(N,k)=(
     >       -(s(N-1,k-1)+s(N-1,k)+s(N-1,k+1))/(3.d0*dh**2.)
     >       -(0.5d0*v(N,k-1)+0.5d0*v(N,k+1)
     >       +0.25d0*v(N-1,k-1)+v(N-1,k)+0.25d0*v(N-1,k+1))/(9.d0)
     >        )*9.d0/2.d0
 5    continue
       v(0,0)=(-(s(1,1))/(3.d0*dh**2.)-(0.5d0*v(1,0)+0.5d0*v(0,1)
     >       +0.25d0*v(1,1))/(9.d0))*9.d0
       v(N,0)=(-(s(N-1,1))/(3.d0*dh**2.)-(0.5d0*v(N-1,0)+0.5d0*v(N,1)
     >       +0.25d0*v(N-1,1))/(9.d0))*9.d0
       v(N,N)=(-0.5d0/dh-(s(N-1,N-1))/(3.d0*dh**2.)-(0.5d0*v(N-1,N)
     >       +0.5d0*v(N,N-1)+0.25d0*v(N-1,N-1))/(9.d0))*9.d0
       v(0,N)=(-0.5d0/dh-(s(1,N-1))/(3.d0*dh**2.)-(0.5d0*v(1,N)
     >       +0.5d0*v(0,N-1)+0.25d0*v(1,N-1))/(9.d0))*9.d0

C     Update Vorticity variable
      do 6 i=1,N-1
      do 6 j=1,N-1
       v(i,j)=beta*(  (
     >        v(i-1,j)+v(i+1,j)+v(i,j-1)+v(i,j+1)
     >       -0.25d0*Re*(s(i,j+1)-s(i,j-1))*(v(i+1,j)-v(i-1,j))
     >       +0.25d0*Re*(s(i+1,j)-s(i-1,j))*(v(i,j+1)-v(i,j-1))
     >                )/4.d0
     >             )+(1.d0-beta)*v(i,j)
 6    continue


C     Check the residuals to see if convergence is achieved. 
C     You can comment out all or some, for a faster run.

c     residual_1 is the residual of the governing equations.
       residual_1_s_A=0.d0
       residual_1_v_A=0.d0
c     residual_2 is the change in variables (indicates the significant digit) in a time step.
       residual_2_s_A=0.d0
       residual_2_v_A=0.d0
c     residual_3 is the normalized change in variables (indicates percent change) in a time step.
       residual_3_s_A=0.d0
       residual_3_v_A=0.d0
      
      do 60 i=1,N-1
      do 60 j=1,N-1
       residual_1_s_B=abs(
     >   (s(i-1,j)-2.d0*s(i,j)+s(i+1,j))/dh**2.
     >  +(s(i,j-1)-2.d0*s(i,j)+s(i,j+1))/dh**2.
     >  +v(i,j)          )
       residual_1_v_B=abs(
     >   (1.d0/Re)*(v(i-1,j)-2.d0*v(i,j)+v(i+1,j))/dh**2.
     >  +(1.d0/Re)*(v(i,j-1)-2.d0*v(i,j)+v(i,j+1))/dh**2.
     >  -(s(i,j+1)-s(i,j-1))/(2.d0*dh)*(v(i+1,j)-v(i-1,j))/(2.d0*dh)
     >  +(s(i+1,j)-s(i-1,j))/(2.d0*dh)*(v(i,j+1)-v(i,j-1))/(2.d0*dh)
     >                   )

       residual_2_s_B=abs(s(i,j)-s_old(i,j))
       residual_2_v_B=abs(v(i,j)-v_old(i,j))

       residual_3_s_B=abs((s(i,j)-s_old(i,j))/s_old(i,j))
       residual_3_v_B=abs((v(i,j)-v_old(i,j))/v_old(i,j))

       residual_1_s_A=max(residual_1_s_A,residual_1_s_B)
       residual_1_v_A=max(residual_1_v_A,residual_1_v_B)

       residual_2_s_A=max(residual_2_s_A,residual_2_s_B)
       residual_2_v_A=max(residual_2_v_A,residual_2_v_B)

       residual_3_s_A=max(residual_3_s_A,residual_3_s_B)
       residual_3_v_A=max(residual_3_v_A,residual_3_v_B)
 60   continue
      
C     Output the residuals
       write(*,*) ' '
       write(*,*) iteration
       write(*,*) residual_1_s_A,residual_1_v_A
       write(*,*) residual_2_s_A,residual_2_v_A
       write(*,*) residual_3_s_A,residual_3_v_A
c     condition to stop iterations
       if((residual_1_s_A.lt.1.d-6).AND.(residual_1_v_A.lt.1.d-6))
     > goto 1000
      
 999  continue
      
C     Record the CPU time at finish
 1000 call cpu_time(time_finish)

       write(*,*) ' '
       write(*,*) 'Convergence is achieved in',iteration,'   iterations'
       write(*,*) 'CPU time=',time_finish-time_start

c     Output to a file        
      open(1,file='out.txt')
      do 1111 i=0,N
      do 1111 j=0,N
       write(1,2222) x(i),y(j),s(i,j),v(i,j)
 1111 continue
 2222 format(f,x,f,x,es,x,es)

      stop
      end
