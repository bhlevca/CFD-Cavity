c...cubic spline approximation from num. recip.
c   to construct wall normal spacing for step configuration
c   y-positions are written on file
c
c  cubic spline approximation of vertical spacing
c  continuous in first and second derivative
c  stretching is delta(dy)/dy
c   

      program main

      implicit none
  
      integer MAXD, MAXP 
      parameter(MAXD=50, MAXP=500)

      real x(MAXD), dx(MAXD), odx(MAXD), y2(MAXD), ystore(MAXP), 
     &     stiffness(MAXP), dsold(MAXP),fako,fak,ds(MAXP), 
     &     xe(MAXD), ystoree(MAXP) ,xstart ,deltax, residual   

      real bdry, xloc, dxloc, height, stretch, dxloc_old       
      
      integer nmax, ne ,ke(MAXD) ,i, nout ,neval,l, k,ny ,j,n,pl 

      character*200 path
      character*200 file
      character     eos  


      print*,' Anzahl der DefinitionsPunkte furs Gitter '
      read*,nmax
      print*,' Eingabe der x-Positionen der Definit. Punkte' 
      do i=1,nmax
        read*,x(i) 
      end do      
      print*,' Eingabe der dx-Werte der Definit. Punkte' 
      do i=1,nmax
        read*,odx(i) 
      end do      

      print*,'Anz.d. x-Pos., d. exakt auf d. Gitter liegen sollen >=2'
      read*,ne
      print*,'Eingabe d. x-Pos. d. exakt auf d. Gitter liegen sollen'
      do i=1,ne
        read*,xe(i) 
      end do      
  
      if (ne.lt.2) then
        print *,' Fehler Anz. kleiner als 2 '
        goto 9999
      end if      
     
      print*,'Faktor fur Gitterpunktzahl (sonst 1.0):'
      read*,fako    
      fak=fako 

      do i=1,200
       path(i:i)=' '
      end do

      print*,' Pfad fur Ergebnisfiles '
      read *,path
c      print *,path    

      do i=1,200
       if (path(i:i).ne." ") pl=i 
      end do            
c      print *,'PL ',pl 
      
c      file='test.fil'
c      path(pl+1:pl+8)=file
c      print *,path 
      
c      open(10,file=path,form='formatted')
c      close(10)  
       
 

 1    do i=1,nmax
        dx(i)=fak*odx(i)
      end do

      bdry=2.e30
      call spline(x,dx,nmax,bdry, bdry, y2)
      print*,' *** spline computed'
      do i=1,nmax
        write(*,910) i,x(i),dx(i),y2(i)
      enddo
  910 format(1x,i8,3(1x,e12.5)) 

      nout=10
      open(nout,file='spline.dat',form='formatted')

      xstart=x(1)
      neval=15
      do l=2,nmax
        deltax = (x(l)-x(l-1))/(neval-1)
        print*,' l  deltax = ',l,deltax
        do i=1,neval
          xloc = xstart + (i-1)*deltax
          call splint(x, dx, y2, nmax, xloc, dxloc)
          write(nout,900) xloc, dxloc
        enddo
        xstart = xloc
      enddo
      close(nout)

c      open(nout,file='spline.loc',form='formatted')
c      do i=1, nmax
c        write(nout,900) x(i), dx(i)
c      enddo
c      close(nout)

c...count numer of points in y
      print*,' *** j   y(j)  dy(j) stretch'
      ny=1
      height=x(nmax)
cFK      xloc=0.
      dxloc=0.
      stretch=0.
cFK      ystore(1)=0.
       xloc=x(1)
       ystore(1)=x(1)

      write(*,920) ny, xloc, dxloc, stretch
   40 continue
      dxloc_old=dxloc
      call splint(x, dx, y2, nmax, xloc, dxloc)
      if(dxloc_old.gt.0.) stretch=dxloc/dxloc_old
      xloc=xloc+dxloc
      ny=ny+1
      ystore(ny)=xloc
c      write(*,920) ny, xloc, dxloc, stretch
      if(xloc.lt.height) goto 40

c
c  Genau einpassen der GitterPunkte aufs gewunschte Gitter 
c

c Global einpassen 
  
      do j=1,ny
        ystore(j)=ystore(1)+(ystore(j)-ystore(1))*(x(nmax)-x(1))/
     *                                      (ystore(ny)-ystore(1))
       print*,j,ystore(j)
      end do

c Suchen der event. zu korrigierenden Bereiche 

      do i=1,ne
       call findnearest(xe(i),ystore,ny,ke(i))
       print*,ke(i)
      end do

c Plausibilitatsprufung

      do i=1,ne-1
        if (ke(i).ge.ke(i+1)) then
          print *,'Punkte fur exakte Gitterlinien zu dicht 
     $ oder nicht geordnet angegeben'
          goto 9999 
        end if  
      end do       

c Verteilung der Gitterzellen

      do i=1,ne-1
       do l=ke(i),ke(i+1)
        ystoree(l)=xe(i)+(ystore(l)-ystore(ke(i)))*(xe(i+1)-xe(i))/
     $                               (ystore(ke(i+1))-ystore(ke(i)))
       end do
      end do   

c
c Bem. zu obigem Alg.: Setze   h=ystore(ke(i+1))  -ystore(ke(i+1)-1) 
c                              H=ystore(ke(i+1)+1)-ystore(ke(i+1)
c
c                              he=ystoree(ke(i+1))  -ystoree(ke(i+1)-1)
c                              He=ystoree(ke(i+1)+1)-ystoree(ke(i+1)
c Dann gilt:
c
c he= h * (xe(i+1)-xe(i)  ) / (ystore(ke(i+1))-ystore(ke(i))) 
c He= H * (xe(i+2)-xe(i+1)) / (ystore(ke(i+2))-ystore(ke(i+1))) 
c
c Nach Vorausseztung ist h=H(1+O(H))
c weiterhin ist  (xe(i+1)-xe(i)  ) / (ystore(ke(i+1))-ystore(ke(i)))   = 1+O(h)
c und analog     (xe(i+2)-xe(i+1)) / (ystore(ke(i+2))-ystore(ke(i+1))) = 1+O(H)
c
c Daraus folgt, dass he=H(1+O(H))*(1+O(h)) = H(1+O(H)) bzw. He(1+O(He)) 
c        
c Die Glattheitseigenschaften sollten also wieder einigermassen erfullt sein.
c

 999  print*,' resdiual distributed= ',residual

      file='gridpoints.dat'
      path(pl+1:pl+14)=file
      open(10,file=path,form='formatted')
      write(10,*) ny
      do j=1,ny
        write(10,940) j,ystoree(j)
      enddo
      close(10)  

      file='gridpoints.gnu'
      path(pl+1:pl+14)=file
      open(10,file=path,form='formatted')
      do j=1,ny
        write(10,940) j,ystoree(j)
      enddo
      close(10)

      file='gridspaces.gnu'
      path(pl+1:pl+14)=file
      open(10,file=path,form='formatted')
      do j=1,ny-1
        write(10,*) ystore(j),ystoree(j+1)-ystoree(j)
      enddo
      close(10)

      print*,' *** y-locations with numbers on file step.ygrid: ny= ',ny

c      open(10,file='grid.dat',form='formatted')
c      write(10,*) ny
c      do j=1,ny
c        write(10,*) ystore(j)
c      enddo
c      close(10)

      print*,' *** y-locations without numbers on file grid.dat: ny= ',ny
      print*,' *** dy(y) fur gnuplot in spline.dat'
     

  920 format(1x,i4,3(1x,f12.5)) 
  940 format(1x,i4,3(f14.7,1x))
  950 format(1x,3(f14.7,1x))

  900 format(1x,3(e12.5,1x))

 9999 end
c.......................................................................
  
      subroutine findnearest(x,y,n,k)

      implicit none 
        
       real    mx,x,y(1:n)
       integer i,k,n  
    
       mx=10000000000 
 
       do i=1,n 
        if (abs(x-y(i)).lt.mx) then
         mx=abs(x-y(i)) 
         k=i   
        end if
       end do
    
      end 
        


      subroutine spline(x,y,n,yp1,ypn,y2)
c
c  input: x(1:n), y(1:n)     tabulated function
c         yp1, ypn           first deriv. at endpoints
c                            if > 1.e30 : natural spline, zero second deriv
c  returns:   y2(1:n)        second deriv. of interpolating at tabulated points
c         
      real x(*), y(*),y2(*)
      integer i,k
cgb      real p,qn,sig,un,u(n)
      real p,qn,sig,un,u(600)

      if(yp1.gt.0.99e30) then
        y2(1)=0.
        u(1) =0.
      else
        y2(1) = -0.5
        u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1)+2.
        y2(i) = (sig-1.)/p
        u(i)=(y(i+1)-y(i))/(x(i+1)-x(i)) 
     &     - (y(i)-y(i-1))/(x(i)-x(i-1))
        u(i)=(6.*u(i)/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo

      if(ypn.gt. 0.99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1) + u(k)
      enddo
      return
      end

c............................................................................
      subroutine splint(xa, ya, y2a, n, x, y)
c
c    input: xa, ya   tabulated function
c           y2a      output from spline
c           x        evaluate spline here
c           y        interpolated value

      real xa(*),ya(*),y2a(*)

      integer klo, khi, k
      real h,b,a

      klo=1
      khi=n

   10 continue
        k=(khi+klo)/2
        if(xa(k).gt.x) then 
          khi=k
        else
          klo=k
        endif
      if( (khi-klo) .gt. 1) goto 10

      h=xa(khi)-xa(klo)
      if (h.eq.0.) print*,' *** bad xa input to splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo) + b*ya(khi) + ((a*a*a-a)*y2a(klo)
     & +(b*b*b-b)*y2a(khi))*(h*h)/6.0
      return
      end
c ---------------------------------------------------------------------
      subroutine springs(stiff, residual , delta,  n)

c purpose: distribute residual on n segments with given stiffness
c          by solving a tridiagonal system
c
c input:  stiffness for n segments (should be non-zero)
c         residual:  amount which will be distributed 
c output: delta contains individual change for each of
c         the n segments; individual changes sum up to residual

      parameter(nmax=300)
      real stiff(*),residual,delta(*)
      real a(nmax),b(nmax),c(nmax),r(nmax),shifts(nmax)

      if(stiff(n).gt.1.e10) then
        print*,' attention: stiffness of last element to high'
        stop
      endif
      do i=2,n-2
        r(i)=0.
        a(i)=-stiff(i)
        b(i)= stiff(i)+stiff(i+1)
        c(i)=-stiff(i+1)
      enddo
      b(1)=stiff(1)+stiff(2)
      c(1)=-stiff(2)
      r(1)=0.
      a(n-1)=-stiff(n-1)
      b(n-1)=stiff(n-1)+stiff(n)
      r(n-1)=residual * stiff(n)
      call tridag(a,b,c,r,shifts,n-1)
      delta(1)=shifts(1)
      shifts(n)=residual
      do i=2,n
        delta(i)=shifts(i)-shifts(i-1)
      enddo
      sum=0.
c      do i=1,n
c        sum=sum+delta(i)
c      enddo
c      print*,' spring: residual, sum segm= ',residual,sum
      return
      end
    
c---------------------------------------------------------------------
      SUBROUTINE tridag(A,B,C,R,U,n)

c  purpose: tridiag. matrix solver from num. recipes

      parameter(nmax=2000)

      real gam(nmax),A(*),B(*),C(*),R(*),U(*)

      if(b(1).eq.0.) then
        print*,' tridag: rewrite matrix'
        stop
      endif
      if(n.gt.nmax) then
        print*,' tridag: enlarge nmax'
        stop
      endif
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
          print*,' error in tridag'
          stop
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      do j=n-1,1,-1
        u(j) =u(j) -gam(j+1)*u(j+1)
      enddo
      return
      end

