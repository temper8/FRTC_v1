      subroutine drivencurrent(outj,sigmaj)
cc******************************************************************
cc   outj(i)  = LH driven current density, MA/m^2
cc   dndt(i)  = d^2Jr1/dt^2/E, MA/m^2/sec^2/(V/m), ~runaway d(el.density)/dt/E
cc   djdt(i)  = dJr2/dt, time drivative of runaway current Jr2, MA/m^2/sec
cc   outjrun(i)  = LH driven runaway current density, MA/m^2
cc   outnerun(i) = runaway electron density/10^19 m^-3
cc******************************************************************
      implicit none
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      real*8 outj(NRD),sigmaj(NRD),afld(NRD),dtau
      integer i,inpt,ispectr,ntau
      real*8,dimension(:),allocatable:: outjp,outjm,ohjp,ohjm
      real*8 dt,zero,eps,cup,cup0,cum,cum0,cp,cm,cp0,cm0,aiint
      parameter(zero=0.d0,eps=1.d-2,ntau=10)
!
      inpt=NA1
      allocate(outjp(inpt),outjm(inpt),ohjp(inpt),ohjm(inpt))
      do i=1,inpt
       afld(i)=UPL(i)/RTOR/GP2 !!variant
      end do
!
!!!!!!!!!!!!! starting LH current calculation !!!!!!!!!!!!!!!!!
      outj=zero
      outjp=zero
      outjm=zero
      ohjp=zero
      ohjm=zero
      cup=zero
      cum=zero
      cp=zero
      cm=zero
      cup0=zero
      cum0=zero
      cp0=zero
      cm0=zero
!
!!positive spectrum:
       ispectr=1
       call lhcurrent(outjp,ohjp,cup,cup0,inpt,ispectr)
       if(cup0.ne.zero) then
        cp0=aiint(ohjp,roc)
        if(cp0.ne.zero) then
         do i=1,inpt
          ohjp(i)=cup0*ohjp(i)/cp0
         end do
        end if
       end if
       if(cup.ne.zero) then
        cp=aiint(outjp,roc)
        if(cp.ne.zero) then
         do i=1,inpt
          outjp(i)=cup*outjp(i)/cp
         end do
        end if
       end if

!!negative spectrum:
       ispectr=-1
       call lhcurrent(outjm,ohjm,cum,cum0,inpt,ispectr)
       if(cum0.ne.zero) then
        cm0=aiint(ohjm,roc)
        if(cm0.ne.zero) then
         do i=1,inpt
          ohjm(i)=cum0*ohjm(i)/cm0
         end do
        end if
       end if
       if(cum.ne.zero) then
        cm=aiint(outjm,roc)
        if(cm.ne.zero) then
         do i=1,inpt
          outjm(i)=cum*outjm(i)/cm
         end do
        end if
       end if

      do i=1,inpt
       outj(i)=outjp(i)+outjm(i)
       sigmaj(i)=zero
       if(dabs(afld(i)).gt.eps) then
        sigmaj(i)=(ohjp(i)+ohjm(i))/afld(i)
       end if
!!!!       write(*,*) i,outj(i)
      end do
!
      dt=tau/dble(ntau) !seconds 
      write(*,*)'time=',time,' dt=',dt
      write(*,*)'cup=',cup,' cp=',cp
      write(*,*)'cum=',cum,' cm=',cm
      write(*,*)'cup0=',cup0,' cp0=',cp0
      write(*,*)'cum0=',cum0,' cm0=',cm0
      write(*,*)'sigma driven current, MA=',cp0+cm0
      write(*,*)'driven current, MA=',cup+cum
      write(*,*)
!
!! calculation of distribution functions at time t1=t+dtau !!
      do i=1,ntau
       write(*,*)'fokkerplanck №',i,'of',ntau
       call fokkerplanck(dt,i,tau,time)
      end do
!
      deallocate(outjp,outjm,ohjp,ohjm)
!
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine lhcurrent(outj,ohj,cuj,cujoh,inpt,ispectr)
!!      implicit real*8 (a-h,o-z)
      implicit none
      real*8 outj(*),ohj(*),cuj,cujoh,curs,curs0,curdir
      real*8 currn,pqe,vt0,fvt,ccur,cfull,cfull0
      real*8 r,pn,fn1,fn2,fnr,fnrr,vt,vto,rh1
      integer nr,klo,khi,ierr,nrr,i,j,inpt,ispectr,ismthout
      common /a0ab/ nr
      real*8 rh,y2dn,y2tm,y2tmi
      common /a0l3/ rh(501),y2dn(501),y2tm(501),y2tmi(501)
      integer inew
      common /cnew/ inew !est !sav2008
      real*8 zv1,zv2,sk,fout
      common/plosh/ zv1(100,2),zv2(100,2),sk(100)
      integer i0,k
      parameter(i0=1002)
      real*8 vij,fij0,fij,dfij,dij,enorm,fst
      common/lh/vij(i0,100),fij0(i0,100,2),fij(i0,100,2),dfij(i0,100,2)
     &,dij(i0,100,2),enorm(100),fst(100)
      real*8,dimension(:),allocatable:: vj,fj,fj0,cur,cur0
     &,currnt,rxx,wrk
      real*8 zero
      parameter(zero=0.d0, ismthout=1)
!
      allocate(vj(i0),fj(i0),fj0(i0),cur(nr),cur0(nr)
     &         ,currnt(nr+2),rxx(nr+2),wrk(nr+2))
c---------------------------------------------------
c initial constants
c---------------------------------------------------
      pqe=4.803e-10
      vt0=fvt(zero)
      ccur=pqe*vt0*0.333d-9
      curdir=-dble(ispectr)
!
      cfull=zero
      cfull0=zero
      k=(3-ispectr)/2
      do j=1,nr
       do i=1,i0
        vj(i)=vij(i,j) !Vpar/Vt
        fj0(i)=fij0(i,j,k)
        fj(i)=fij(i,j,k)-fij0(i,j,k)
!!        fj(i)=zero
!!        if((vj(i)-zv1(j,k))*(vj(i)-zv2(j,k)).le.zero) then
!!         fj(i)=fij(i,j,k)-fij0(i,j,k)
!!        end if
       end do
       r=dble(j)/dble(nr+1)
       if(inew.eq.0) then !vardens
        pn=fn1(r,fnr)
       else
        pn=fn2(r,fnr,fnrr)
       end if
       vt=fvt(r)
       vto=vt/vt0
       call currlhcd(i0,vj,fj,fj0,curs,curs0)
       cur(j)=curs*pn*ccur*curdir*vto  !Ampere/cm2
       cfull=cfull+cur(j)*sk(j)
       cur0(j)=curs0*pn*ccur*curdir*vto  !Ampere/cm2
       cfull0=cfull0+cur0(j)*sk(j)
!!!       tok(j)=cur(j)*sk(j) !Ampere
!!!       write(*,88) dble(j),cur(j)*sk(j)
      end do
      cuj=cfull*1d-6   !driven current, MA
      cujoh=cfull0*1d-6   !driven current, MA
!!      write(*,*)
!!      write(*,*)'ccur',ccur,' curdir=',curdir,' nr=',nr
!!      write(*,*)'cu_out, MA=',cu_out,' cfull, A=',cfull
!!           close(111)
c      pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      currn=cur(1)                   ! Jstoped, A/cm^2
      currnt(1)=currn*1.d-2          ! Jstoped, MA/m^2
      rxx(1)=zero
      do j=1,nr
       rxx(j+1)=dble(j)/dble(nr+1)
       currn=cur(j)                   ! Jstopped, A/cm^2
       currnt(j+1)=currn*1.d-2        ! Jstoped, MA/m^2
      end do
      nrr=nr+2
      rxx(nrr)=1.d0
      currnt(nr+2)=zero
!
      if(ismthout.ne.0) then
       do i=1,nrr
        wrk(i)=currnt(i)
       end do
       call fsmoth4(rxx,wrk,nrr,currnt)
      end if
!
      rh(1)=rh1
      if(rh(inpt).gt.1d0) rh(inpt)=1.d0
      do j=1,inpt
       call lock2(rxx,nrr,rh(j),klo,khi,ierr)
       if(ierr.ne.0) then
        write(*,*)'lock2 error in current profile for ASTRA'
        write(*,*)'ierr=',ierr,' j=',j,' rh(j)=',rh(j)
        write(*,*)'rxx(1)=',rxx(1),' rxx(nrr)=',rxx(nrr)
        pause
       end if
       call linf(rxx,currnt,rh(j),fout,klo,khi)
       outj(j)=fout
      end do
      rh(1)=zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      currn=cur0(1)                   ! Jstoped, A/cm^2
      currnt(1)=currn*1.d-2          ! Jstoped, MA/m^2
      rxx(1)=zero
      do j=1,nr
       rxx(j+1)=dble(j)/dble(nr+1)
       currn=cur0(j)                   ! Jstopped, A/cm^2
       currnt(j+1)=currn*1.d-2        ! Jstoped, MA/m^2
      end do
      nrr=nr+2
      rxx(nrr)=1.d0
      currnt(nr+2)=zero
!
      if(ismthout.ne.0) then
       do i=1,nrr
        wrk(i)=currnt(i)
       end do
       call fsmoth4(rxx,wrk,nrr,currnt)
      end if
!
      rh(1)=rh1
      if(rh(inpt).gt.1d0) rh(inpt)=1.d0
      do j=1,inpt
       call lock2(rxx,nrr,rh(j),klo,khi,ierr)
       if(ierr.ne.0) then
        write(*,*)'#2 lock2 error in current profile for ASTRA'
        write(*,*)'ierr=',ierr,' j=',j,' rh(j)=',rh(j)
        write(*,*)'rxx(1)=',rxx(1),' rxx(nrr)=',rxx(nrr)
        pause
       end if
       call linf(rxx,currnt,rh(j),fout,klo,khi)
       ohj(j)=fout
      end do
      rh(1)=zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      deallocate(vj,fj,fj0,cur,cur0,currnt,rxx,wrk)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine currlhcd(i0,v,f,f0,curs,curs0)
      implicit none
      integer i0,k
      real*8 v(*),f(*),f0(*),curs,curs0
      real*8 vl,vr,fl,fr,zero
      parameter(zero=0.d0)
      curs=zero
      curs0=zero
      do k=1,i0-1
       vl=v(k)
       vr=v(k+1)
       fl=f(k)
       fr=f(k+1)
       curs=curs+(fl*vl+fr*vr)/2d0*(vr-vl)
       fl=f0(k)
       fr=f0(k+1)
       curs0=curs0+(fl*vl+fr*vr)/2d0*(vr-vl)
      end do
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fokkerplanck(dtstep,pachka,tay,tim1)
      implicit none
      real*8 t,dtstep,dtau
      integer nr,ispl,ip,im,ii,ibeg, iunit, iunit1, step, tc, koltoch
      common /a0ab/ nr
      real*8 ynzm0,pm0,plaunp,plaunm,fmaxw,d,dDdX,part,part1
      common/grillspektr/ ynzm0(1001),pm0(1001),ispl
     &,plaunp,plaunm,ip,im
      integer i0,iptnew
      parameter(i0=1002)
      real*8 vij,fij0,fij,dfij,dij,enorm,fst,curtime
      common/lh/vij(i0,100),fij0(i0,100,2),fij(i0,100,2),dfij(i0,100,2)
     &,dij(i0,100,2),enorm(100),fst(100)
      real*8 dijk,vrjnew
      common/t01/dijk(101,100,2),vrjnew(101,100,2),iptnew
      integer n,i,j,it,nt,k,pachka
      real*8 xend,h,shift,ybeg,yend,tend,dt,dff,tcur,tay,tim1
      real*8,dimension(:),allocatable:: y,x,xx,f,ddD1,D2,y1
      real*8,dimension(:),allocatable:: vj,fj,dfj,ft,D1,dD1,dD2
      real*8,dimension(:),allocatable:: a,b,c,a0,b0,c0
      real*8 znak,alfa2,zero,dt0,h0,eps,r,fvt,fout1,fout2
      common/ef/ alfa2
      real*8 calls
      common/firstcall/calls
      real*8 d0,B1,C1,w
      real*8 score1
      real*16 tau1,spacing,curtime0,h1
      double precision :: minmod,minmod1
      external B1,C1,w
      integer jindex,kindex,klo,khi,ierr,i2,i1,klo1,khi1
      common/dddql/ d0,jindex,kindex
      parameter(tau1=3.000990745207882E-002)
      parameter(zero=0.d0,dt0=1.d-1,h0=1.d-1,eps=1.d-7)
!      real*8 tcur
      common/testf/ tcur
!       
      do k=1,2
       kindex=k
       znak=2.d0*dble(k)-3.d0 !! k=1 znak=-1, in abc Efield with "+"
!
!!       do j=1,nr
!!        fij0(i,j,k)=fmaxw(vij(i,j),znak*enorm(j),dff)
!!       end do
!!       if(1.gt.0) go to 1
!
       d0=zero
       do j=1,nr
       
c       do i=1,iptnew
c       write(*,*) dijk(i,j,k), vrjnew(i,j,k)
c       enddo
c       pause
c       do i=1, iptnew
c       !if (iptnew.ne.0) then
c       write(*,*) dijk(i,j,1)
c       !endif
c       enddo
c       pause
        jindex=j
        dtau=dtstep*fst(j)
        nt=1
        if(dtau.gt.dt0) then
         nt=1+dtau/dt0
        end if
        dt=dtau/nt
	
        r=dble(j)/dble(nr+1)
        xend=3.d10/fvt(r)
!!        xend=vij(i0,j)
        n=xend/h0-1
        h=xend/dble(n+1)
        if(h.gt.h0) then
         n=n+1
         h=xend/dble(n+1)
        end if
        allocate(y(n+2),x(n+2),xx(n+1),a(n+1),b(n+1),c(n+1),f(n+1)
     &,y1(n+1))
!!!!!! grid !!!!!!!!!
!!       shift=h*0.1d0 !0.01d0
        do i=1,n+2
         x(i)=h*dble(i-1) !+shift
        end do
!	write(*,*)'xend=',xend,'n=',n,'h=',h
        do i=1,n+1
         xx(i)=h/2.d0+h*dble(i-1) !+shift
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(D1(n+1))

        D1=0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(vj(i0),fj(i0))
        
        do i=1,i0
         vj(i)=vij(i,j)
         fj(i)=fij0(i,j,k)
        end do
        
        
        do i=1,n+2
         call lock(vj,i0,x(i),klo,khi,ierr)
         if(ierr.eq.1) then
          write(*,*)'lock error #1 in function fokkerplanck'
          write(*,*)'j=',j,' v=',x(i+1)
          write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
          pause
          stop
         end if
         call linf(vj,fj,x(i),y(i),klo,khi)
        end do
        deallocate(fj)
        ybeg=fij0(1,j,k)  !boundary conditions
        yend=fij0(i0,j,k)
        alfa2=znak*enorm(j)
!!!solve problem
        
        do i=1,n+1
            y1(i)=y(i+1)
        end do
        
        do it=1,nt
            call abccoef(a,b,c,f,y1,dt,n+1,ybeg,yend,x,xx,h,D1)
            call tridag(a,b,c,f,y1,n+1)
            
        
        do i=1,n+1
        if(y1(i).lt.0.d0) then
        if((y1(i) + epsilon(y1(i))).gt.0.d0) then
        y1(i)=0.d0
        else
        write(*,*) 'y(i)=',y1(i),' lt negative epsilon=',epsilon(y1(i))
        pause
        stop
        endif
        endif
        enddo
        
        end do

      do i=1,n+1
        y(i+1)=y1(i)
      end do

        do i=1,i0
         if(vij(i,j).lt.xend) then
          call lock(x,n+2,vij(i,j),klo,khi,ierr)
          if(ierr.eq.1) then
           write(*,*)'lock error #2 in finction fokkerplanck'
           write(*,*)'j=',j,' vij=',vij(i,j)
           write(*,*)'x(1)=',x(1),' x(n+2)=',x(n+2)
           pause
           stop
           end if
          call linf(x,y,vij(i,j),fij0(i,j,k),klo,khi)
         else
          fij0(i,j,k)=zero
         end if
        end do

c        call integral(1,i0,vij(:,j),fij0(:,j,k),fout1)
c        if(k.eq.1) then
c            open(8,file='RESULT2019/part01.dat',position="append")
c            write(8,*) tim1, fout1
c            close(8)
c        elseif (k.eq.2) then
c            open(9,file='RESULT2019/part02.dat',position="append")
c            write(9,*) tim1, fout1
c            close(9)
c        endif

        deallocate(y,x,xx,a,b,c,f,D1,y1)

        allocate(fj(i0),dfj(i0))
        do i=1,i0
         fj(i)=fij0(i,j,k)
         dfj(i)=zero
        end do


        do i=1,i0
         if(i.eq.1) then
          dfj(i)=zero
         else if(i.eq.i0) then
          dfj(i)=(fj(i)-fj(i-1))/vij(2,j)
         else
          dfj(i)=0.5d0*(fj(i+1)-fj(i-1))/vij(2,j)
         end if
        end do
	
        ii=0
        ibeg=0

c        do i=i0-1,1,-1
c         if(dfj(i).gt.zero) then
cc          write(*,*) '#1 positive derivs'
cc          write(*,*) '#1 df>0: i,j,k=',i,j,k
cc          write(*,*) '#1 dfj(i),i,j,k=',dfj(i),i,j,k
cc          write(*,*)
c          fij0(i,j,k)=fij0(i+1,j,k)
c          ii=i
c         end if
c         if(fij0(i,j,k).lt.fij0(i+1,j,k)) then 
c          fij0(i,j,k)=fij0(i+1,j,k)
c          ii=i
c         end if
c       end do 
!     
        ibeg=ii

c        if(ibeg.gt.0) then
c         call integral(ibeg,i0,vj,fj,fout1)
c         do i=1,i0
c         fj(i)=fij0(i,j,k)
c         end do
c         call integral(ibeg,i0,vj,fj,fout2)
c         do i=ibeg,i0
c          fij0(i,j,k)=fj(i)*fout1/fout2
c         end do
!!         write(*,*)'#1 j,k,ibeg=',j,k,ibeg
!!         write(*,*)'#1 v(ibeg)=',vj(ibeg),' f1/f2=',fout1/fout2
c        end if

!!!!	write(*,*)'time=',tau1


        deallocate(vj,fj,dfj)
       end do
!!!!!!!!!!!!!!!!!
!
1      continue
!
       d0=1.d0
       do j=1,nr
        jindex=j
        dtau=dtstep*fst(j)
        nt=1
        if(dtau.gt.dt0) then
         nt=1+dtau/dt0
        end if
        dt=dtau/nt
        r=dble(j)/dble(nr+1)
        xend=3.d10/fvt(r)
!!        xend=vij(i0,j)
        n=xend/h0-1
        h=xend/dble(n+1)
        if(h.gt.h0) then
         n=n+1
         h=xend/dble(n+1)
        end if
        allocate(y(n+2),x(n+2),xx(n+1),a(n+1),b(n+1),c(n+1),f(n+1)
     &,y1(n+1))
!!!!!! grid !!!!!!!!!
!!       shift=h*0.1d0 !0.01d0
        do i=1,n+2
         x(i)=h*dble(i-1) !+shift
        end do
        do i=1,n+1
         xx(i)=h/2.d0+h*dble(i-1) !+shift
        enddo

	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(D1(n+1),D2(101))

        D2=dijk(:,j,k)

        do i=1,n+1
        if(xx(i).gt.vrjnew(iptnew,j,k)) then
         D1(i)=zero
        else
        call lock(vrjnew(:,j,k),iptnew,xx(i),klo,khi,ierr)
        if(ierr.eq.1) then
         write(*,*)'lock error in finction d(x)'
         write(*,*)'j=',j,' v=',xx(i)
         write(*,*)'vj(1)=',vrjnew(1,j,k),' vj(i0)='
     & ,vrjnew(iptnew,j,k)
         pause
         stop
        end if
        call polint(vrjnew(klo,j,k),dijk(klo,j,k),2,xx(i),D1(i),ierr)
        
c        D1(i)=dij(klo,j,k)
        
        if(D1(i).lt.0d0) then
            write(*,*)'diff is  negative; check drivencurrent 558'
            write(*,*) 'xx(i)=', xx(i), 'D1(i)=', D1(i)
            pause
            stop
        endif
        
        endif
        enddo

c1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(vj(i0),fj(i0))
        do i=1,i0
         vj(i)=vij(i,j)
         fj(i)=fij(i,j,k)
        end do




        do i=1,n+2
         call lock(vj,i0,x(i),klo,khi,ierr)
         if(ierr.eq.1) then
          write(*,*)'lock error #3 in finction fokkerplanck'
          write(*,*)'j=',j,' v=',x(i+1)
          write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
          pause
          stop
         end if
         call linf(vj,fj,x(i),y(i),klo,khi)
        end do
        
        deallocate(fj)


        ybeg=fij(1,j,k)  !boundary conditions
        yend=fij(i0,j,k)
        alfa2=znak*enorm(j)
!!!!!!!         START ABCCOEF SOLVING           !!!!!!!

        do i=1,n+1
            y1(i)=y(i+1)
        end do
        
        do it=1,nt

         call abccoef(a,b,c,f,y1,dt,n+1,ybeg,yend,x,xx,h,D1)
         call tridag(a,b,c,f,y1,n+1)


        
        do i=1,n+1
        if(y1(i).lt.0.d0) then
        if((y1(i) + epsilon(y1(i))).gt.0.d0) then
        y1(i)=0.d0
        else
        write(*,*)'y(i)=',y1(i),' lt negative epsilon=',epsilon(y1(i))
        pause
        stop
        endif
        endif
        enddo
           
        enddo

      do i=1,n+1
        y(i+1)=y1(i)
      end do
        
        do i=1,i0
         if(vij(i,j).lt.xend) then
          call lock(x,n+2,vij(i,j),klo,khi,ierr)
          if(ierr.eq.1) then
           write(*,*)'lock error #4 in finction fokkerplanck'
           write(*,*)'j=',j,' vij=',vij(i,j)
           write(*,*)'x(1)=',x(1),' x(n+2)=',x(n+2)
           pause
           stop
          end if
          call linf(x,y,vij(i,j),fij(i,j,k),klo,khi)
         else
          fij(i,j,k)=zero
         end if
         
        enddo
        
        
!        if((tim1 > 0.03109d0).and.(pachka.eq.10)) then
!c        call integral(1,i0,vij(:,j),fij(:,j,k),fout1)
!        if((k.eq.1).and.(j.eq.nr)) then
!            open(8,file='RESULT2019/xyhdt1.dat',position="append")
!            open(11,file='RESULT2019/xxd1.dat',position="append")            
!c            write(8,10) tim1, fout1
!            do i=1,n+2
!            write(8,10) x(i), y(i), h, dt
!            enddo
!            
!            do i=1,n+1
!            write(11,13) xx(i), D1(i), alfa2
!            enddo
!            
!            close(8)
!            close(11)
!        elseif((k.eq.2).and.(j.eq.nr)) then
!            open(9,file='RESULT2019/xyhdt2.dat',position="append")
!            open(12,file='RESULT2019/xxd2.dat',position="append")
!c            write(9,10) tim1, fout1
!
!            do i=1,n+2
!            write(9,10) x(i), y(i), h, dt
!            enddo
!            
!            do i=1,n+1
!            write(12,13) xx(i), D1(i), alfa2
!            enddo
!            close(9)
!            close(12)
!        endif
!        
!        pause
!        endif
10      format(4F20.16)
13      format(3F20.16)



        deallocate(y,x,xx,a,b,c,f,D1,D2,y1)
        allocate(fj(i0),dfj(i0))

        do i=1,i0
         fj(i)=fij(i,j,k)
         dfj(i)=zero
        end do

        do i=1,i0
         if(i.eq.1) then
          dfj(i)=zero
         else if(i.eq.i0) then
          dfj(i)=(fj(i)-fj(i-1))/vij(2,j)
         else
          dfj(i)=0.5d0*(fj(i+1)-fj(i-1))/vij(2,j)
         end if
        end do
        ii=0
        ibeg=0

c	do i=i0-1,1,-1
c         if(dfj(i).gt.zero) then
cc          write(*,*) '#2 positive derivs'
cc          write(*,*) '#2 df>0: i,j,k=',i,j,k
cc          write(*,*) '#2 dfj(i),i,j,k=',dfj(i),i,j,k
cc          write(*,*) '#2 fij=',fij(i,j,k)
cc          write(*,*)
c          fij(i,j,k)=fij(i+1,j,k)
c          dfij(i,j,k)=dfij(i+1,j,k)
c          ii=i
c         end if
c         if(fij(i,j,k).lt.fij(i+1,j,k)) then
c          fij(i,j,k)=fij(i+1,j,k)
c          dfij(i,j,k)=dfij(i+1,j,k)
c          ii=i
c         end if
c        end do
        ibeg=ii

c        if(ibeg.gt.0) then
c         call integral(ibeg,i0,vj,fj,fout1)
c         do i=1,i0
c          fj(i)=fij(i,j,k)
c          dfj(i)=dfij(i,j,k)
c         end do
c         call integral(ibeg,i0,vj,fj,fout2)
c         do i=ibeg,i0
c          fij(i,j,k)=fj(i)*fout1/fout2
c          dfij(i,j,k)=dfj(i)*fout1/fout2
c         end do
c!!         write(*,*)'#2 j,k,ibeg=',j,k,ibeg
c!!         write(*,*)'#2 v(ibeg)=',vj(ibeg),' f1/f2=',fout1/fout2
c        end if


        deallocate(vj,fj,dfj)

        enddo
        
       end do



      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
      subroutine ABCcoef(A,B,C,f,Y,k,n,ybeg,yend,x,xx,h,df)
      implicit none
      integer i,n,jindex,kindex
      real*8 f(n),Y(n),k,df(n),d0
      real*8 ybeg,yend,x(*),h,xx(*)
      real*8 C1,B1,z,w,r,dlt,alfa2
      real*8 tmp1,tmp2,tmp3,A(n),B(n),C(n)
      common/dddql/ d0,jindex,kindex
      common/ef/ alfa2 
      external C1,w,B1,dlt


        r=k/h
       do i=1,n-1
      
      tmp1=dlt(xx(i),h,df(i))*B1(xx(i))
      A(i)=-r*(C1(xx(i),df(i))/h-tmp1)

      tmp2=C1(xx(i+1),df(i+1))/h-dlt(xx(i+1),h,df(i+1))*B1(xx(i+1))
      tmp3=(1.d0-dlt(xx(i),h,df(i)))*B1(xx(i))
      B(i)=r*(tmp2+tmp3+C1(xx(i),df(i))/h)+1.d0
	
      tmp1=(1.d0-dlt(xx(i+1),h,df(i+1)))*B1(xx(i+1))
      C(i)=-r*(tmp1+C1(xx(i+1),df(i+1))/h)

            f(i)=Y(i)
      enddo


 
        f(1)=f(1)-A(1)*ybeg
        f(n)=Y(n)
        A(n)=-r*(C1(xx(n),df(n))/h-dlt(xx(n),h,df(n))*B1(xx(n)))
        B(n)=r*((1.d0-dlt(xx(n),h,df(n)))*B1(xx(n))+C1(xx(n),df(n))/h)+1.d0
        C(n)=0.d0


 
        !f(1)=f(1)-A(1)*ybeg
        !f(n)=f(n)-C(n)*yend !yend in either way=0 all the time

c        C(1)=0.d0
c        B(1)=1.d0
c        A(N)=0.d0
c        B(N)=1.d0

        
c        if (kindex.eq.1) then
c        tmp1=(1.d0-dlt(xx(n-1),h,df(n-1)))*B1(xx(n-1))
c        B(N)=r*(tmp1+C1(xx(n-1),df(n-1))/h)+1.d0
        
c        tmp1=dlt(xx(n-1),h,df(n-1))*B1(xx(n-1))
c        A(n)=-r*(C1(xx(n-1),df(n-1))/h-tmp1)

c        tmp1=(1.d0-dlt(xx(1),h,df(1)))*B1(xx(1))
c        C(1)=-r*(tmp1+C1(xx(1),df(1))/h)

c        tmp1=dlt(xx(1),h,df(1))*B1(xx(1))
c        tmp2=(1-dlt(xx(1),h,df(1)))*B1(xx(1))+C2(xx(1),df(1))/h
c        B(1)=r*(C1(xx(1),df(1))/h-tmp1)+1.d0
cc        elseif (kindex.eq.2) then
c        tmp1=(1.d0-dlt(xx(n-1),h,df(n-1)))*B1(xx(n-1))
cc        tmp2=C1(xx(n-1)+h,0.d0)/h-dlt(xx(n-1)+h,h,0.d0)*B1(xx(n-1)+h)
c        B(N)=r*(tmp1+C1(xx(n-1),df(n-1))/h)+1.d0
c        C(1)=0.d0
c        B(1)=1.d0
c        f(1)=ybeg
c        else
c        write(*,*) 'error boundary cond. drvcrrnt 830'
c        endif
c        tmp1=dlt(x(1)+h/2.d0,h,df(1))*B1(x(1)+h/2.d0)
c        B(1)=r*(C1(x(1)+h/2.d0,df(1))/h-tmp1)+1.d0
c        tmp1=(1.d0-dlt(x(n)-h/2.d0,h,df(n-1)))*B1(x(n)-h/2.d0)
c        B(N)=r*(tmp1+C1(x(n)-h/2.d0,df(n-1))/h)+1.d0
        end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        real*8 function B1(xx)
        implicit none
        real*8 xx,alfa2,beta
        common/ef/ alfa2 
            B1=-alfa2+1.d0/(xx*xx)
        end function

        
        real*8 function C1(xx,dif)
        implicit none
        real*8 xx,dif
            C1=dif+1.d0/(xx*xx*xx)
        end function
        
        
        real*8 function w(xx,h,dif)
        implicit none
        real*8 xx,h,B1,C1,dif
           w=h*B1(xx)/C1(xx,dif)
        end function
        
        
        real*8 function dlt(xx,h,dif)
        implicit none
        real*8 xx,h,B1,C1,w,dif
c        if(w(xx,h,dif).gt.700.d0) then
c        dlt=1.d0/w(xx,h,dif)
c        elseif(w(xx,h,dif).lt.-730.d0) then
c        dlt=1.d0/w(xx,h,dif)+1.d0
c        else
        dlt=1.d0/w(xx,h,dif)-1.d0/(dexp(w(xx,h,dif))-1.d0)
c        endif
        end function


c        real*8 function z(xx,h,dif)
c        implicit none
c        real*8 xx,h,w,dif
c        if(w(xx,h,dif).lt.0d0) then
c            z=w(xx,h,dif)/(dexp(w(xx,h,dif))-1d0)
c        elseif(w(xx,h,dif).ge.0d0) then
c            z=dexp(-w(xx,h,dif))*w(xx,h,dif)/(1d0-dexp(-w(xx,h,dif)))
c        endif
c        end function

!
!      subroutine abccoef(a,b,c,f,y,dt,n,ybeg,yend,x,h
!     *,dif,ddif)
!      implicit none
!      integer i,n,jindex,kindex 	
!      real*8 a(n),b(n),c(n),f(n),y(n),dt,dif(n),ddif(n)
!      real*8 ybeg,yend,x(*),h,temp
!      real*8 a1,b1,c1,tmp1,tmp2,dDdX,tcur,d0
!      external a1,b1,c1,dDdX
!      common/dddql/ d0,jindex,kindex
!      common/b1dif/ temp

!      do i=1,n
!    
!      a(i)=-(2d0*dt*a1(x(i),dif(i))-dt*h*b1(x(i),dif(i),ddif(i)))
!    
!      b(i)=4d0*h**2+4d0*dt*a1(x(i),dif(i))-2d0*h**2*dt*c1(x(i))
!   
!      c(i)=-(2d0*dt*a1(x(i),dif(i))+dt*h*b1(x(i),dif(i),ddif(i)))

 
!          if(i.eq.1) then
!      f(1)=ybeg
!             else if(i.eq.n) then
!      f(n)=yend
!             else
!      tmp1=(4d0*h**2-4d0*dt*a1(x(i),dif(i))
!     * +2d0*h**2*dt*c1(x(i)))*y(i)
!      tmp2=(2d0*dt*a1(x(i),dif(i))
!     *-dt*h*b1(x(i),dif(i),ddif(i)))*y(i-1)
!      f(i)=(2d0*dt*a1(x(i),dif(i))
!     *+dt*h*b1(x(i),dif(i),ddif(i)))*y(i+1)+tmp1+tmp2
!             end if
!      enddo
!	c(1)=0d0
!	b(1)=1
 !       a(n)=0d0
!	b(n)=1
!	end
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!FUNCTION a1

!      real*8 function a1(x,dif)
!      implicit none
!      real*8 x,t,tlock,d,dif
!      !common /time/ tlock
!
!	a1=dif+1d0/x**3
 !     end

!FUNCTION b1

!	real*8 function b1(x,dif,ddif)
 !     implicit none
 !     real*8 x,t,tlock,dDdX,d,h,alfa2
!      real*8 dif,ddif,temp
!      !common /time/ tlock
!      common/ef/ alfa2 
!      common/b1dif/ temp

!	b1=(1d0/x**2-3d0/x**4)-alfa2+ddif

 !     end

!FUNCTION c1

!	real*8 function c1(x)
 !     implicit none
 !     real*8 x,t,d
 !     !common /time/ tlock
!	c1=-2d0*1d0/x**3
 !     end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine tridag(a,b,c,r,u,n)
      implicit none
      integer n,nmax
	real*8 r(n),u(n)
      real*8 a(n),b(n),c(n),eps
      parameter (nmax=1000000,eps=1.d-200)
      integer j
      double precision bet,gam(nmax)


      if(b(1).eq.0.d0)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.d0)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue

	return

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine polint(xa,ya,n,x,y,dy)
      implicit none
      integer n,nmax
      double precision dy,x,y,xa(n),ya(n)
      parameter (nmax=10)
      integer i,m,ns
      double precision den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
          dift=abs(x-xa(i))
          if (dift.lt.dif) then
              ns=i
              dif=dift
          endif
          c(i)=ya(i)
          d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
          do 12 i=1,n-m
              ho=xa(i)-x
              hp=xa(i+m)-x
              w=c(i+1)-d(i)
              den=ho-hp
              if(den.eq.0.d0)pause 'failure in polint'
              den=w/den
              d(i)=hp*den
              c(i)=ho*den
12        continue
          if (2*ns.lt.n-m)then
              dy=c(ns+1)
          else
              dy=d(ns)
              ns=ns-1
          endif
          y=y+dy
13    continue
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine trapzd(func,a,b,s,n)
      implicit none
      integer n
      double precision a,b,s,func
      external func
      integer it,j
      double precision del,sum,tnm,x
      if (n.eq.1) then
          s=0.5d0*(b-a)*(func(a)+func(b))
      else
          it=2**(n-2)
          tnm=it
          del=(b-a)/tnm
          x=a+0.5d0*del
          sum=0.d0
          do 11 j=1,it
             sum=sum+func(x)
              x=x+del
11        continue
          s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      double precision function aiint(arr,yr)
c aiint:	integral {0,r} of current density
c 	aiint=integral {0,r} (arr/ipol**2)dv*ipol/(gp2*ro)
c 	only radial dependent array may be a parameter of the function
c 	examples:
c    out_iint(cu)	!radial profile of toroidal current
c    out_iint(cd,ro)    !toroidal driven current inside {0,ro}
c    out_iint(cub)      !total toroidal current =iint(cu,roc); (=ipl)
c			(pereverzev 23-oct-99)
      implicit none
      integer j,jk
      double precision arr(*),yr,dr,ya
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      aiint=0.d0
      if(yr.le.0.d0) return
      jk=yr/hro+1.d0-1.d-4
      if(jk.gt.na) jk=na
      ya=0.d0
      do j=1,jk
       aiint=aiint+ya
       ya=arr(j)*rho(j)/(g33(j)*ipol(j)**3)
       dr=yr-jk*hro+hro
      end do
      if(jk.ge.na) then
       dr=min(dr,0.5d0*(hro+hroa))
      endif
      aiint=gp2*ipol(jk)*(hro*aiint+dr*ya)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine integral11(f,n,h1,fout,fbeg,fend)
      integer i,n
      Real*8 h1, summ, fout, f(n),fbeg,fend

      fout=0.d0
      summ=0.d0
      do i=1,n+2
          if(i.eq.1) then
                summ=summ+5.d0*fbeg/12.d0
          elseif(i.eq.2) then
                summ=summ+13.d0*f(1)/12.d0
          elseif(i.eq.(n+1)) then
                summ=summ+13.d0*f(n)/12.d0
          elseif(i.eq.(n+2)) then
                summ=summ+5.d0*fend/12.d0
          else
                summ=summ+f(i)
          endif
      end do
      fout=summ*h1
      end

      subroutine part_output(x,f,t,n,k,fout)
      integer n,k
      real*8 x,f,t,fout
      
      call integral(1,n,x,f,fout1)

      if(k.eq.1) then
            open(8,file='RESULT2019/part1.dat',position="append")
            write(8,10) t, fout1
            close(8)
      elseif (k.eq.2) then
            open(9,file='RESULT2019/part2.dat',position="append")
            write(9,10) t, fout1
            close(9)
      endif
        
10      format(2F20.16)
      end