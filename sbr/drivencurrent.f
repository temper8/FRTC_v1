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
       call fokkerplanck(dt,time,i)
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
      subroutine fokkerplanck(dtstep,time,nomer)
      implicit none
      real*8 t,dtstep,dtau,d5,check1,check2,check3
      integer nr,ispl,ip,im,ii,ibeg,nomer,vivod,iunit
      common /a0ab/ nr
      real*8 ynzm0,pm0,plaunp,plaunm,fmaxw,time,pachka
      common/grillspektr/ ynzm0(1001),pm0(1001),ispl
     &,plaunp,plaunm,ip,im
      integer i0
      parameter(i0=1002)
      real*8 vij,fij0,fij,dfij,dij,enorm,fst
      common/lh/vij(i0,100),fij0(i0,100,2),fij(i0,100,2),dfij(i0,100,2)
     &,dij(i0,100,2),enorm(100),fst(100)
      integer n,i,j,it,nt,k
      real*8 xend,h,shift,ybeg,yend,tend,dt,dff
      real*8,dimension(:),allocatable:: y,x,xx,xxm,xxp,a,b,c,f
      real*8,dimension(:),allocatable:: vj,fj,dfj,d1,d2,d3,givi
      real*8 znak,alfa2,zero,dt0,h0,eps,r,fvt,fout1,fout2
      common/ef/ alfa2
      real*8 calls
      common/firstcall/calls
      real*8 d0
      integer jindex,kindex,klo,khi,ierr,klo1,khi1
      integer klo2,klo3,khi2,khi3,ierr1,ierr2,ierr3
      common/dddql/ d0,jindex,kindex
      parameter(zero=0.d0,dt0=0.1d0,h0=0.1d0,eps=1.d-7)
!
      do k=1,2
       kindex=k
       znak=2.d0*dble(k)-3.d0
!
!!       do j=1,nr
!!        fij0(i,j,k)=fmaxw(vij(i,j),znak*enorm(j),dff)
!!       end do
!!       if(1.gt.0) go to 1
!
       d0=zero
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
        allocate(y(n),x(n+2),xx(n+1),a(n),b(n),c(n),f(n))
!!!!!! grid !!!!!!!!!
!!       shift=h*0.1d0 !0.01d0
        do i=1,n+2
         x(i)=h*dble(i-1) !+shift
        end do
        do i=1,n+1
         xx(i)=h/2.d0+h*dble(i-1) !+shift
        end do
        allocate(vj(i0),fj(i0))
        do i=1,i0
         vj(i)=vij(i,j)
         fj(i)=fij0(i,j,k)
        end do
        do i=1,n
         call lock(vj,i0,x(i+1),klo,khi,ierr)
         if(ierr.eq.1) then
          write(*,*)'lock error #1 in finction fokkerplanck'
          write(*,*)'j=',j,' v=',x(i+1)
          write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
          pause
          stop
         end if
         call linf(vj,fj,x(i+1),y(i),klo,khi)
        end do
        deallocate(fj)
        ybeg=fij0(1,j,k)  !boundary conditions
        yend=zero
        alfa2=znak*enorm(j)
!!!!!!!!!!!!!EVALUATING DIFFUSION!!!!!!!!!!!!!!!!!!

	allocate(d1(n+1),d2(n+1),d3(n+1))
	do i=1,n+1
	d1(i)=0d0
	d2(i)=0d0
	d3(i)=0d0
	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do it=1,nt
         call abccoef(a,b,c,f,y,dt,n,ybeg,yend,xx,h,d1,d2,d3)
         call tridag(a,b,c,f,y,n)
!!         t=dt*dble(it)
        end do
        deallocate(d1,d2,d3)
	pachka=0.35
	if(time*nomer.gt.pachka)then
	open(iunit,file='lhcd/pressF.dat',position="append")
	do vivod=1,n
!!	write(iunit,*)time,vivod,f(vivod)
	end do
	end if
        allocate(fj(n+2))
        fj(1)=ybeg
        fj(n+2)=yend
        do i=1,n
         fj(i+1)=y(i)
        end do
        do i=2,i0-1
         if(vij(i,j).lt.xend) then
          call lock(x,n+2,vij(i,j),klo,khi,ierr)
          if(ierr.eq.1) then
           write(*,*)'lock error #2 in finction fokkerplanck'
           write(*,*)'j=',j,' vij=',vij(i,j)
           write(*,*)'x(1)=',x(1),' x(n+2)=',x(n+2)
           pause
           stop
           end if
          call linf(x,fj,vij(i,j),fij0(i,j,k),klo,khi)
         else
          fij0(i,j,k)=zero
         end if
        end do
        deallocate(fj)
        deallocate(y,x,xx,a,b,c,f)
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
        do i=i0-1,1,-1
         if(dfj(i).gt.zero) then
c          write(*,*) '#1 positive derivs'
c          write(*,*) '#1 df>0: i,j,k=',i,j,k
c          write(*,*) '#1 dfj(i),i,j,k=',dfj(i),i,j,k
c          write(*,*)
          fij0(i,j,k)=fij0(i+1,j,k)
          ii=i
         end if
         if(fij0(i,j,k).lt.fij0(i+1,j,k)) then 
          fij0(i,j,k)=fij0(i+1,j,k)
          ii=i
         end if
        end do
        ibeg=ii
!
        if(ibeg.gt.0) then
         call integral(ibeg,i0,vj,fj,fout1)
         do i=1,i0
          fj(i)=fij0(i,j,k)
         end do
         call integral(ibeg,i0,vj,fj,fout2)
         do i=ibeg,i0
          fij0(i,j,k)=fj(i)*fout1/fout2
         end do
!!         write(*,*)'#1 j,k,ibeg=',j,k,ibeg
!!         write(*,*)'#1 v(ibeg)=',vj(ibeg),' f1/f2=',fout1/fout2
        end if
!
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
        allocate(y(n),x(n+2),xx(n+1),a(n),b(n),c(n),f(n))
	allocate(xxm(n+1),xxp(n+1))
!!	write(*,*)'OSHIBKA2, j=',j
!!!!!! grid !!!!!!!!!
!!       shift=h*0.1d0 !0.01d0
        do i=1,n+2
         x(i)=h*dble(i-1) !+shift
        end do
        do i=1,n+1
         xx(i)=h/2.d0+h*dble(i-1) !+shift
!	 xxm(i)=xx(i)-h/2d0
!	 xxp(i)=xx(i)+h/2d0
        end do
	if (nomer.gt.9) then 
	open(iunit,file='RESULT2019/xxx.dat',position="append")
!	check1=xx(4)-h/2d0-xx(4)+h/2d0
!	check2=xx(5)+h/2d0-xxp(5)
!	check3=xx(5)-h/2d0-xxm(5)
!	do i=1,n
!	write(iunit,*)xx(1)-h/2d0-xxx(1),xx(1)+h/2d0-xxx(2)
!	write(iunit,*)xx(2)-h/2d0-xxx(2),xx(2)+h/2d0-xxx(4)
!	write(iunit,*)xx(2)-xxx(3),xx(1)-xxx(1)
!	write(iunit,*)check1,check2,check3
!	end do
	end if
	close(iunit)
        allocate(vj(i0),fj(i0))
        do i=1,i0
         vj(i)=vij(i,j)
         fj(i)=fij(i,j,k)
        end do
        do i=1,n
         call lock(vj,i0,x(i+1),klo,khi,ierr)
         if(ierr.eq.1) then
          write(*,*)'lock error #3 in finction fokkerplanck'
          write(*,*)'j=',j,' v=',x(i+1)
          write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
          pause
          stop
         end if
         call linf(vj,fj,x(i+1),y(i),klo,khi)
        end do
        deallocate(fj)
        ybeg=fij(1,j,k)  !boundary conditions
        yend=zero
        alfa2=znak*enorm(j)
!!!!!!!!!!!!!EVALUATING DIFFUSION!!!!!!!!!!!!!!!!!!
!!	write(*,*)'OSHIBKA, j=',j
	if (nomer.gt.9) then 
	open(iunit,file='RESULT2019/dddd.dat',position="append")
	end if
        allocate(d1(n+1),d2(n+1),d3(n+1))	
	do i=1,n+1
	call lock(vj,i0,xx(i),klo1,khi1,ierr1)
	call lock(vj,i0,xx(i)-h/2d0,klo2,khi2,ierr2)
	call lock(vj,i0,xx(i)+h/2d0,klo3,khi3,ierr3)
      	if(ierr1.eq.1) then
	write(*,*)'lock error in finction d2(x)'
	write(*,*)'j=',j,' v=',xx(i)
	write(*,*)'klo1=',klo1,'khi1=',khi1,'i=',i
	write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
	pause
	stop
      	end if
      	if(ierr2.eq.1) then
	write(*,*)'lock error in finction d2(x)'
	write(*,*)'j=',j,' v=',xxm(i)
	write(*,*)'klo2=',klo2,'khi2=',khi2,'i=',i
	write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
	pause
	stop
      	end if
      	if(ierr3.eq.1) then
	write(*,*)'lock error in finction d2(x)'
	write(*,*)'j=',j,' v=',xxp(i)
	write(*,*)'klo3=',klo3,'khi3=',khi3,'i=',i
	write(*,*)'vj(1)=',vj(1),' vj(i0)=',vj(i0)
	pause
	stop
      	end if
	d1(i)=dij(klo1,j,k)
	d2(i)=dij(klo2,j,k)
	d3(i)=dij(klo3,j,k)	
	if (nomer.gt.9) then 
!	write(iunit,*)d1(i),d2(i),d3(i)
	end if
	end do
	if (nomer.gt.9) then 
	close(iunit)
	end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!solve problem
        do it=1,nt
         call abccoef(a,b,c,f,y,dt,n,ybeg,yend,xx,h,d1,d2,d3)
         call tridag(a,b,c,f,y,n)
!!         t=dt*dble(it)
        end do
        deallocate(d1,d2,d3)
        allocate(fj(n+2))
        fj(1)=ybeg
        fj(n+2)=yend
        do i=1,n
         fj(i+1)=y(i)
        end do
	if (nomer.gt.9) then 
	open(iunit,file='RESULT2019/distribution.dat',position="append")
	do i=1,n
!	write(iunit,*)i,fj(i+1)
	end do
	end if
	close(iunit)
        do i=2,i0-1
         if(vij(i,j).lt.xend) then
          call lock(x,n+2,vij(i,j),klo,khi,ierr)
          if(ierr.eq.1) then
           write(*,*)'lock error #4 in finction fokkerplanck'
           write(*,*)'j=',j,' vij=',vij(i,j)
           write(*,*)'x(1)=',x(1),' x(n+2)=',x(n+2)
           pause
           stop
          end if
          call linf(x,fj,vij(i,j),fij(i,j,k),klo,khi)
         else
          fij(i,j,k)=zero
         end if
        end do
        deallocate(fj)

        deallocate(y,x,xx,xxm,xxp,a,b,c,f)
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
        do i=i0-1,1,-1
         if(dfj(i).gt.zero) then
c          write(*,*) '#2 positive derivs'
c          write(*,*) '#2 df>0: i,j,k=',i,j,k
c          write(*,*) '#2 dfj(i),i,j,k=',dfj(i),i,j,k
c          write(*,*) '#2 fij=',fij(i,j,k)
c          write(*,*)
          fij(i,j,k)=fij(i+1,j,k)
          dfij(i,j,k)=dfij(i+1,j,k)
          ii=i
         end if
         if(fij(i,j,k).lt.fij(i+1,j,k)) then
          fij(i,j,k)=fij(i+1,j,k)
          dfij(i,j,k)=dfij(i+1,j,k)
          ii=i
         end if
        end do
        ibeg=ii
!
        if(ibeg.gt.0) then
         call integral(ibeg,i0,vj,fj,fout1)
         do i=1,i0
          fj(i)=fij(i,j,k)
          dfj(i)=dfij(i,j,k)
         end do
         call integral(ibeg,i0,vj,fj,fout2)
         do i=ibeg,i0
          fij(i,j,k)=fj(i)*fout1/fout2
          dfij(i,j,k)=dfj(i)*fout1/fout2
         end do
!!         write(*,*)'#2 j,k,ibeg=',j,k,ibeg
!!         write(*,*)'#2 v(ibeg)=',vj(ibeg),' f1/f2=',fout1/fout2
        end if
!
        deallocate(vj,fj,dfj)

       end do
      end do
!
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine abccoef(a,b,c,f,y,dt,n,ybeg,yend,xx,h,d1,d2,d3)
      implicit none
      integer i,n,iunit,iunit2
      real*8 a(n),b(n),c(n),f(n),y(n),d1(n+1),d2(n+1),d3(n+1)
      real*8 a1(n),b1(n),c1(n),f1(n),a2(n),b2(n),c2(n),f2(n)
      real*8 dt,kinv,rs,rmink,rplusk,q,qf,r1,rmink2,rplusk2,kinv2
      real*8 ybeg,yend,xx(*),h,r,kappa,sum,bmin,bplus,sum2,sum3,sum4
      real*8 dc,as(n+1),k,k2,d
      external kinv,rs,rmink,rplusk,q,kinv2,rmink2,rplusk2,d

      sum=(kinv(xx(1)-h/2d0,d2(1))+kinv(xx(1)+h/2d0,d3(1)))*h/2d0
      as(1)=h/sum

      sum=(kinv(xx(2)-h/2d0,d2(2))+kinv(xx(2)+h/2d0,d3(2)))*h/2d0
      as(2)=h/sum

      r=h/2d0*dabs(rs(xx(1)+h/2d0,d3(1)))/k(xx(1)+h/2d0,d3(1))
      kappa=1d0/(1d0+r)
      sum=(rmink(xx(1),d1(1))+rmink(xx(2),d1(2)))*h/2d0
      bmin=sum/h

      sum=(rplusk(xx(1),d1(1))+rplusk(xx(2),d1(2)))*h/2d0
      bplus=sum/h

      sum=qf(xx(2))-qf(xx(1))
      dc=sum/h

      a(1)=as(1)*(kappa/h**2-bmin/h)
      c(1)=as(2)*(kappa/h**2+bplus/h)
      b(1)=-(1d0/dt+a(1)+c(1)+dc)
      f(1)=-y(1)/dt-a(1)*ybeg
	open(iunit,file='RESULT2019/RAZNICA.dat',position="append")
      do i=2,n
!       sum=1d0      !(kinv(xx(i+1)-h/2d0,d2(i+1)))+kinv(xx(i+1)+h/2d0,d3(i+1)))*h/2d0
!      sum2=(rplusk(xx(i),d1(i))+rplusk(xx(i+1),d1(i+1)))*h/2d0
      !sum2=sum2*h/2d0
!      sum3=(rplusk2(xx(i),d1(i))+rplusk2(xx(i+1),d1(i+1)))*h/2d0
      !sum3=sum3*h/2d0
!	sum4=sum3-sum2
!        if(dabs(sum4).gt.0d0) then
!	write(iunit,*)sum2,sum3,sum4
!	end if
       sum=(kinv(xx(i+1)-h/2d0,d2(i+1))+kinv(xx(i+1)+h/2d0,d3(i+1)))
       sum=sum*h/2d0
       as(i+1)=h/sum
       r=h/2d0*dabs(rs(xx(i)+h/2d0))/k(xx(i)+h/2d0,d3(i))
       kappa=1d0/(1d0+r)
       sum=(rmink(xx(i),d1(i))+rmink(xx(i+1),d1(i+1)))*h/2d0
       bmin=sum/h
       sum=(rplusk(xx(i),d1(i))+rplusk(xx(i+1),d1(i+1)))*h/2d0
       bplus=sum/h
       sum=qf(xx(i+1))-qf(xx(i))
       dc=sum/h
       a(i)=as(i)*(kappa/h**2-bmin/h)
       c(i)=as(i+1)*(kappa/h**2+bplus/h) 
       b(i)=-(1d0/dt+a(i)+c(i)+dc) 
       f(i)=-y(i)/dt

!       sum=(kinv2(xx(i+1)-h/2d0,d2(i+1))+kinv2(xx(i+1)+h/2d0,d3(i+1)))
!       sum=sum*h/2d0
!       as(i+1)=h/sum
!       r=h/2d0*dabs(rs(xx(i)+h/2d0))/k2(xx(i)+h/2d0,d3(i))
!       kappa=1d0/(1d0+r)
!       sum=(rmink2(xx(i),d1(i))+rmink2(xx(i+1),d1(i+1)))*h/2d0
!       bmin=sum/h
!       sum=(rplusk2(xx(i),d1(i))+rplusk2(xx(i+1),d1(i+1)))*h/2d0
!       bplus=sum/h
!       sum=qf(xx(i+1))-qf(xx(i))
!       dc=sum/h
!       a1(i)=as(i)*(kappa/h**2-bmin/h)
!       c1(i)=as(i+1)*(kappa/h**2+bplus/h) 
!       b1(i)=-(1d0/dt+a(i)+c(i)+dc) 
!       f1(i)=-y(i)/dt

!       a2(i)=a1(i)-a(i)
!       c2(i)=c1(i)-c(i)
!       b2(i)=b1(i)-b(i)
!       f2(i)=f1(i)-f(i)
!	open(iunit2,file='RESULT2019/abcr.dat',position="append")
!        if(dabs(c2(i)).gt.0d0) then
!	write(iunit,*)f2(i),f1(i),f(i)
!	end if
!	close(iunit2)
      end do
	close(iunit)
      f(n)=f(n)-c(n)*yend
      a(1)=0d0
      c(n)=0d0
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function rplusk(x,dif)
      implicit none
      integer iunit
      real*8 x,k,rs,dif,d,razn
      rplusk=0.5d0*(rs(x)+dabs(rs(x)))/k(x,dif)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function rplusk2(x,dif)
      implicit none
      integer iunit
      real*8 x,k2,rs,dif,d,razn
      rplusk2=0.5d0*(rs(x)+dabs(rs(x)))/k2(x,dif)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function rmink(x,dif)
      implicit none
      integer iunit
      real*8 x,k,rs,dif,d,razn
      rmink=0.5d0*(rs(x)-dabs(rs(x)))/k(x,dif)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function rmink2(x,dif)
      implicit none
      integer iunit
      real*8 x,k2,rs,dif,d,razn
      rmink2=0.5d0*(rs(x)-dabs(rs(x)))/k2(x,dif)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function rs(x)
      implicit none
      real*8 x
      real*8 alfa2
      common/ef/ alfa2
       rs=1d0/x**2-alfa2
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function q(x)
      implicit none
      real*8 x
      q=2d0/x**3
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function qf(x)
      implicit none
      real*8 x
      qf=-1d0/x**2
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function k(x,dif)
      implicit none
      integer iunit
      real*8 x,dif,d,razn
!	razn=dif-d(x)
!	open(iunit,file='RESULT2019/kkk.dat',position="append")
!	if(dabs(razn).gt.0d0) then
!	write(iunit,*)razn
!	end if
!	close(iunit)

      k=dif+1d0/x**3
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function k2(x,dif)
      implicit none
      integer iunit
      real*8 x,dif,d,razn
!	razn=dif-d(x)
!	open(iunit,file='RESULT2019/kkk.dat',position="append")
!	if(dabs(razn).gt.0d0) then
!	write(iunit,*)razn
!	end if
!	close(iunit)

      k2=d(x)+1d0/x**3
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function kinv(x,dif)
      implicit none
      integer iunit
      real*8 x,dif,razn,d
      kinv=x**3/(dif*x**3+1d0)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function kinv2(x,dif)
      implicit none
      integer iunit
      real*8 x,dif,razn,d,kino
      kinv2=x**3/(d(x)*x**3+1d0)
!	kino=x**3/(dif*x**3+1d0)
!	razn=kino-kinv2
!	open(iunit,file='RESULT2019/kino.dat',position="append")!
!	if(dabs(razn).gt.0d0) then
!	write(iunit,*)razn
!	end if
!	close(iunit)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ddc(diffusion)
      implicit none
      common/testf/ tcur

      integer ntau,tc,koltoch,i,j,k,klo
      real*8 curtime,tcur,zero,dij
      real*16 tau0,spacing,curtime0
      parameter(tau0=3.000990745207882E-002, zero=0.d0)
      common/lh/dij(1002,100,2)
      real*8 b,b1,b2,d,diffusion
!      real*8,dimension(:),allocatable:: diffusion
      integer i1,iunit6
      b1=0
      b2=60
      j=10
      k=1
!      allocate (diffusion(500))
! write(*,*)'time=',tcur
! do tc=1,koltoch
!       spacing=0.008/koltoch
!       curtime=tau0+spacing*tc!-0.0002
!      curtime0=curtime+0.000000001
!      if((tcur-curtime)*(tcur-curtime0).lt.zero) then
!       if((tcur-0.0301)*(tcur-0.0302).lt.zero) then
!       open(iunit6,file='RESULT2019/ddc.dat',position="append")
!      do i1=1,500
!      b=b1+(b2/500)*(i1-1)
!      diffusion(i1)=d(b)
!       do i=1,1001
!        write(iunit6,*) i, diffusion
!       end do
!      write(iunit6,*)
!       close(iunit6)
!       end if
! end do
! !     deallocate(diffusion)
      end














      real*8 function d(x)
      implicit none
      integer i0
      parameter(i0=1002)
      real*8 vij,fij0,fij,dfij,dij,enorm,fst
      common/lh/vij(i0,100),fij0(i0,100,2),fij(i0,100,2),dfij(i0,100,2)
     &,dij(i0,100,2),enorm(100),fst(100)
      real*8,dimension(:),allocatable:: vvj,ddj
      integer klo,khi,ierr
      real*8 d0,zero,x
      integer jindex,kindex,k,j,i
      common/dddql/ d0,jindex,kindex
      parameter(zero=0.d0)
      d=zero
      if(d0.eq.zero) return
      j=jindex
      if(x.ge.vij(i0,j)) return
      k=kindex
      allocate(vvj(i0),ddj(i0))
!      write(*,*)'function d(x): k=',k,' j=',j
      do i=1,i0
       vvj(i)=vij(i,j)
       ddj(i)=dij(i,j,k)
      end do
      call lock(vvj,i0,x,klo,khi,ierr)
      if(ierr.eq.1) then
       write(*,*)'lock error in finction d2(x)'
       write(*,*)'j=',j,' v=',x
       write(*,*)'vj(1)=',vvj(1),' vj(i0)=',vvj(i0)
       pause
       stop
      end if
      d=ddj(klo)
!!	call ddc(d)
!!!      call linf(vvj,ddj,x,d,klo,khi)
!
!      write(*,*)'klo=',klo,' khi=',khi
!      write(*,*)'vj(klo)=',vvj(klo),' vj(khi)=',vvj(khi)
!      write(*,*)'x=',x,' d=',d
!
      deallocate(vvj,ddj)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine tridag(a,b,c,r,u,n)
      implicit none
      integer n,nmax
      double precision a(n),b(n),c(n),r(n),u(n)
      parameter (nmax=1000000)
      integer j
      double precision bet,gam(nmax)
      if(b(1).eq.0.d0)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.d0) then
	write(*,*)'b(j)=',b(j),'a(j)=',a(j),'gam(j)=',gam(j)
	pause 'tridag failed'
	end if
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine qromb(func,a,b,ss)
      implicit none
      integer jmax,jmaxp,k,km
      double precision a,b,func,ss,eps
      external func
      parameter (eps=1.d-6, jmax=200, jmaxp=jmax+1, k=5, km=k-1)
cu    uses polint,trapzd
      integer j
      double precision dss,h(jmaxp),s(jmaxp)
      h(1)=1.d0
      do 11 j=1,jmax
          call trapzd(func,a,b,s(j),j)
          if (j.ge.k) then
              call polint(h(j-km),s(j-km),k,0.d0,ss,dss)
              if (abs(dss).le.eps*abs(ss)) return
          endif
          s(j+1)=s(j)
          h(j+1)=0.25d0*h(j)
11    continue
      pause 'too many steps in qromb'
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