      subroutine saveprofiles
      implicit none
      real*8 fn,polin,polin1
      external fn,polin,polin1
      integer i,k,j,inpt,inpt3,ipsy,ipsy1,klo,khi,ierr
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      real*8 drhodr(NRD),delta(NRD),ell(NRD),gamm(NRD),amy(NRD)
     &,afld(NRD),rha(NRD)
      integer itend0,kvv,ncheb,ntet,nnz,ispl,ilhdata,nspl,ncoef
      integer nmaxm,maxstep2,maxstep4,nr,ni1,ni2,niterat,im,ip
      integer ipri,iw,ismth,ismthalf,ismthout,inew,itor,ipoll
      real*8 coeffs(10),rm,r0,z0,cltn,zero,p_in
      real*8 rmin,rmax,sitet,cotet,xb1,yb1,xb2,yb2
      real*8 freq,xmi1,zi1,xmi2,zi2,dni2,xmi3,zi3,dni3
      real*8 rh,con,tem,temi,zeff,y2dn,y2tm,y2tmi,y2zeff
      real*8 cdl,cly,cgm,cmy,ynzm,pm,anz,apz,share
      real*8 energy,dra,dble,dsign
      real*8 cleft,cright,cdel,rbord,pchm0,pabs0,pgiter
      real*8 abtor,fpol,fdf,dfmy,hmin1
      real*8 zplus,zminus,ynzm0,pm0,rrange,eps,hdrob
      real*8 chebne,chebdne,chebddne
      real*8 xlog,zalfa,xmalfa,dn1,dn2,factor,xlogj

      common/a00/ xlog,zalfa,xmalfa,dn1,dn2,factor
      common /a0k/ cdl(10),cly(10),cgm(10),cmy(10),ncoef
      common /a0l3/ rh(501),y2dn(501),y2tm(501),y2tmi(501)
      common /a0l4/ con(501),tem(501),temi(501),nspl
      common /a0l5/ zeff(501),y2zeff(501)
      common /a0a1/ ynzm(1001),pm(1001),nmaxm(4)
      common /a0ab/ nr
      common /a0abcd/ ipri
      common /a0bcd/ eps
      common /a0bd/ rrange,hdrob
      common /a0cd/ rbord,maxstep2,maxstep4
      common /a0cdm/ hmin1
      common /a0ef1/ r0,z0,rm,cltn
      common/b0/ itend0
      common /cnew/ inew
      common/physpar/ freq,xmi1,zi1,xmi2,zi2,dni2,xmi3,zi3,dni3
      common/alfaspar/ energy,dra,kvv
      common/numpar/ cleft,cright,cdel,pchm0,pabs0,pgiter
      common/numpar1/ ni1,ni2,niterat
      common/optpar/ iw,ismth,ismthalf,ismthout
      common/grillpar/ zplus,zminus,ntet,nnz
      real*8 rh1,znak_tor,znak_pol
      common/left/ abtor,rh1,znak_tor,znak_pol,itor,ipoll
      real*8 ynzmp,pmp,ynzmm,pmm,plaunp,plaunm
      common/grillmp/ ynzmp(1001),pmp(1001),ynzmm(1001),pmm(1001)
      common/grillspektr/ ynzm0(1001),pm0(1001),ispl
     &,plaunp,plaunm,ip,im
      common/ne_cheb/chebne(50),chebdne(50),chebddne(50),ncheb
      real*8 efld(100),r,vmax,fvt,funmaxwell,fmaxw
      real*8 pme,pqe,pi4,c0,zff,zefff,fnr,fnrr
      real*8 pn,fn1,fn2,gst,dens,tmp,ft,vt,vclt
      integer i0
      parameter(i0=1002)
      real*8 vij,fij0,fij,dfij,dij,enorm,znak,fst
      common/lh/vij(i0,100),fij0(i0,100,2),fij(i0,100,2),dfij(i0,100,2)
     &,dij(i0,100,2),enorm(100),fst(100)
      real*8 calls
      common/firstcall/calls
      parameter(zero=0.d0,ipsy=5)
      save share

      ncoef=ipsy
      inpt=NA1          ! ASTRA radial grid number
      nspl=inpt
      p_in=dble(QLH)    ! input LH power, MW
      do i=1,inpt
!!!       rhj(i)=RHO(i)/ROC
       rh(i)=AMETR(i)/ABC
       rha(i)=RHO(i)/ABC  !/ABC instead of /ROC is not a mistake!
       delta(i)=(SHIF(1)-SHIF(i))/ABC  !FRTC Shafr. shift. defin.
       ell(i)=ELON(i)
       gamm(i)=rh(i)*TRIA(i)
!!variant       afld(i)=ULON(i)/RTOR/GP2
       afld(i)=UPL(i)/RTOR/GP2 !!variant
       !write(*,*)'afld=',afld(i),'Upl=',UPL(i)
       con(i)=NE(i)
       tem(i)=TE(i)
       temi(i)=TI(i)
       zeff(i)=ZEF(i)
      end do
      rh(inpt)=1.d0
      rh1=rh(1)          !saving the first ASTRA radial grid element
      rh(1)=zero         !shifting the first element to zero
      rha(1)=zero        !shifting the first element to zero
      delta(1)=zero      !putting delta(rh=0.)=0.
      gamm(1)=zero       !putting gamm(rh=0.)=0.

      abtor=1.d4*BTOR*RTOR/(RTOR+SHIF(1)) !B_tor_(magnetic axis), Gauss
      rm=1.d2*ABC                       !minor radius in mid-plane, cm
      r0=1.d2*(RTOR+SHIF(1))     !x-coordinate of the magnetic axis, cm
      z0=1.d2*UPDWN              !z-coordinate of the magnetic axis, cm

!!!!!!!!!!!!!! spline approximation of plasma profiles !!!!!!!!!!!!!!!!
      ipsy1=ipsy-1
!
cccc   shift as a function of "minor radius":
       call approx(rh,delta,inpt,polin1,ipsy1,coeffs)
       cdl(1)=zero
       do k=2,ipsy
        cdl(k)=coeffs(k-1)
       end do

cccc   triangularity as a function of "minor radius":
       call approx(rh,gamm,inpt,polin1,ipsy1,coeffs)
       cgm(1)=zero
       do k=2,ipsy
        cgm(k)=coeffs(k-1)
       end do

cccc   ellipticity as a function of "minor radius":
       call approx(rh,ell,inpt,polin,ipsy,cly)

cccc   "poloidal magnetic field":
       call diff(rh,rha,inpt,drhodr)
       do i=2,inpt
        amy(i)=1.d4*BTOR*MU(i)*rha(i)*drhodr(i)
       end do
       amy(1)=zero
!! amy=(btor/q)*rho*(drho/dr) is a function of "minor radius" r=rh(i).
!! Poloidal magnetic field: B_pol=amy(r)*sqrt(g22/g), where g is
!! determinant of 3D metric tensor and g22 is the (22) element of
!! the tensor, normalized on ABC^4 and ABC^2, correspondingly.
!!
!!  Polinomial approximation of the amy(r):
!!!      if(calls.eq.zero) then
       inpt3=inpt-3
       call approx(rh,amy,inpt3,polin1,ipsy1,coeffs)
       cmy(1)=zero
       do k=2,ipsy
        cmy(k)=coeffs(k-1)
       end do
!!!      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(calls.eq.zero) then
       call get_unit(ilhdata)
       if(ilhdata.eq.0) then
        write(*,*)'no free units up to 299 to read lhdata.dat'
        pause
        stop
       end if
 !!!      open(ilhdata,file='lhcd/gaus25.dat')
 !!!      open(ilhdata,file='lhcd/lhdataFT2_05p.dat')
       open(ilhdata,file='lhcd/ray_tracing.dat')
!!!      open(ilhdata,file='lhdatanew_Npar_2.0.dat')
!!!!!!!!!!!!!  read  physical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!
       read(ilhdata,*)
       read(ilhdata,*) freq
       read(ilhdata,*) xmi1
       read(ilhdata,*) zi1
       read(ilhdata,*) xmi2
       read(ilhdata,*) zi2
       read(ilhdata,*) dni2
       read(ilhdata,*) xmi3
       read(ilhdata,*) zi3
       read(ilhdata,*) dni3
!!!!!!!!!!!!  read parameters for alphas calculation !!!!!!!!!!!!!!!!
       read(ilhdata,*)
       read(ilhdata,*) itend0
       read(ilhdata,*) energy
       read(ilhdata,*) factor
       read(ilhdata,*) dra
       read(ilhdata,*) kvv
!!!!!!!!!!!!  read  numerical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       read(ilhdata,*)
       read(ilhdata,*) nr
       read(ilhdata,*) hmin1
       read(ilhdata,*) rrange
       read(ilhdata,*) eps
       read(ilhdata,*) hdrob
       read(ilhdata,*) cleft
       read(ilhdata,*) cright
       read(ilhdata,*) cdel
       read(ilhdata,*) rbord
       read(ilhdata,*) pchm0
       read(ilhdata,*) pabs0
       read(ilhdata,*) pgiter
       read(ilhdata,*) ni1
       read(ilhdata,*) ni2
       read(ilhdata,*) niterat
       read(ilhdata,*) nmaxm(1)
       read(ilhdata,*) nmaxm(2)
       read(ilhdata,*) nmaxm(3)
       read(ilhdata,*) nmaxm(4)
       read(ilhdata,*) maxstep2
       read(ilhdata,*) maxstep4
!!!!!!!!!!!!  read  options !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       read(ilhdata,*)
       read(ilhdata,*) ipri
       read(ilhdata,*) iw
       read(ilhdata,*) ismth
       read(ilhdata,*) ismthalf
       read(ilhdata,*) ismthout
       read(ilhdata,*) inew
       read(ilhdata,*) itor     !Btor direction in right-hand {drho,dteta,dfi}
       read(ilhdata,*) ipoll     !Bpol direction in right-hand {drho,dteta,dfi}

!!!!!!!!!!!!  read grill parameters and input LH spectrum !!!!!!!!!!!!
       read(ilhdata,*)
       read(ilhdata,*) zplus
       read(ilhdata,*) zminus
       read(ilhdata,*) ntet
       read(ilhdata,*) nnz
!!!!!!!!!!!!  read positive spectrum !!!!!!!!
       read(ilhdata,*)
       do i=1,10000
        read (ilhdata,*) anz,apz
        if(apz.eq.-88888.d0) then
         share=anz 
         go to 10
        end if
        ynzmp(i)=anz
        pmp(i)=apz
        ip=i
       end do
10     continue
       if(ip.gt.1001) then
        pause 'too many points in positive spectrum'
        stop
       end if
!!!!!!!!!!!!  read negative spectrum !!!!!!!!
       read(ilhdata,*)
       do i=1,10000
        read (ilhdata,*,end=11) ynzmm(i),pmm(i)
        im=i
       end do
11     continue
       close(ilhdata)
       if(im.gt.1001) then
        pause 'too many points in negative spectrum'
        stop
       end if
      end if
!
      plaunp=p_in*share !input power in positive spectrum
      plaunm=p_in*(1.d0-share) !input power in negative spectrum
!!!      write(*,*)'plaunp=',plaunp,' plaunm=',plaunm
!!!      pause
!
      call splne(rh,con,nspl,y2dn)
      call splne(rh,tem,nspl,y2tm)
      call splne(rh,zeff,nspl,y2zeff)
      call splne(rh,temi,nspl,y2tmi)
      if(inew.ne.0) then
       ncheb=20
       call chebft1(zero,1.d0,chebne,ncheb,fn)
       call chder(zero,1.d0,chebne,chebdne,ncheb)
       call chder(zero,1.d0,chebdne,chebddne,ncheb)
      end if
!
      znak_tor=dsign(1.d0,dble(itor))
      abtor=znak_tor*dabs(abtor)
      fpol=fdf(1.d0,cmy,ncoef,dfmy)
      znak_pol=dsign(1.d0,dble(ipoll))*dsign(1.d0,fpol)
      do i=1,ncoef
       cmy(i)=znak_pol*cmy(i)
      end do
!
      pi4=16.d0*datan(1.d0)
      pme=9.11e-28
      pqe=4.803e-10
      c0=dsqrt(pi4*pqe**2/pme)
      xlog=16.d0+dlog(16.d0)
!
      do j=1,nr
       r=dble(j)/dble(nr+1)
       call lock(rh,nspl,r,klo,khi,ierr)
       if(ierr.eq.1) then
        write(*,*)'lock error in saveprofiles, Efield'
        write(*,*)'j=',j,' rh(j)=',rh(j),' r=',r
        pause
        stop
       end if
       call linf(rh,afld,r,efld(j),klo,khi)
       if(inew.eq.0) then !vardens
        pn=fn1(r,fnr)
       else
        pn=fn2(r,fnr,fnrr)
       end if
       vt=fvt(r)
       tmp=ft(r)/0.16d-8  !Te,  KeV
       dens=pn/1.d+13     !10^13 cm^-3
       xlogj=dlog(5.1527d7*tmp*16.d0*dsqrt(tmp)/dsqrt(dens))
       enorm(j)=(3.835d0/xlogj)*efld(j)*tmp/dens
       enorm(j)=enorm(j)*5.d0/(5.d0+zefff(r))
!!       fst(j)=pn*xlogj*c0**4/pi4/vt**3
       fst(j)=((5.d0+zefff(r))/5.d0)*pn*xlogj*c0**4/pi4/vt**3
       end do
!
 !     open(96, file='lhcd/out/difsave.dat', position='append')
      if(calls.eq.zero) then
       do k=1,2
        znak=2.d0*dble(k)-3.d0
        do j=1,nr
         r=dble(j)/dble(nr+1)
         vclt=3.d10/fvt(r)
         vmax=2.d0*vclt
         do i=1,i0
          vij(i,j)=dble(i-1)*vmax/dble(i0-1)
          fij(i,j,k)=fmaxw(vij(i,j),znak*enorm(j),dfij(i,j,k))
!          write(96,*) i,fij(i,j,k)
!!!!          fij(i,j,k)=fmaxw(vij(i,j),zero,dfij(i,j,k))
          if(vij(i,j).ge.vclt) then
           fij(i,j,k)=zero
           dfij(i,j,k)=zero
          end if
          fij0(i,j,k)=fij(i,j,k)
          dij(i,j,k)=zero
         end do
        end do
       end do
      end if
      calls=1.d0
!      close(96)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function funmaxwell(v,dfunmaxwell)
      implicit none
      real*8 v,dfunmaxwell,arg,pi2sqrt
      parameter(pi2sqrt=2.506628274631d0)
      arg=-0.5d0*v**2
      funmaxwell=dexp(arg)/pi2sqrt
      dfunmaxwell=-v*funmaxwell
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      double precision function fmaxw(v,alfa2,dfmaxw)
      implicit none
      real*8 v,alfa2,dfmaxw
      real*8 arg,alfa,api,b,psiq,f,df,erfcc
      real*8 pi2sqrt,pisqrt,zero
      parameter(pi2sqrt=2.506628274631d0,pisqrt=1.77245385090552d0)
      parameter(zero=0.d0)
      if(alfa2.le.zero) then
       arg=-0.5d0*v**2*(1.d0-0.5d0*alfa2*v**2)
       fmaxw=dexp(arg)/pi2sqrt
       dfmaxw=-v*(1.d0-alfa2*v**2)*fmaxw
      else
       alfa=dsqrt(alfa2)
       api=2.d0*alfa*dexp(-0.25d0/alfa2)/pisqrt
       b=2.d0-erfcc(0.5d0/alfa)+api
       f=psiq(v,alfa2)
       fmaxw=(f+api)/b/pi2sqrt
       df=-v*((1.d0-alfa2*v**2)*f+api)
       dfmaxw=df/b/pi2sqrt
      end if
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      double precision function psiq(v,alfa2)
!!! psiq=exp(ksiV**2)*erfcc(ksiV)*exp(-0.25/alfa2)
      implicit none
      double precision v,alfa2,df
      double precision x,t,z,f,asymp,alfa,q,u
      double precision zero,zmax,pisqrt
      parameter(zero=0.d0,zmax=10.d0,pisqrt=1.77245385090552d0)
      alfa=dsqrt(alfa2)
      q=-0.25d0/alfa2
      x=0.5d0*(alfa*v**2-1.d0/alfa)
      z=abs(x)
      if(z.gt.zmax) then !asymptotics
       f=dexp(q)*(1.d0-0.5d0/z**2+0.75d0/z**4-15.d0/8.d0/z**6)/z/pisqrt
      else
       t=1.d0/(1.d0+0.5d0*z)
       f=t*exp(q-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+t*
     * (.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t*
     * (1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
      end if
      if(x.lt.zero) then
       u=-0.5d0*v**2+0.25d0*alfa2*v**4 !u=x**2-0.25d0/alfa2
       f=2.d0*dexp(u)-f
      end if
      psiq=f
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function erfcc(x)
      implicit none
      double precision erfcc,x
      double precision t,z,zero,zmax,pisqrt
      parameter(zero=0.d0,zmax=10.d0,pisqrt=1.77245385090552d0)
      z=abs(x)
      if(z.gt.zmax) then !asymptotics
       erfcc=(1.d0-0.5d0/z**2+0.75d0/z**4-15.d0/8.d0/z**6)/z/pisqrt
       erfcc=exp(-z*z)*erfcc
      else
       t=1.d0/(1.d0+0.5d0*z)
       erfcc=t*exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+t*
     * (.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t*
     * (1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
      end if
      if(x.lt.zero) erfcc=2.d0-erfcc
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

