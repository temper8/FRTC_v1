      subroutine lhcd2017(outpe)
cc******************************************************************
cc   outj(i)  = LH driven current density, MA/m^2
cc   outpe(i) =LH power density (Landau+coll.) deposited into electrons, MW/m^3
cc   outpec(i) = LH power density (collisions) deposited into electrons, MW/m^3
cc   outpef(i) = LH power dens. dep. into el. by fast wave, MW/m^3
cc   dndt(i)  = d^2Jr1/dt^2/E, MA/m^2/sec^2/(V/m), ~runaway d(el.density)/dt/E
cc   djdt(i)  = dJr2/dt, time drivative of runaway current Jr2, MA/m^2/sec
cc   outpa(i)  = LH power density deposited into alphas, MW/m^3
cc   outda(i)  = (Na/Ne) relative density of alpha-particles
cc******************************************************************
      implicit none
      integer i,k,iview,ipsy,ipsy1,inpt,iunit
      integer klo,khi,ierr,inpt2,ispectr
      real*8 abtor,rm,x0,z0,rh1,zero,p_in,pe_p,pe_m,c_p,c_m
      real*8 polin,polin1,vint,tcur
      common/testf/ tcur
      external polin,polin1
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      real*8 outpe(NRD)
      real*8,dimension(:),allocatable:: outpep,outpem
      real*8,dimension(:),allocatable:: con,tem,temi,azef,afld
      real*8,dimension(:),allocatable:: rh,rha,drhodr,delta,ell,gamm,amy
      real*8,dimension(:),allocatable:: cdl,cly,cgm,cmy,coeffs
      parameter(zero=0.d0,ipsy=5)
!
cc*********************************************************************
cc   ipsy = number of polinomial decomposition coefficients
cc           used for interpolation of Zakharov's moments.
cc*********************************************************************

cc*********************************************************************
cc    Co-ordinates used in ray-tracing:
cc         (x-x0)/rm=r*cos(teta)-delta-gamma*sin^2(teta) !sav2008 -gamma
cc         (z-z0)/rm=ell*r*sin(teta)
cc    Definitions:
cc    (x0,z0) - magnetic axis position, centimeters
cc    rm      - minor radius in mid-plane, cenrimeters
cc    r(rho_ASTRA),delta(r),gamma(r),ell(r) - dimensionless functions
cc    rho_ASTRA=sqrt(Phi_tor/GP/BTOR)
cc    Interval for r:  0.<= r <=1.
cc*********************************************************************

      tcur=time
      inpt=NA1          ! ASTRA radial grid number
      outpe=zero
      p_in=dble(QLH)    ! input LH power, MW
!!!      if(p_in.le.zero) return

      allocate(rh(inpt),rha(inpt),drhodr(inpt),con(inpt),tem(inpt))
      allocate(temi(inpt),azef(inpt))
      allocate(delta(inpt),ell(inpt),gamm(inpt),amy(inpt))
      allocate(cdl(ipsy),cly(ipsy),cgm(ipsy),cmy(ipsy),coeffs(ipsy))

      do i=1,inpt
       rh(i)=AMETR(i)/ABC
       rha(i)=RHO(i)/ABC  !/ABC instead of /ROC is not a mistake!
       delta(i)=(SHIF(1)-SHIF(i))/ABC  !FRTC Shafr. shift. defin.
       ell(i)=ELON(i)
       gamm(i)=rh(i)*TRIA(i)
       con(i)=dble(NE(i))
       tem(i)=dble(TE(i))
       temi(i)=dble(TI(i))
       azef(i)=dble(ZEF(i))
!!!variant       afld(i)=ULON(i)/RTOR/GP2
!!!       afld(i)=UPL(i)/RTOR/GP2 !!variant
      end do

      rh1=rh(1)          !saving the first ASTRA radial grid element
      rh(1)=zero         !shifting the first element to zero
      rha(1)=zero        !shifting the first element to zero
      delta(1)=zero      !putting delta(rh=0.)=0.
      gamm(1)=zero       !putting gamm(rh=0.)=0.

      abtor=1.d4*BTOR*RTOR/(RTOR+SHIF(1)) !B_tor_(magnetic axis), Gauss
      rm=1.d2*ABC                       !minor radius in mid-plane, cm
      x0=1.d2*(RTOR+SHIF(1))     !x-coordinate of the magnetic axis, cm
      z0=1.d2*UPDWN              !z-coordinate of the magnetic axis, cm

      ipsy1=ipsy-1

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
       inpt2=inpt-3
       call approx(rh,amy,inpt2,polin1,ipsy1,coeffs)
       cmy(1)=zero
       do k=2,ipsy
        cmy(k)=coeffs(k-1)
       end do

!!!!!!!!!!!!! starting ray-tracing !!!!!!!!!!!!!!!!!!!!!
      allocate(outpep(inpt),outpem(inpt))

!!positive spectrum:
      ispectr=1
      pe_p=zero
      outpep=zero
      call ourlhcd2017(ispectr,p_in,inpt,ipsy,rm,x0,z0,abtor,
     &  tem,con,temi,azef,rh,rh1,cdl,cly,cgm,cmy, !input
     &  outpep,pe_p)     !output
       if(pe_p.ne.zero) then
        c_p=vint(outpep,roc)
        if(c_p.ne.zero) then
         do i=1,inpt
          outpep(i)=pe_p*outpep(i)/c_p
         end do
        end if
       end if

!!negative spectrum:
      ispectr=-1
      pe_m=zero
      outpem=zero
      call ourlhcd2017(ispectr,p_in,inpt,ipsy,rm,x0,z0,abtor,
     &  tem,con,temi,azef,rh,rh1,cdl,cly,cgm,cmy, !input
     &  outpem,pe_m)     !output
       if(pe_m.ne.zero) then
        c_m=vint(outpem,roc)
        if(c_m.ne.zero) then
         do i=1,inpt
          outpem(i)=pe_m*outpem(i)/c_m
         end do
        end if
       end if
!
      do i=1,inpt
       outpe(i)=outpep(i)+outpem(i)
      end do
!
      deallocate(outpep,outpem)
      deallocate(rh,rha,drhodr,con,tem,temi,azef)
      deallocate(delta,ell,gamm,amy,cdl,cly,cgm,cmy,coeffs)
 !!     pause

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ourlhcd2017(ispec,p_in,ngrid,ncoe,arm,ar0,az0,btt
     & ,ate,xne,ati,azef,rhj,rh1,acdl,acly,acgm,acmy,  ! input data
     & outpe,pe_out) ! output data
      implicit real*8 (a-h,o-z)
      real*8 outpe,pe_out !,outpec,outpef,outpa,outda
      dimension outpe(*),rhj(*),ate(*),xne(*),azef(*),ati(*)
      dimension acdl(ncoe),acly(ncoe),acgm(ncoe),acmy(ncoe)
!!!!    maximum grid sizes:
!!!!    rho=100, v_par_electrons=100, v_perp_ions=50
      external obeom,ploshad,fn
      dimension galfa(50,100),vpmin(100),vcva(100)
     &,pd2(100),pd2a(100),pd2b(100),pdprev1(100),pdprev2(100)
     &,source(100),sour(100)
     &,rxx(102),pwe(102),wrk(102)
     &,ynzm0(1001),pm0(1001),yn2z(1001),powinp(1001)
      dimension vmid(100),vz1(100),vz2(100),ibeg(100),iend(100)
      parameter(mpnt=10000)
      common/refl/nrefj(mpnt)
      common /a0a1/ ynzm(1001),pm(1001),nmaxm(4)
      common /a0a2/ tet1,tet2
      common /a0a4/ plost,pnab
      common /a0ab/ nr
      common /a0abcd/ ipri
      common /a0bcd/ eps
      common /a0bcp/ tin
      common /a0bd/ rrange,hdrob
      common /a0befr/ pi,pi2
      common /a0cd/ rbord,maxstep2,maxstep4
      common /a0cdm/ hmin1
      common /a0ef1/ r0,z0,rm,cltn
      common /bcef/ ynz,ynpopq
      common /a0ef2/ btor,ww
      common /a0ef3/ xmi,c0,c1,cnye,cnyi,xsz,vt0
      common /a0gh/ pabs
      common /a0ghp/ vlf,vrt,dflf,dfrt
      common/plosh/ zv1(100,2),zv2(100,2),sk(100)
      common /a0i2/ vk(100),pchm,pme
      common /a0i3/ dql(101,100),pdl(100),vzmin(100),vzmax(100)
      common /a0i4/ fcoll(100),dens(100),eta(100)
      common /asou/ rsou(102),sou(102),npta
      common /a0i5/ vperp(50,100),cnstal,zza,zze,valfa,kv
      common /a0k/ cdl(10),cly(10),cgm(10),cmy(10),ncoef
      common /a0l3/ rh(501),y2dn(501),y2tm(501),y2tmi(501)
      common /a0l4/ con(501),tem(501),temi(501),nspl
      common /a0l5/ zeff(501),y2zeff(501)
      common/gridv/vgrid(101,100),dfundv(101,100),nvpt
      common/vvv1/dq1(101,100),dq2(101,100),pdc(100),pda(100),ppv1,ppv2
      common/findsigma/dncount(101,100)
      common /vvv2/ psum4
      common /vvv3/ pdfast(100)
      common /arr/ dgdu(50,100),kzero(100)
      common /alph/ dqi0(50,100)
      common/a00/ xlog,zalfa,xmalfa,dn1,dn2,factor
      common /ag/ inak,lenstor,lfree
      common/b0/ itend0
      common /maxrho/ rmx_n,rmx_t,rmx_z,rmx_ti
      common /cnew/ inew !est !sav2008
      common/ne_cheb/chebne(50),chebdne(50),chebddne(50),ncheb
      parameter(zero=0.d0, one=1.d0)
      parameter(cnst1=0.2965924106d-6)  ! cnst1=(m_e/m_p)**2, CGS
      parameter(cnst2=0.359680922d-35)  ! cnst2=(m_e/e)**2,  CGS
      integer i0,ispectr
      parameter(i0=1002)
      real*8 vij,fij0,fij,dfij,dij,enorm,fst,kofpar,timecof
      common/lh/vij(i0,100),fij0(i0,100,2),fij(i0,100,2),dfij(i0,100,2)
     &,dij(i0,100,2),enorm(100),fst(100)
      real*8 dijk,vrjnew
      common/t01/dijk(101,100,2),vrjnew(101,100,2),iptnew
      real*8,dimension(:),allocatable:: vvj,vdfj
      integer kpt1,kpt3
      parameter(kpt1=20,kpt3=20)
      double precision vrj(101),dj(101),djnew(1001)
      double precision dj2(101),d2j(101)
      common/testf/ tcur
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      if(p_in.eq.zero) then
       do j=1,100
        do k=1,2
         do i=1,i0
          dij(i,j,k)=zero
          dijk(i,j,k)=zero
         end do
        end do
       end do
       return
      end if
      ispectr=ispec
!
      pi=4.d0*datan(1.d0)
      pi2=2.d0*pi
      pi4=4.d0*pi
      piq=dsqrt(pi)

      lfree=1
      ncoef=ncoe
      nspl=ngrid
      rm=arm
      r0=ar0
      z0=az0
      btor=btt
      do i=1,ncoef
       cdl(i)=acdl(i)
       cly(i)=acly(i)
       cgm(i)=acgm(i)
       cmy(i)=acmy(i)
      end do
!
       call get_unit(iunit)
       if(iunit.eq.0) then
        write(*,*)'no free units up to 299'
        pause
        stop
       end if
 !!!      open(iunit,file='lhcd/lhdataFT2_05m.dat')
 !!!      open(iunit,file='lhcd/gaus25.dat')
       open(iunit,file='lhcd/ray_tracing.dat')
!!!!!!!!!!!!!  read  physical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       read(iunit,*)
       read(iunit,*) freq
       read(iunit,*) xmi1
       read(iunit,*) zi1
       read(iunit,*) xmi2
       read(iunit,*) zi2
       read(iunit,*) dni2
       read(iunit,*) xmi3
       read(iunit,*) zi3
       read(iunit,*) dni3
!!!!!!!!!!!!!  read parameters for alphas calculation !!!!!!!!!!!!!!!!!!!
       read(iunit,*)
       read(iunit,*) itend0
       read(iunit,*) energy
       read(iunit,*) factor
       read(iunit,*) dra
       read(iunit,*) kv

!!!!!!!!!!!!!  read  numerical parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       read(iunit,*)
       read(iunit,*) nr
       read(iunit,*) hmin1
       read(iunit,*) rrange
       read(iunit,*) eps
       read(iunit,*) hdrob
       read(iunit,*) cleft
       read(iunit,*) cright
       read(iunit,*) cdel
       read(iunit,*) rbord
       read(iunit,*) pchm
       read(iunit,*) pabs0
       read(iunit,*) pgiter
       read(iunit,*) ni1
       read(iunit,*) ni2
       read(iunit,*) niterat
       read(iunit,*) nmaxm(1)
       read(iunit,*) nmaxm(2)
       read(iunit,*) nmaxm(3)
       read(iunit,*) nmaxm(4)
       read(iunit,*) maxstep2
       read(iunit,*) maxstep4

!!!!!!!!!!!!!  read  options !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       read(iunit,*)
       read(iunit,*) ipri
       read(iunit,*) iw
       read(iunit,*) ismth
       read(iunit,*) ismthalf
       read(iunit,*) ismthout
       read(iunit,*) inew

       read(iunit,*) itor     !Btor direction in right-hand {drho,dteta,dfi}
       read(iunit,*) ipol     !Bpol direction in right-hand {drho,dteta,dfi}

       znak_tor=dsign(1.d0,dble(itor))
       btor=znak_tor*dabs(btor)
       fpol=fdf(1.d0,cmy,ncoef,dfmy)
       znak_pol=dsign(1.d0,dble(ipol))*dsign(1.d0,fpol)
       do i=1,ncoef
        cmy(i)=znak_pol*cmy(i)
       end do

!!!!!!!!!!!!!  read grill parameters and input LH spectrum !!!!!!!!!!!!
       read(iunit,*)
       read(iunit,*) zplus
       read(iunit,*) zminus
       read(iunit,*) ntet
       read(iunit,*) nnz
       read(iunit,*)
!
       if(ispectr.eq.1) then !read positive spectrum
        do i=1,10000
         read (iunit,*) anz,apz
         if(apz.eq.-88888.d0) then
          plaun=p_in*anz !input power in positive spectrum
          if(plaun.eq.zero) then
           dij(:,:,1)=zero
           dijk(:,:,1)=zero
           close(iunit)
           return
          end if
          go to 10
         end if
         ynzm0(i)=anz
         pm0(i)=apz
         i1=i
        end do
       else if(ispectr.eq.-1) then !read negative spectrum
        apz=zero
        do while(apz.ne.-88888.d0)
         read (iunit,*) anz,apz
        end do
        read(iunit,*)
        plaun=p_in*(1.d0-anz) !input power in negative spectrum
        if(plaun.eq.zero) then
         dij(:,:,2)=zero
         dijk(:,:,2)=zero
         close(iunit)
         return
        end if
        do i=1,10000
         read (iunit,*,end=10) ynzm0(i),pm0(i)
         i1=i
        end do
       else
        write(*,*)'wrong ispectr=',ispectr
        pause
        stop
       end if
10     ispl=i1
      close(iunit)
      if(ispl.gt.4001) stop 'too many points in spectrum'

!!!!!!!!!!!!! test !!!!!
      write(*,*)'time=',tcur
      if((tcur-0.030d0)*(tcur-0.036d0).lt.zero) then
       call get_unit(iunit)
       if(iunit.eq.0) then
        write(*,*)'#test: no free units up to 299'
        pause
        stop
       end if
       if(ispectr.eq.1) then
        open(iunit,file='test/fundfunp.dat')
         write(iunit,*)'time=',tcur
         write(iunit,*)
       else
        open(iunit,file='test/fundfunm.dat')
         write(iunit,*)'time=',tcur
         write(iunit,*)
       end if
       k=(3-ispectr)/2
       do j=1,nr
        do i=1,i0
         write(iunit,88) vij(i,j),fij(i,j,k),dfij(i,j,k)
     &                  ,fij0(i,j,k),dble(j)
        end do
        write(iunit,*)
       end do
       close(iunit)
      end if
88    format(1x,10(e14.7,1x))


!!!!!!!!!!!!! checking initial parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(kv.gt.50) kv=50
      if(nr.gt.100) nr=100
      if(ni1.eq.0) ni1=20
      if(ni2.eq.0) ni2=20
      if(ni1+ni2.gt.100) then
       ni1=60
       ni2=40
      end if
      if(nnz*ntet.gt.10000) then
       nnz=250
       ntet=40
       pause 'nnz and ntet changed, because nnz*ntet>10000'
      end if
!!!!!!!!!!!!!!! spline approximation of plasma profiles !!!!!!!!!!!!!!!!
      rhj(nspl)=1.d0
      do j=1,nspl
       rh(j)=rhj(j)
       con(j)=xne(j)
       tem(j)=ate(j)
       zeff(j)=azef(j)
       if(itend0.gt.0) then
        temi(j)=ati(j)
       end if
      end do

      call splne(rh,con,nspl,y2dn)
      call splne(rh,tem,nspl,y2tm)
      call splne(rh,zeff,nspl,y2zeff)
      if(itend0.gt.0) then
       call splne(rh,temi,nspl,y2tmi)
      end if
!
      if(inew.ne.0) then
       ncheb=20
       call chebft1(zero,1.d0,chebne,ncheb,fn)
       call chder(zero,1.d0,chebne,chebdne,ncheb)
       call chder(zero,1.d0,chebdne,chebddne,ncheb)
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      hr=1.d0/dble(nr+1)
      xly=fdf(one,cly,ncoef,xlyp)
      arg1=(zplus-z0)/(xly*rm)
      arg2=(zminus-z0)/(xly*rm)
      if(dabs(arg1).lt.1.d0) then
       tet1=dasin(arg1)      ! upper grill corner poloidal coordinate
      else
       tet1=0.5d0*pi         ! upper grill corner poloidal coordinate
      end if
      if(dabs(arg2).lt.1.d0) then
       tet2=dasin(arg2)      ! lower grill corner poloidal coordinate
      else
       tet2=-0.5d0*pi        ! lower grill corner poloidal coordinate
      end if
c---------------------------------------------------
c initial constants
c---------------------------------------------------
      talfa=3.5d0    ! alpha particles' birth energy, MeV
      zalfa=2.d0     ! alpha particles' electrical charge
      xmalfa=4.d0    ! alpha particles' atomic mass
      tin=1d-7
      clt=3.0d+10
      pme=9.11e-28
      pqe=4.803e-10
      xlog=16.d0+dlog(16.d0)
      c0=dsqrt(pi4*pqe**2/pme)
      c1=pqe/pme/clt
      xsgs=1d+13
      xwtt=1d-7
c------------------------------------------------------------
c calculate constants
c---------------------------------------
      dn1=1d0/(zi1+dni2*zi2+dni3*zi3)
      dn2=dni2*dn1
      dn3=dni3*dn1
      sss=zi1**2*dn1/xmi1+zi2**2*dn2/xmi2+zi3**2*dn3/xmi3
      xmi=1836.d0/sss
      cnstvc=(.75d0*piq*sss/1836.d0)**(1.d0/3.d0)
      ww=freq*pi2*1.0d+09
      cnye=xlog/pi4
      cnyi=dsqrt(2d0)/(3d0*piq) !%for Vt=sqrt(Te/m)
      vt0=fvt(zero)
!!!!!!!!      ptkev=ft(zero)/0.16d-8  !Te in keV
      cltn=clt/vt0
      xsz=clt/ww/rm
      ccur=pqe*vt0*0.333d-9
!!      ccurnr=pqe*pqe*0.333d-9/pme
      rrange=rrange*hr
      iw0=iw

      valfa=1.d9*dsqrt(1.91582d0*talfa/xmalfa)
cccc  valfa (cgs units) = birth velocity
      zza=cnst1*(zalfa/xmalfa/valfa)**2*(clt/valfa)**3/pi
      zze=cnst2*2.d9*freq
      cnstal=(dsqrt(cnst1)/xmalfa/pi)
     &        *(zalfa*vt0/valfa)**2*clt/valfa
      vpmax=dsqrt(energy/talfa)
cccc  "vpmax" in valfa velocity units !
c--------------------------------------------------------
c find volums and surfaces
c--------------------------------------------------------
      vk0=pi2*hr*rm**3
      sk0=hr*rm**2
      nrr=nr+2
      rxx(1)=zero
      rxx(nrr)=one
      do j=1,nr
       rxx(j+1)=hr*dble(j)
       vk(j)=vk0*gaussint(obeom,zero,pi2,rxx(j+1),eps)
       sk(j)=sk0*gaussint(ploshad,zero,pi2,rxx(j+1),eps)
      end do
!!!!!!!!!!!!!!!!!!!!!!!!
c--------------------------------------------------------
c  approximation of input LH spectrum
c--------------------------------------------------------
      call splne(ynzm0,pm0,ispl,yn2z)
      innz=100*ispl
      dxx=(ynzm0(ispl)-ynzm0(1))/innz
      xx2=ynzm0(1)
      yy2=pm0(1)
      pinp=0d0
       do i=1,innz
        xx1=xx2
        yy1=yy2
        xx2=xx1+dxx
        call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
        dpw=.5d0*(yy2+yy1)*(xx2-xx1)
        pinp=pinp+dpw
       end do

      dpower=pinp/dble(nnz)
      xx2=ynzm0(1)
      yy2=pm0(1)
      pwcurr=zero
      ptot=zero
      do i=1,nnz-1
       xx0=xx2
11     continue
        xx1=xx2
        yy1=yy2
        xx2=xx1+dxx
        call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
        dpw=.5d0*(yy2+yy1)*(xx2-xx1)
         if(pwcurr+dpw.gt.dpower) then
          xx2=xx1+dxx*(dpower-pwcurr)/dpw
          call splnt(ynzm0,pm0,yn2z,ispl,xx2,yy2,dynn)
          dpw=.5d0*(yy2+yy1)*(xx2-xx1)
          pwcurr=pwcurr+dpw
         else
          pwcurr=pwcurr+dpw
          go to 11
         end if
       ynzm(i)=.5d0*(xx2+xx0)
       pm(i)=pwcurr
       ptot=ptot+pwcurr
       pwcurr=zero
      end do
      ynzm(nnz)=.5d0*(ynzm0(ispl)+xx2)
      pm(nnz)=pinp-ptot
      pnorm=plaun*xsgs/(pinp*ntet)
      pmax=-1d+10
      do i=1,nnz
       call splnt(ynzm0,pm0,yn2z,ispl,ynzm(i),powinp(i),dynn)
       pm(i)=pm(i)*pnorm
       if (pm(i).gt.pmax) pmax=pm(i)
       ynzm(i)=dble(ispectr)*ynzm(i) !sav2009
      end do
c       call get_unit(iunit)
c       if(iunit.eq.0) then
c        write(*,*)'no free units up to 299'
c        pause
c        stop
c       end if
c       if(ispectr.eq.1) then
c        open(iunit,file='lhcd/out/used_spectrP.dat')
c       else if(ispectr.eq.-1) then
c        open(iunit,file='lhcd/out/used_spectrM.dat')
c       end if
c       do i=1,nnz
c        write(iunit,1008) ynzm(i),powinp(i)
c       end do
c       write(iunit,*)
c       close(iunit)
1008   format (1x,10(e14.7,3x))

      pabs=pabs0*pmax/1.d2
      ppv1=zero
      ppv2=zero
      pnab=zero
      plost=zero
      psum4=zero
      anb=zero
      fuspow=zero
      o_da=zero
c-------------------------------------------
c find velocity limits and initial dfdv
c--------------------------------------------
      ipt1=kpt1+1
      ipt2=ni1+ni2
      ipt=ipt1+ni1+ni2+kpt3
      if(ipt.gt.101) then
       write(*,*)'ipt >101'
       pause'stop program'
       stop
      end if
      nvpt=ipt

      do j=1,nr                  ! begin 'rho' cycle
       r=hr*dble(j)
!!!!sav2008       pn=fn(r)
!!       pn=fn1(r,fnr)
!!       pn=fn2(r,fnr,fnrr) !sav2008
       if(inew.eq.0) then !vardens
        pn=fn1(r,fnr)
       else
        pn=fn2(r,fnr,fnrr)
       end if
       dens(j)=pn
       vt=fvt(r)
       vto=vt/vt0
       wpq=c0**2*pn
       whe=dabs(btor)*c1
       v=wpq/ww**2
       u1=whe/ww
       u=u1**2
       e1=1d0-v*(1d0/xmi-1d0/u)
       e2=v/u1
       e3=v
       tmp=ft(r)/0.16d-8 !Te, keV
       cn1=dsqrt(50d0/tmp)  !sav2008
        if(itend0.gt.0) then
         eta(j)=1d0-v
         vcva(j)=cnstvc*vt*dsqrt(2d0)/valfa
         vpmin(j)=2.0d0*dsqrt(tmp/(-eta(j)))
222      continue
         dvperp=(vpmax-vpmin(j))/dble(kv-1)
            if(dvperp.le.zero) then
             vpmax=1.3d0*vpmax
             go to 222
            end if
         do k=1,kv
          vperp(k,j)=vpmin(j)+dble(k-1)*dvperp
         end do
         fcoll(j)=.5d-13*dens(j)*zalfa**2*xlog/xmalfa/tmp**1.5d0
         ddens=dn1*dens(j)
         tdens=dn2*dens(j)
         tt=fti(r)**0.33333d0    ! (ti, keV)^1/3
         source(j)=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
         anb=anb+source(j)*vk(j)
        end if
       cn2=dsqrt(dabs(e1))+e2/dsqrt(e3) !sav2008
c       vz1(j)=cleft*cltn/cn1  !Vpar/Vt0
c       vz2(j)=cright*cltn/cn2  !Vpar/Vt0
c       if(vz2(j).gt.0.9d0*cltn) vz2(j)=0.9d0*cltn
c       v1=vz1(j)/vto !Vpar/Vt(rho)
c       v2=vz2(j)/vto !Vpar/Vt(rho)

       vmax=cltn/vto
       v1=4.d0  !Vpar/Vt(rho)
       v2=10.d0 !cright*cltn/cn2 !10.d0 !Vpar/Vt(rho)
       if(v2.ge.vmax) v2=0.5d0*vmax
       if(v1.ge.v2) v1=v2-2.d0
       call gridvel(v1,v2,vmax,0.5d0,ni1,ni2,ipt1,kpt3,vrj)
       vz1(j)=v1*vto !Vpar/Vt0
       vz2(j)=v2*vto !Vpar/Vt0
       if(vz2(j).gt.0.9d0*cltn) vz2(j)=0.9d0*cltn
       do i=1,ipt
        vgrid(i,j)=vrj(i)*vto
       end do
      end do                     ! end 'rho' cycle

!!!!!!!!!read data !!!!!!!!!!!!
       allocate(vvj(i0),vdfj(i0))
       k=(3-ispectr)/2
       do j=1,nr
        r=hr*dble(j)
        vt=fvt(r)
        vto=vt/vt0
        do i=1,i0
         vvj(i)=vij(i,j)
         vdfj(i)=dfij(i,j,k) !=dfundv(i,j)*vto**2
        end do
        do i=1,ipt
         vrj(i)=vgrid(i,j)/vto   !Vpar/Vt
         call lock(vvj,i0,vrj(i),klo,khi,ierr)
         if(ierr.eq.1) then
          write(*,*)'lock error in read distribution function'
          write(*,*)'j=',j,'i0=',i0
          write(*,*)'vvj(1)=',vvj(1),' vvj(i0)=',vvj(i0)
          write(*,*)'i=',i,' vrj(i)=',vrj(i),' vmax=',cltn/vto
          write(*,*)
          pause'next key = stop'
          stop
         end if
         call linf(vvj,vdfj,vrj(i),dfout,klo,khi)
         dfundv(i,j)=dfout/vto**2
         if(dfundv(i,j).gt.zero) dfundv(i,j)=zero
        end do
       end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(itend0.gt.0) then  ! begin alpha-source renormalisation
       fuspow=anb*talfa*1.6022d-19
       anb0=anb
       anb=zero
       do j=1,nr
        r=hr*dble(j)
        if(r.le.dra) then
         tt=fti(zero)**0.33333d0
        else
         tt=fti(r-dra)**0.33333d0    ! (shifted ti, kev)^1/3
        end if
        ddens=dn1*dens(j)
        tdens=dn2*dens(j)
        sour(j)=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
        anb=anb+sour(j)*vk(j)
       end do
       aratio=anb0/anb
       rsou(1)=zero
       sou(1)=aratio*sour(1)
       do j=1,nr
        r=hr*dble(j)
        rsou(j+1)=r
        sou(j+1)=aratio*sour(j)
        if(j.eq.nr) sssour=source(j)
        source(j)=sou(j+1)
       end do
       npta=nr+2
       rsou(npta)=1.d0
       sou(npta)=aratio*sour(nr)
      end if
c------------------------------------
c set initial values of arrays
c------------------------------------
!
      dland=zero
      dcoll=zero
      perpn=zero
      dalf=zero
      vel=zero
      jrad=0
      iww=0
      tatai=zero
      xnpar=zero
      izz=zero
!
      pdl=zero
      pdc=zero
      pda=zero
      pdfast=zero
      pdprev1=zero
      pdprev2=zero
      tok=zero
      cur=zero
      pd2=zero
      pd2a=zero
      pd2b=zero
      dql=zero
      dq1=zero
      dq2=zero
      dncount=zero
      vzmin=cltn
      vzmax=-cltn
      kzero=kv
      nrefj=0

      if(itend0.gt.0) then
        do j=1,nr           ! begin 'rho' cycle
          do i=1,50
           dqi0(i,j)=zero
          end do
         call alphas(dqi0,vperp,j,kv,galfa)
        end do              ! end 'rho' cycle
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!sign of driven current in right coordinate system {dro,dteta,dfi}:
!!!!!curdir=+1.0 for current drive in positive direction "dfi"
!!!!!curdir=-1.0 for current drive in negative direction "dfi"
!!!!!spectrum Nz>0 is along dfi>0 and Nz<0 is along dfi<0
!!!!!it is also OK if Npar is used instead of Nz, but for Btor>0, that is along dfi>0
!!      curdir=-dble(ispectr)
!!!!!!!!!!!!!!! begin iterations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      q_rest=plaun
      iterat=0
      nvmin=1 !minimum counted events at a given radius rho
80    continue
      call manager(iterat,iw0,nnz,ntet)
c-----------------------------------------------
c  find achieved radial points jbeg-jend
c----------------------------------------------
      jbeg=1
      jend=0
      do j=1,nr
       nvach=0
       do i=1,nvpt
        nvach=nvach+dncount(i,j)
       end do
       if(nvach.lt.nvmin) then
        if(jend.eq.0) jbeg=jbeg+1
       else
        jend=j
       end if
      end do
      if(jend.eq.0.or.jbeg.ge.jend) then
       write(*,*)'failure: jbeg=',jbeg,' jend=',jend 
       pause
       stop
      end if

      do j=1,nr
       pdl(j)=pdl(j)*xwtt
       pdc(j)=pdc(j)*xwtt
       pda(j)=pda(j)*xwtt
       pdfast(j)=pdfast(j)*xwtt
       pwe(j+1)=(pdl(j)+pdc(j))/vk(j)
      end do
      pwe(1)=pwe(2)
      pwe(nr+2)=zero

!!   find nevyazka
!!----------------------------
      psum1=zero
      psum2=zero
      pchg=zero
      pchg1=zero
      pchg2=zero
       do j=1,nr
        dpw1=pdl(j)+pdc(j)
        dpw2=pda(j)
        psum1=psum1+dpw1**2
        psum2=psum2+dpw2**2
        pchg1=pchg1+(dpw1-pdprev1(j))**2
        pchg2=pchg2+(dpw2-pdprev2(j))**2
        pdprev1(j)=dpw1
        pdprev2(j)=dpw2
       end do
      if(psum1.ne.zero) pchg=pchg1/psum1 !sav2008
      if(psum2.ne.zero) pchg=pchg+pchg2/psum2
c----------------------------------------
c     calculate total current and power
c----------------------------------------
      cppl=zero
      cppc=zero
      cppa=zero
      cppf=zero
      do j=1,nr
       cppl=cppl+pdl(j)
       cppc=cppc+pdc(j)
       cppa=cppa+pda(j)
       cppf=cppf+pdfast(j)
      end do
      ol=cppl*1d-6
      oc=cppc*1d-6
      oa=cppa*1d-6
      of=cppf*1d-6
!!!!!!!!! prepare to the next iteration !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iterat=iterat+1
      q_abs=ol+oc+oa
      q_rest=pnab/xsgs
      q_cond=zero
      if(q_abs.ne.zero) q_cond=0.5d0*q_rest/q_abs
!!!      if(q_cond.le.pabs0.and.pchg.lt.pgiter)
      if(iterat.gt.5.and.q_cond.le.pabs0.and.pchg.lt.pgiter)
     & go to 110
      pchgprev=pchg
      if(ipri.gt.1) then
       write(*,*)
       write(*,*) 'ITERATION:'
       write(*,*) 'iteration=',iterat
       write (*,*) 'ispectr=',ispectr
       write(*,*) 'P_launched, MW=',plaun
       write(*,*) 'P_landau, MW=',ol
       write(*,*) 'P_coll, MW=',oc
       write(*,*) 'P_alph, MW=',oa
       write(*,*) 'Alphas power, MW=',fuspow
       write(*,*) 'P_fast (landau+coll), MW=',of
       write(*,*) 'P_lost, MW=',plost/xsgs
       write(*,*) 'P_not accounted, MW=',pnab/xsgs
       write(*,*) 'P_landau (strong absorption), MW=',ppv1/xsgs
       write(*,*) 'P_landau (weak absorption), MW=',ppv2/xsgs
       write(*,*) 'P_turns, MW=', psum4/xsgs
       write(*,*) 'efficiency, I(MA)/P(MW)=',oi/plaun !sav2008
       call integral(1,nspl,rh,con,avedens) !sav2010
       write (*,*) '<Ne>, m^-3=',avedens*1.d19,' R, m=',r0*1.d-2
       eta_eff=1.d17*avedens*r0*oi/plaun
       write (*,*) 'eta_eff=<Ne>*R*I/P, A/(W*m^2)=',eta_eff !sav2010
       write(*,*) 'nevyazka=', pchg
      end if
      if(iterat.le.niterat) then
c-------------------------------------------
c   recalculate f' for a new mesh
c-------------------------------------------
      k=(3-ispectr)/2
      do j=1,nr
       r=hr*dble(j)
       vt=fvt(r)
       vto=vt/vt0
       if(iterat.gt.0) then
        v1=dmin1(vzmin(j),vz1(j))
        v2=dmax1(vzmax(j),vz2(j))
       else
        v1=vzmin(j)
        v2=vzmax(j)
       end if
       vmax=cltn/vto
       vp1=v1/vto
       vp2=v2/vto
       call gridvel(vp1,vp2,vmax,cdel,ni1,ni2,ipt1,kpt3,vrj)
       do i=1,i0
        vvj(i)=vij(i,j)
        vdfj(i)=dfij(i,j,k) !=dfundv(i,j)*vto**2
       end do
       do i=1,ipt
        call lock(vvj,i0,vrj(i),klo,khi,ierr)
        if(ierr.eq.1) then
!!!         if(vrj(i).gt.vvj(i0)) exit
         write(*,*)'lock error in new v-mesh'
         write(*,*)'j=',j,' i0=',i0
         write(*,*)'vvj(1)=',vvj(1),' vvj(i0)=',vvj(i0)
         write(*,*)'i=',i,' vrj(i)=',vrj(i)
         write(*,*)
         pause'next key = stop'
         stop
        end if
        call linf(vvj,vdfj,vrj(i),dfout,klo,khi)
        vgrid(i,j)=vrj(i)*vto
        dfundv(i,j)=dfout/vto**2
        if(dfundv(i,j).gt.zero) dfundv(i,j)=zero
       end do
       vz1(j)=v1
       vz2(j)=v2
      end do
!!----------------------------
        ppv1=zero
        ppv2=zero
        psum4=zero
        pnab=zero
        plost=zero
        dql=zero
        dq1=zero
        dq2=zero
        dncount=zero
        vzmin=cltn
        vzmax=-cltn
        pdl=zero
        pdc=zero
        pda=zero
        pdfast=zero
        if(itend0.gt.0) then
         dqi0=zero
        end if
       goto 80
      end if
c------------------------------------------
c save results
c------------------------------------------
110   continue

      if(ipri.gt.0) then
       write (*,*)
       write (*,*) 'RAY-TRACING RESULTS:'
       write (*,*) 'ispectr=',ispectr
       write (*,*) 'P_launched, MW=',plaun
       write (*,*) 'P_landau, MW=',ol
       write (*,*) 'P_coll, MW=',oc
       write (*,*) 'P_alph, MW=',oa
       write (*,*) 'Alphas power, MW=',fuspow
       write (*,*) 'P_fast (landau+coll), MW=',of
       write (*,*) 'P_lost, MW=',plost/xsgs
       write (*,*) 'P_not accounted, MW=',pnab/xsgs
       write (*,*) '~~~~~~~~~~~~~~'
       write (*,*) 'P_landau (strong absorption), MW=',ppv1/xsgs
       write (*,*) 'P_landau (weak absorption), MW=',ppv2/xsgs
       write (*,*) 'P_turns, MW=', psum4/xsgs
       write (*,*) '~~~~~~~~~~~~~~'
       write (*,*) 'efficiency, I(MA)/P(MW)=',oi/plaun !sav2008
       call integral(1,nspl,rh,con,avedens) !sav2010
       write (*,*) '<Ne>, m^-3=',avedens*1.d19,' R, m=',r0*1.d-2
       eta_eff=1.d17*avedens*r0*oi/plaun
       write (*,*) 'eta_eff=<Ne>*R*I/P, A/(W*m^2)=',eta_eff !sav2010
       write (*,*) 'nevyazka=', pchg
       write (*,*) 'iterations=',iterat
       write (*,*) '-------------------------------------------'
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      k=(3-ispectr)/2
      do j=1,nr
       r=hr*dble(j)
       vt=fvt(r)
       vto=vt/vt0
       vmax=cltn/vto
       zff=(5d0+zefff(r))/5d0
       cnyfoc=zff*c0**4*cnye
       if(inew.eq.0) then !vardens
        pn=fn1(r,fnr)
       else
        pn=fn2(r,fnr,fnrr)
       end if
       dconst=vt0/(1.d-10*cnyfoc*pme*pn**2) !divided by 10^-10 here 
!!!!!!!!                         and multiplied by 10^-10 in dfind()
!!!old       dconst=vt0/(cnyfoc*pme*pn**2)
!!!        dj(i)=dql(i,j)*dconst*vto !D_normir
       do i=1,ipt
        vrj(i)=vgrid(i,j)/vto      !Vpar/Vt
        dj(i)=dql(i,j)*dconst*vto  !D_normir
      vrjnew(i,j,k)=vrj(i)
	dijk(i,j,k)=dj(i)
	iptnew=ipt 
      end do
       do i=1,i0
        if(vij(i,j).ge.vmax) then
         ddout=zero
        else
         call lock(vrj,ipt,vij(i,j),klo,khi,ierr)
         if(ierr.eq.1) then
          write(*,*)'lock error in output dql'
          write(*,*)'j=',j,'ipt=',ipt
          write(*,*)'vrj(1)=',vrj(1),' vrj(ipt)=',vrj(ipt)
          write(*,*)'i=',i,' v=',vij(i,j),' vmax=',vmax
          write(*,*)
          pause'next key = stop'
          stop
         end if
!!!         call linf(vrj,dj,vij(i,j),ddout,klo,khi)
!!         if(ddout.le.1.d0) ddout=zero
         ddout=dj(klo)
        end if
        dij(i,j,k)=ddout
       end do
       zv1(j,k)=vrj(ipt1)
       zv2(j,k)=vrj(ni1+ni2+ipt1)
      end do
        call view(tcur,1,nnz,ntet)  !writing trajectories into a file
!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(ismthout.ne.0) then
       do i=1,nrr
        wrk(i)=pwe(i)
       end do
       call fsmoth4(rxx,wrk,nrr,pwe)
      end if
!
      rh(1)=rh1
      if(rh(nspl).gt.1.d0) rh(nspl)=1.d0
      do j=1,nspl
       call lock2(rxx,nrr,rh(j),klo,khi,ierr)
       if(ierr.ne.0) then
        write(*,*)'lock2 error in profiles for ASTRA'
        write(*,*)'ierr=',ierr,' j=',j,' rh(j)=',rh(j)
        write(*,*)'rxx(1)=',rxx(1),' rxx(nrr)=',rxx(nrr)
        pause
       end if
       call linf(rxx,pwe,rh(j),fout,klo,khi)
       outpe(j)=fout
      end do
      pe_out=ol+oc
      rh(1)=zero
!
      deallocate(vvj,vdfj)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine manager(iterat,iw0,nnz,ntet)
      implicit real*8 (a-h,o-z)
      parameter(length=5000000, mpnt=10000)
      dimension dland(length),dcoll(length),perpn(length),dalf(length)
      dimension vel(length),jrad(length),iww(length),tetai(length)
      dimension xnpar(length),izz(length)
      common/agh/xnpar,vel,dland,dcoll,dalf,perpn,tetai,jrad,iww,izz
      dimension mbeg(mpnt),mend(mpnt),mbad(mpnt),rbeg(mpnt) !sav2008
      dimension tetbeg(mpnt),xnrbeg(mpnt),xmbeg(mpnt),yn3beg(mpnt)
      common/viewdat/mbeg,mend,mbad,rbeg,tetbeg,xnrbeg,xmbeg,yn3beg
      dimension iznzap(mpnt),iwzap(mpnt),irszap(mpnt)
      dimension rzap(mpnt),tetzap(mpnt),xmzap(mpnt),yn3zap(mpnt)
      common /a0a1/ ynzm(1001),pm(1001),nmaxm(4)
      common /a0a2/ tet1,tet2
      common /a0a4/ plost,pnab
      common /a0ab/ nr
      common /a0abcd/ ipri
      common /abc/ rzz,tetzz,xmzz,iznzz,iwzz,irszz
      common /abcd/ irs
      common /abcde/ izn,iw
      common /abcdg/ iabsorp
      common /abefo/ yn3
      common /a0befr/ pi,pi2
      common /acg/ pow
      common /a0gh/ pabs
      common /aef2/ icall1,icall2
      common /ag/ inak,lenstor,lfree
      common/refl/nrefj(mpnt)
      parameter(zero=0.d0)
      data mbad /mpnt*0/
      lenstor=length
      htet=zero
      hr=1.d0/dble(nr+1) !sav2008
      if(ntet.ne.1) htet=(tet2-tet1)/(ntet-1)
      irs=1
      iout=0
      mbeg(1)=1
      itr=0
      inak=0
      nnj=0
      do n_it=0,3
       nnj=nnj+nmaxm(n_it+1)
      end do
      maxref=nnj
      if(iterat.lt.3) nmax0=nmaxm(iterat+1)
      if(iterat.ge.3) nmax0=nmaxm(4)
       if(ipri.gt.1) then
        write(*,1001) iterat+1
        write(*,1002)
       end if
      ibad=0
c--------------------------------------
c begin outer loop on teta
c--------------------------------------
      do itet=1,ntet
        nref=0
        nbad1=0
        nbad2=0
        icall1=0
        icall2=0
        tetin=tet1+htet*(itet-1)
c--------------------------------------
c begin inner loop on nz
c--------------------------------------
        do inz=1,nnz
          itr=itr+1
!ipri          if(ipri.eq.4)  write(23,*)
          if(iterat.eq.0) then
c-----------------------------------------
c    find initial radius for a trajectory
c    on the 1th iteration
c-----------------------------------------
            yn=ynzm(inz) !sav2008, yn is introduced
            pow=pm(inz)
            irs=1
            iw=iw0
            rin=rini(xmin,tetin,xnr,yn,hr,ifail)
            if(ifail.eq.1) then
             if(ipri.gt.1) write (*,*) 'error: no roots'
             iabsorp=-1
             inak0=inak
             go to 10
            end if
            rbeg(itr)=rin !sav2008
            tetbeg(itr)=tetin !sav2008
            xnrbeg(itr)=xnr !sav2008
            xmbeg(itr)=xmin !sav2008
            yn3beg(itr)=yn3 !sav2008
          else
            if(mbad(itr).ne.0) then
              plost=plost+pm(inz)
              go to 31
            end if
            ib=mbeg(itr)
            ie=mend(itr)
            powexit=pm(inz)
            dltpow=pabs
            call dqliter(dltpow,ib,ie,hr,powexit,iout)
             if(nmax0.eq.0) then
              ib=mbeg(itr)
              ie=mend(itr)
              pow1=powexit
              pgamma=1.d0-pow1/pm(inz)
              powexit=pow1/pgamma
              dltpow=powexit-pow1+pabs
              call dqliter(dltpow,ib,ie,hr,powexit,iout)
              powexit=powexit-dltpow+pabs
              if(powexit.lt.zero) powexit=zero
              go to 30
             end if
	      if(iout.eq.0) then
              go to 30
            else
              tetin=tetzap(itr)
              xmin=xmzap(itr)
              rin=rzap(itr)
              yn3=yn3zap(itr)
              pow=powexit
              irs=irszap(itr)
              iw=iwzap(itr)
              izn=iznzap(itr)
              jrad(ie+1)=1
              dland(ie+1)=lfree
              inak=lfree-1
            end if
          end if
c---------------------------------------
c initial parameters for a trajectory
c---------------------------------------
          xm=xmin
          rstart=rin !sav2008
          tet=tetin
          nmax=nmax0
          iabsorp=0
          inak0=inak
c-------------------------------------
c call ray tracing
c-------------------------------------
          call traj(xm,tet,rstart,nmax,nb1,nb2,itet,inz) !sav2009
          nbad1=nbad1+nb1
          nbad2=nbad2+nb2
          nrefj(itr)=nrefj(itr)+nmax
          powexit=pow
          nref=nref+nmax
10        if(iabsorp.lt.0) then
c-------------------------------------
c    encounted problems
c-------------------------------------
            if (inak.eq.lenstor-1) then
              write (*,*) 'fix maximal length'
              nmax0=0
              do i=1,4
                nmaxm(i)=0
              end do
              iout=1
              goto20
            end if
             if (ipri.gt.1) then
              tetin0=tet1+htet*(itet-1)
              write (*,111) tetin0,ynzm(inz)
111           format(1x,'traj. with tet0=',f10.5,1x
     &                            ,', Ninput=',f10.5,1x,'failed')
             end if
            mbad(itr)=1
            plost=plost+pow
            inak=inak0
            mend(itr)=inak-1
            goto 30
          end if
c---------------------------------------
c remember end point of trajectory
c---------------------------------------
          rzap(itr)=rzz
          tetzap(itr)=tetzz
          xmzap(itr)=xmzz
          yn3zap(itr)=yn3
          iznzap(itr)=iznzz
          iwzap(itr)=iwzz
          irszap(itr)=irszz
          if(iterat.eq.0) then
            if(itr.gt.1) mbeg(itr)=mend(itr-1)+2
            mend(itr)=inak
            jrad(mend(itr)+1)=0
            lfree=mend(itr)+2
            inak=lfree-1
          end if
20        continue
          if(iout.ne.0) then
           dcoll(ie+1)=inak
           jrad(inak+1)=0
           lfree=inak+2
          end if
          if(nrefj(itr).gt.maxref.and.pow.gt.pabs) then !forced absorp
           if(pow.ge.pm(inz)) go to 30 !sav2008
           ib=mbeg(itr)
           ie=mend(itr)
           pow1=pow
           pgamma=1.d0-pow1/pm(inz)
           powexit=pow1/pgamma
           dltpow=powexit-pow1+pabs
           call dqliter(dltpow,ib,ie,hr,powexit,iout)
           powexit=powexit-dltpow+pabs
           if(powexit.lt.zero) powexit=zero
          end if
30       continue
         pnab=pnab+powexit
31       continue
        end do
        if(ipri.gt.1)
     &  write(*,1003)itet,icall1,icall2,nref,lfree-1,nbad1,nbad2
      end do
1001  format (30x,i4,' iteration')
1002  format (6x,'n',5x,'call2',6x,'call4',6x,'nrefl',4x,'last',
     *5x,'bad2',5x,'bad4')
1003  format (3x,i4,2(1x,i10),2x,i7,2x,i8,2(1x,i7),2(2x,i7))
1004  format(1x,i8)
1005  format(1x,i5)
1006  format (e14.7)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dql1 !sav2008
      implicit real*8 (a-h,o-z)
      parameter(length=5000000)
      real*8 radth,fvt
      dimension dland(length),dcoll(length),perpn(length),dalf(length)
      dimension vel(length),jrad(length),iww(length),tetai(length)
      dimension xnpar(length),izz(length)
      common/agh/xnpar,vel,dland,dcoll,dalf,perpn,tetai,jrad,iww,izz
      dimension an1(length),an2(length)
      common /xn1xn2/ an1,an2
      common /a0gh/ pabs
      common /vth/ vthc(length),poloidn(length)
      common /a0ghp/ vlf,vrt,dflf,dfrt
      common /abcdg/ iabsorp
      common /acg/ pow
      common /ag/ inak,lenstor,lfree
      common /bcg/ hrad
      common /bg/ im4
      common /ceg/ ipow,jfoundr
      common /a0ab/ nr
      common /eg1/ vfound,ifound
      common /eg2/ pdec1,pdec2,pdec3,pdecv,pdecal,dfdv,icf1,icf2
      common /eg3/ cf1,cf2,cf3,cf4,cf5,cf6
      common /dg/ pintld4,pintcl4,pintal4
      common /a0i5/ vperp(50,100),cnstal,zza,zze,valfa,kv
      common/b0/ itend0
      parameter(clt=3.d10,zero=0.d0)
      powpr=pow
      iabsorp=0
      hdis=hrad
      vz=vfound
      i=ifound
      if(i.eq.0) i=1
      j=jfoundr
      refr=cf1
      tet_i=cf2
      npoloid=cf6
      xparn=cf3
      xan1=cf4
      xan2=cf5
      ifast=icf1
      idir=icf2
      dek3=zero
      dfsr=(vlf*dflf+vrt*dfrt)/2d0*(vrt-vlf)
      vsr=(vrt+vlf)*(vrt-vlf)/2d0
c--------------------------------------
c   find power
c--------------------------------------
      if(im4.eq.1) then
!!       pintld=-pintld4*dfdv
       pintld=dabs(pintld4*dfdv)
       pintcl=dabs(pintcl4)
       if(itend0.gt.0) then
        argum=clt/(refr*valfa)
        dek3=zatukh(argum,j,vperp,kv)
       end if
       pintal=dabs(pintal4*dek3)
       dcv=pintld4/vsr
      else
       pintld=dabs(pdec1*hdis)
       pintcl=dabs(pdec2*hdis)
       pintal=dabs(pdec3*hdis)
       dcv=pdecv*hdis/vsr
      end if
      if(pabs.ne.zero) then
       powd=pow*dexp(-2d0*pintld)
       powccc=dexp(-2d0*pintcl)
       powcol=powd*powccc
       powal=powcol*dexp(-2d0*pintal)
       pow=powal
      end if
      if(pow.le.pabs) iabsorp=1
      pil=pintld
      pic=pintcl
      pia=pintal
      call dfind(j,i,vz,powpr,pil,pic,pia,dfsr,dcv
     &                         ,refr,vlf,vrt,ifast)
c-----------------------------------
c      memorize trajectory
c----------------------------------
      inak=inak+1
      if(inak.eq.lenstor) then
       write(*,*)'storage capacity exceeded !'
       iabsorp=-1
       inak=lenstor-1
       return
      end if
      vel(inak)=vz
      perpn(inak)=refr
      poloidn(inak)=npoloid
      tetai(inak)=tet_i
      radth=dble(j)/dble(31)
      vthc(inak)=3.d10/fvt(radth)
      iww(inak)=ifast
      izz(inak)=idir
      xnpar(inak)=xparn
      an1(inak)=xan1
      an2(inak)=xan2
      if(im4.eq.1) then
       jrad(inak)=-j
       dland(inak)=pintld4
       dcoll(inak)=pintcl4
       dalf(inak)=pintal4
       im4=0
       return
      end if
      jrad(inak)=j
      dland(inak)=pdecv
      dalf(inak)=pdecal
      if(ipow.ne.1) dcoll(inak)=powccc
      if(ipow.eq.1) dcoll(inak)=1d0
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dqliter(dltpow,ib,ie,h,powexit,iout) !sav2008
      implicit real*8 (a-h,o-z)
      parameter(length=5000000)
      dimension dland(length),dcoll(length),perpn(length),dalf(length)
      dimension vel(length),jrad(length),iww(length),tetai(length)
      dimension xnpar(length),izz(length)
      common/agh/xnpar,vel,dland,dcoll,dalf,perpn,tetai,jrad,iww,izz
      dimension an1(length),an2(length)
      common /xn1xn2/ an1,an2
      common /a0ghp/ vlf,vrt,dflf,dfrt
      common /a0i5/ vperp(50,100),cnstal,zza,zze,valfa,kv
      common /vvv2/ psum4
      common/b0/ itend0
      parameter(clt=3.d10,zero=0.d0)
      pow=powexit
      pdec1=zero
      pdec1z=zero
      pdec3=zero
      pdec3z=zero
      pdecv=zero
      pintld=zero
      pintal=zero
10    continue
      iout=0
      do i=ib,ie
c-----------------------------------
c restore memorized decrements and
c integrate power equation
c------------------------------------
       v=vel(i)
       jr=jrad(i)
       refr=perpn(i)
       ifast=iww(i)
       dek3=zero
       if(itend0.gt.0) then
        argum=clt/(refr*valfa)
        dek3=zatukh(argum,abs(jr),vperp,kv)
       end if
!!!old variant
!!!       call raspr(v,abs(jr),iv,df)
!!!       if(iv.eq.0) iv=1
!!!!!!!!!!!!!!!!!!!!!!!!!!
       call distr(v,abs(jr),iv,df)
!!       dfsr=v*df*(vrt-vlf)
!!       vsr=v*(vrt-vlf)
       dfsr=(vlf*dflf+vrt*dfrt)/2d0*(vrt-vlf) !sav2008
       vsr=(vrt+vlf)*(vrt-vlf)/2d0 !sav2008
       if(jr.lt.0) then !case of turn
        jr=-jr
!variant        pintld=-dland(i)*df
!!        pintld=-dland(i)*(dflf+dfrt)/2d0
        pintld=dabs(dland(i)*(dflf+dfrt)/2d0)
        pdec2=dexp(-2d0*dcoll(i))
        pintal=dabs(dalf(i)*dek3)
        pcurr=pdec2*dexp(-2d0*pintld-2d0*pintal)
        psum4=psum4+pow*(1d0-pcurr)
        dcv=dland(i)/vsr
       else
        pdec2=dcoll(i)
        pdecv=dland(i)
!!        pdec1=-pdecv*df
        pdec1=dabs(pdecv*df)
        pdec3=dabs(dalf(i)*dek3)
        pintld=(pdec1+pdec1z)/2d0*h
        pintal=(pdec3+pdec3z)/2d0*h
        pdec1z=pdec1
        pdec3z=pdec3
        dcv=pdecv*h/vsr
       end if
       powpr=pow
       if(dltpow.ne.zero) then
        powd=pow*dexp(-2d0*pintld)
        powcol=powd*pdec2
        powal=powcol*dexp(-2d0*pintal)
        pow=powal
       end if
       pil=pintld
       pic=.5d0*dabs(dlog(pdec2))
       pia=pintal
       call dfind(jr,iv,v,powpr,pil,pic,pia,dfsr,dcv
     &                           ,refr,vlf,vrt,ifast)
       if(pow.lt.dltpow) then
        powexit=pow
        return
       end if
      end do
      jchek=jrad(ie+1)
c-------------------------------------------
c  check whether trajectory has continuation
c---------------------------------------------
      if(jchek.eq.0) then
       iout=1
       powexit=pow
       return
      else
       ib=idnint(dland(ie+1))
       ie=idnint(dcoll(ie+1))
       goto 10
      end if
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dfind(j,i,v,powpr,pil,pic,pia,df,decv
     &                             ,refr,vlf,vrt,ifast)
      implicit real*8 (a-h,o-z)
      common /a0i2/ vk(100),pchm,pme
      common /a0i3/ dql(101,100),pdl(100),vzmin(100),vzmax(100)
      common /a0i4/ fcoll(100),dens(100),eta(100)
      common /a0i5/ vperp(50,100),cnstal,zza,zze,valfa,kv
      common/vvv1/dq1(101,100),dq2(101,100),pdc(100),pda(100),ppv1,ppv2
      common /vvv3/ pdfast(100)
      common /alph/ dqi0(50,100)
      common/b0/ itend0
      common /a0ef1/ r0,z0,rm,cltn
      common /a0befr/ pi,pi2
      common/findsigma/dncount(101,100)
      parameter(zero=0.d0,clt=3.0d10,tiny=1.d-100)

      if(v.gt.cltn) return
      if(pil.gt.zero) then
       if(v.lt.vzmin(j)) vzmin(j)=v
       if(v.gt.vzmax(j)) vzmax(j)=v
      end if
      pchgl=zero
      pchgc=zero
      pchga=zero
      denom=pil+pic+pia
      powlandau=1.d0-dexp(-2.d0*pil)
      powdamped=1.d0-dexp(-2.d0*denom)
      domin=powpr*powdamped
      if(denom.ne.zero) then
!!       pchgl=powpr*(1.d0-dexp(-2d0*pil))
!!       pchgc=powpr*dexp(-2d0*pil)*dabs(-2d0*pic)
!!       pchga=powpr*dexp(-2d0*pil)*dabs(-2d0*pia)
       fff=domin/denom
       pchgl=dabs(pil*fff)
       pchgc=dabs(pic*fff)
       pchga=dabs(pia*fff)
      end if
      dd=zero
      if(pil.eq.zero) go to 1 !no Landau absorption
      if(powlandau.gt.pchm) then !strong absorption
       ppv1=ppv1+pchgl
       if(dabs(df).gt.tiny) then
        dd=dabs(-pchgl/vk(j)/(df*1.d10))
        dncount(i,j)=dncount(i,j)+1.d0
       else
        dd=zero
       end if
       dq1(i,j)=dq1(i,j)+dd
      else  ! weak absorption
       ppv2=ppv2+pchgl
       dd=dabs(2.d0*decv*powpr*1.d-10/vk(j))
       dncount(i,j)=dncount(i,j)+1.d0
       dq2(i,j)=dq2(i,j)+dd
      end if

1     continue
      dql(i,j)=dql(i,j)+dd
      pdl(j)=pdl(j)+pchgl
      pdc(j)=pdc(j)+pchgc
      pda(j)=pda(j)+pchga

      if(ifast.eq.-1) pdfast(j)=pdfast(j)+pchgl+pchgc+pchga
      if(itend0.gt.0) then
        parn=cltn/v
        dvz=vrt-vlf
        dnpar=cltn*dvz/v**2
        weight=(refr**2-eta(j))**2/(refr**2*parn**3)
!!!        adde=zze*(dd/dens(j))*weight
!!!        e2perp(i,j)=e2perp(i,j)+adde
        addd=zza*(dd/dens(j))*weight/fcoll(j)/refr**3
        arg=clt/(refr*valfa)
         do k=1,kv
           if(vperp(k,j).gt.arg) then
            hevis=dsqrt((vperp(k,j)-arg)*(vperp(k,j)+arg))
            adda=addd*hevis
            dqi0(k,j)=dqi0(k,j)+adda*dnpar
           end if
         end do
      end if
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine traj(xm0,tet0,xbeg,nmax,nb1,nb2,nomth,nomnz) !sav2009
      implicit real*8 (a-h,o-z)
      external extd4
      dimension ystart(2),yy(4)
      common /a0ab/ nr
      common /a0abcd/ ipri
      common /a0bcd/ eps
      common /a0bcp/ tin
      common /a0bd/ rrange,hdrob
      common /a0befr/ pi,pi2
      common /a0ef1/ r0,z0,rm,cltn
 !!     common /a0k/ cdl(10),cly(10),cgm(10),cmy(10),ncoef
      common /abc/ rzz,tetzz,xmzz,iznzz,iwzz,irszz
      common /abcd/ irs
      common /abcde/ izn,iw
      common /abcdg/ iabsorp
      common /bcg/ hrad
      common /bcef/ ynz,ynpopq
      common /bdeo/ ivar
      common /be1/ xnr1,xnr2,xnr3,xnr4
      common /be2/ ider
      common /beo/ iroot
      common /bg/ im4
      common /a0k/ cdl(10),cly(10),cgm(10),cmy(10),ncoef
      integer nomth,nomnz
      parameter (pgdop=0.02d0,hmin=0.d-7) !sav2008, old hmin=1.d-7
      eps0=eps
      rrange0=rrange
      hdrob0=hdrob
      iroot=1
      nrefl=0
      ider=1
      im4=0
      nb1=0
      nb2=0
      irep=0
      tet=tet0
      xm=xm0
      hr=1.d0/dble(nr+1) !sav2008
      hrad=hr
c---------------------------------------
c find saving point and define
c parameters, depending on direction
c---------------------------------------

10    irf1=idnint(xbeg/hr)
      if (dabs(irf1*hr-xbeg).lt.tin)  then
        xsav=hr*irf1
      else
        irf=int(xbeg/hr)
        if (irs.eq.1)  xsav=hr*irf
        if (irs.eq.-1) xsav=hr*(irf+1)
      end if
      xend=.5d0-.5d0*irs+tin*irs
      if (ipri.gt.2) write (*,*) 'xbeg-xend',xbeg,xend
      hsav=-hr*irs
      h1=hsav
c---------------------------------------
c solve eqs. starting from xbeg
c---------------------------------------
      ystart(1)=tet
      ystart(2)=xm
      call driver2(ystart,xbeg,xend,xsav,hmin,h1)
      tet=ystart(1)
      xm=ystart(2)
      ib2=0
         rnew2=ystart(3)
         cotet=dcos(tet)
         sitet=dsin(tet)
         xdl=fdf(rnew2,cdl,ncoef,xdlp)
         xly=fdf(rnew2,cly,ncoef,xlyp)
         xgm=fdf(rnew2,cgm,ncoef,xgmp)
         xx=-xdl+rnew2*cotet-xgm*sitet**2
         zz=rnew2*xly*sitet
         xxx=(r0+rm*xx)/1d2
         zzz=(z0+rm*zz)/1d2 
!      open(33,file='lhcd/out/dots.dat',position="append")
!      write(33,*)xxx, zzz, nomth, nomnz
!      close(33)
c---------------------------------------
c absorption
c---------------------------------------
      if(iabsorp.ne.0) then
        if(ipri.gt.2) write (*,*)'in traj() iabsorp=',iabsorp
        nmax=nrefl
        return
      end if
      if (xend.eq.xbeg) nb1=nb1+1
!sav2008 20    continue

c--------------------------------------------------------
c  pass turning point
c----------------------------------------------------------
      irs0=irs
      ider=0
      call disp2(xend,xm,tet,xnr,prt,prm)
      ider=1
      ynz0=ynz
40    yy(1)=tet
      yy(2)=xm
      yy(3)=xend
      yy(4)=xnr
      x1=0d0
      x2=1d+10
      rexi=xend
      call driver4(yy,x1,x2,rexi,hmin,extd4)
      if(iabsorp.eq.-1) return !failed to turn
      tetnew=yy(1)
      xmnew=yy(2)
      rnew=yy(3)
      xnrnew=yy(4)
         cotet=dcos(tetnew)
         sitet=dsin(tetnew)
         xdl=fdf(rnew,cdl,ncoef,xdlp)
         xly=fdf(rnew,cly,ncoef,xlyp)
         xgm=fdf(rnew,cgm,ncoef,xgmp)
         xx=-xdl+rnew*cotet-xgm*sitet**2
         zz=rnew*xly*sitet
         xxx=(r0+rm*xx)/1d2
         zzz=(z0+rm*zz)/1d2 
!      open(33,file='lhcd/out/dots.dat',position="append")
!      write(33,*)xxx, zzz, nomth, nomnz
!      close(33)
      if(ipri.gt.2) write (*,*) 'from r=',rexi,'to r=',rnew

c---------------------------------------
c find mode
c---------------------------------------
      iroot=3
      ider=0
      xnrv=xnrnew
      call disp2(rnew,xmnew,tetnew,xnrv,prt,prm)
      ider=1
      iroot=1
!ipric      if (ipri.gt.2) then
!ipric       write (*,*)'nr check, r=',rnew,' tet=',tetnew
!ipric       write (*,*)'iw=',iw,' izn=',izn
!ipric       write (*,*) xnrnew,xnr1
!ipric       write (*,*) xnr2,xnr3,xnr4
!ipric       pause
!ipric      end if
      pg1=dabs(xnrnew-xnr1)
      pg2=dabs(xnrnew-xnr2)
      pg3=dabs(xnrnew-xnr3)
      pg4=dabs(xnrnew-xnr4)
      pg=dmin1(pg1,pg2,pg3,pg4)
      if(dabs(pg/xnrnew).gt.pgdop) then
c---------------------------------------------
c bad accuracy, continue with 4 equations
c--------------------------------------------
        ib2=ib2+1
        nb2=nb2+1
        if (ib2.gt.4) then
          if (ipri.gt.1) write (*,*) 'error: cant leave 4 eqs'
          iabsorp=-1
          return
        end if
        eps=eps/5d0
        rrange=rrange*2d0
        hdrob=hdrob*2d0
        goto40
      end if
c-------------------------------------
c          change wave type
c-------------------------------------
      if (pg.ne.pg1) then
        if (pg.eq.pg2) izn=-izn
        if (pg.eq.pg3) iw=-iw
        if (pg.eq.pg4) iw=-iw
        if (pg.eq.pg4) izn=-izn
      end if
      if (irs0.ne.irs) nrefl=nrefl+1
      xbeg=rnew
      tet=tetnew
      xm=xmnew
      im4=1
      eps=eps0
      rrange=rrange0
      hdrob=hdrob0
      if(nrefl.lt.nmax) go to 10
      rzz=xbeg
      tetzz=tet
      xmzz=xm
      iznzz=izn
      iwzz=iw
      irszz=irs
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine driver2(ystart,x1,x2,xsav,hmin,h1) !sav2008
      implicit real*8 (a-h,o-z)
      external extd2
      common /a0abcd/ ipri
      common /a0ab/ nr
      common /a0bcd/ eps
      common /a0bcp/ tin
      common /a0befr/ pi,pi2
      common /a0cd/ rbord,maxstep2,maxstep4
      common /a0cdm/ hmin1
      common /a0ef1/ r0,z0,rm,cltn
      common /abc/ rzz,tetzz,xmzz,iznzz,iwzz,irszz
      common /abcd/ irs
      common /abcde/ izn,iw
      common /abcdg/ iabsorp
      common /bcef/ ynz,ynpopq
      common /bcg/ hrad
      common /cefn/ iconv,irefl
      common /ceg/ ipow,jfoundr
      common /cmn/ ind
      dimension ystart(2)
      parameter(zero=0.d0,nvar=2)
      dimension yscal(nvar),y(nvar),dydx(nvar),yold(nvar),dyold(nvar)
      x=x1
      h=dsign(h1,x2-x1)
      ind=0
      ipow=-1
      xold=x
      hsav=hrad*irs
      hdid=zero
      do i=1,nvar
       y(i)=ystart(i)
       yold(i)=y(i)
      end do
c-----------------------------
c            start moving
c-----------------------------
      do nstp=1,maxstep2
c---------------------------------------
c netpoint control
c---------------------------------------
       dstsav=dabs(x-xsav)
       if(dstsav.lt.tin) then
        ipow=ipow+2
        jfoundr=idnint(x/hrad)
        if(jfoundr.le.0) jfoundr=1
        if(jfoundr.gt.nr) jfoundr=nr
       end if
       call extd2(x,y,dydx)
       irep=0
       if(iconv+irefl.ne.0) then
!!       if(iconv+irefl.ne.0.or.ynz.lt.0.d0) then
c---------------------------------------------
c made step to nontransparent zone-return back
c----------------------------------------------
        x=xold
        do ii=1,nvar
         y(ii)=yold(ii)
         dydx(ii)=dyold(ii)
        end do
        irep=1
        h=hdid/2
        hdid=h
        ipow=0
        if(dabs(h).lt.hmin1) then
         ind=3
         go to 20
        end if
        ynz=ynz0
        go to 10
       end if
c--------------------------------------
c memorize step data
c--------------------------------------
       xold=x
       do i=1,nvar
        dyd=dabs(dydx(i))
        yscal(i)=dabs(y(i))+dabs(h*dyd)+1.d-30/(1.d0+dyd)
        yold(i)=y(i)
        dyold(i)=dydx(i)
       end do
       ynz0=ynz
       if(ipow.gt.0) then !integrate power equation
        call dql1
        if(iabsorp.eq.1) then !absorption
         rzz=x
         tetzz=y(1)
         xmzz=y(2)
         iznzz=izn
         iwzz=iw
         irszz=irs
         return
        end if
        if(iabsorp.eq.-1) return !problem
        ipow=0
        xsav=xsav-hsav
       end if
c--------------------------------------
c choose step size
c--------------------------------------
       dst3=(x-xsav)*(x+h-xsav)
       if(dst3.lt.zero.and.irep.eq.0) h=xsav-x
       if(x.gt.rbord.and.h.gt.zero) then
        ind=2
        go to 20
       end if
10     dst1=(x-rbord)*(x+h-rbord)
       dst2=x*(x+h)
       if((dst1.lt.zero.and.irs.eq.-1).or.dst2.lt.zero) then
        h=h/2.d0
        if(dabs(h).lt.hmin1) then
         ind=4
         go to 20
        end if
        go to 10
       end if
c--------------------------------------
c find solution at x=x+hdid
c---------------------------------------
       ynz0=ynz
       call difeq(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,extd2)
20     continue
       if(ind.ne.0) then !exit
        xsav=xsav+hsav
        x2=x
        do i=1,nvar
         ystart(i)=y(i)
        end do
        ynz=ynz0
        return
       end if
c---------------------------------------
       if(dabs(hnext).lt.hmin) then
        if(ipri.gt.1) write(*,*) 'exit driver2: step is too small'
        go to 40
       end if
       h=hnext
      end do
c---------------------------------------
      if (ipri.gt.1) write (*,*) 'error in driver2: too many steps'
40    iabsorp=-1
      return
1001  format (10(e14.7,1x))
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine driver4(ystart,x1,x2,rexi,hmin,derivs)
      implicit real*8 (a-h,o-z)
      external derivs
      common /a0befr/ pi,pi2
      common /a0abcd/ ipri
      common /a0bcd/ eps
      common /a0cdm/ hmin1
      common /a0bd/ rrange,hdrob
      common /a0cd/ rbord,maxstep2,maxstep4
      common /abcd/ irs
      common /abcde/ izn,iw
      common /abcdg/ iabsorp
      common /bdeo/ ivar
      common /a0ef1/ r0,z0,rm,cltn
      common /bcef/ ynz,ynpopq
      common /df/ pdec14,pdec24,pdec34,idec
      common /dg/ pintld4,pintcl4,pintal4
      parameter(zero=0.d0,hbeg=1.d-4,iturns=1,maxat=3,nvar=4) !sav2008
      dimension ystart(4),yscal(nvar),y(nvar),dydx(nvar),yold(nvar)
      ipr1=0
      iat=0
      x=zero
      eps1=eps
      hdrob1=hdrob
      rbord1=rbord
      hdid=zero
      pintld4=zero
      pintcl4=zero
      pintal4=zero
      pdec14zz=zero
      pdec24zz=zero
      pdec34zz=zero
      xold=x
      do i=1,nvar
        y(i)=ystart(i)
        yold(i)=y(i)
      end do
      rmm=1d+10*irs
!sav2008
!old      rexi1=rexi+rrange
!old      rexi2=rexi-rrange
!old      if(rexi1.gt.0.95d0) rexi1=1.d10
!old      if(rexi2.lt.0.05d0) rexi2=-1.d10
!est      if(rexi1.gt.0.9d0) rexi1=1.1d0
10    continue
c--------------------------------------
c start integration
c--------------------------------------
      do nstp=1,maxstep4
        idec=iturns
        call derivs(x,y,dydx)
        idec=0
        pintld4=pintld4+dabs((pdec14+pdec14zz)/2d0*hdid)
        pintcl4=pintcl4+dabs((pdec24+pdec24zz)/2d0*hdid)
        pintal4=pintal4+dabs((pdec34+pdec34zz)/2d0*hdid)
        pdec14zz=pdec14
        pdec24zz=pdec24
        pdec34zz=pdec34
       if(nstp.eq.1) then
        h=hbeg
!!var        if(dabs(dydx(3)).ne.zero) h=dabs(hmin1/dydx(3))/hdrob1
        if(dabs(dydx(3)).ne.zero) h=0.5d0*dabs(rrange/dydx(3))/hdrob1
       end if
20     continue
       if(y(3).ge.rbord1.and.dydx(3).gt.zero) then
c--------------------------------------
c forced reflection from periphery
c--------------------------------------
         ivar=3
         izn=-izn
         call disp2(y(3),y(2),y(1),xnr,prt,prm)
         if(ivar.eq.-1) then !out of dispersion curve - restart
          do i=1,nvar
           y(i)=ystart(i)
          end do
          x=zero
          iat=iat+1
          if(iat.gt.maxat) then
           if(ipri.gt.1) write (*,*)'turn in driver4 failed'
           goto 40
          end if
          eps1=eps1/2.d0
          hdrob1=hdrob1*2.d0
          ivar=0
          goto 10
         end if
        irs=-irs
        y(4)=xnr
        call derivs(x,y,dydx)
        if(dydx(3).gt.zero.and.ipri.gt.1) then
         write(*,*)'Unsuccesful turn: r, drds=',y(3),dydx(3)
        end if
        ivar=0
        iat=0
       end if
!sav2008       if((y(3).gt.rexi1.or.y(3).lt.rexi2)) then  ! exit
   !!    if(dabs(y(3)-rexi).gt.rrange.or.nstp.eq.maxstep4) then  ! exit !sav2008
       if(dabs(y(3)-rexi).gt.rrange) then  ! exit !sav2008
        if(dydx(3).gt.zero) irs=-1
        if(dydx(3).lt.zero) irs=1
        if(dydx(3).eq.zero) then !sav2008
         write(*,*)'exception dr/ds=0 in driver4'
         pause 'zmi na pedal'
         go to 1
        end if
        x2=x
        x1=rmm
        do i=1,nvar
         ystart(i)=y(i)
        end do
        return
       end if
1      continue
c---------------------------------------
c remember old values
c---------------------------------------
       xold=x
       do i=1,nvar
        dyd=dabs(dydx(i))
        yscal(i)=dabs(y(i))+dabs(h*dyd)+1.d-30/(1.d0+dyd)+1.d-30
        yold(i)=y(i)
       end do
       if (y(3)*irs.lt.rmm*irs) rmm=y(3)
30     continue
!!var       call rkqc(y,dydx,nvar,x,h,eps1,yscal,hdid,hnext,derivs)
       call rkqs(y,dydx,nvar,x,h,eps1,yscal,hdid,hnext,derivs)
       if(y(3).ge.1.d0) then  ! crossed plasma boundary
        do ii=1,nvar
         y(ii)=yold(ii)
        end do
        x=xold
        ipr1=ipr1+1
        if (ipr1.lt.maxat) then
         h=h/3.d0
         goto 30
        end if
        rbord1=y(3)-1.d-4
        goto 20
       end if
       ipr1=0
       if(dabs(hnext).lt.hmin) then
        if(ipri.gt.1) write(*,*)'error in dr4: step is too small'
        goto 40
       end if
       h=hnext
      end do
      if(ipri.gt.1) write(*,*)'error in dr4: too many steps.'
      if(ipri.gt.1) write(*,*)'tet=',y(1),'xm=',y(2),'xend=',y(3)
40    iabsorp=-1
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine disp2(pa,yn2,ptet,xnro,prt,prm)
      implicit real*8 (a-h,o-z)
      common /a0befr/ pi,pi2
      common /a0ef1/ r0,z0,rm,cltn
      common /a0ef2/ btor,ww
      common /a0ef3/ xmi,c0,c1,cnye,cnyi,xsz,vt0
      common /abcde/ izn,iw
      common /bcef/ ynz,ynpopq
      common /abefo/ yn3
      common /aef2/ icall1,icall2
      common /bdeo/ ivar
      common /be1/ xnr1,xnr2,xnr3,xnr4
      common /be2/ ider
      common /beo/ iroot
      common /cefn/ iconv,irefl
      common /ceg/ ipow,jfoundr
      common /eg1/ vfound,ifound
      common /eg2/ pdec1,pdec2,pdec3,pdecv,pdecal,dfdv,icf1,icf2
      common /eg3/ cf1,cf2,cf3,cf4,cf5,cf6
      common/a00/ xlog,zalfa,xmalfa,dn1,dn2,factor
      common /a0i5/ vperp(50,100),cnstal,zza,zze,valfa,kv
      common/b0/ itend0
      common /a0k/ cdl(10),cly(10),cgm(10),cmy(10),ncoef
      common /cnew/ inew !est !sav2008
      common/fj/dhdm,dhdnr,dhdtet,dhdr,ddn,dhdn3,dhdv2v,dhdu2u
      common/direct/znakstart
      common/metrika/g11,g12,g22,g33,gg,g,si,co
      common/fjham/ham
      parameter(zero=0.d0,one=1.d0,two=2.d0)
      parameter(clt=3.d10)
      iconv=0
      irefl=0
      if(pa.ge.one.or.pa.le.zero) goto70
      icall1=icall1+1
      xdl=fdf(pa,cdl,ncoef,xdlp)
      xly=fdf(pa,cly,ncoef,xlyp)
      xgm=fdf(pa,cgm,ncoef,xgmp)
      xmy=fdf(pa,cmy,ncoef,xmyp)
      xlyv=xlyp*pa+xly
      cotet=dcos(ptet)
      sitet=dsin(ptet)
      dxdr=-xdlp+cotet-xgmp*sitet**2
      dxdt=-(pa+two*xgm*cotet)*sitet
      dzdr=xlyv*sitet
      dzdt=xly*pa*cotet
      x0=r0/rm-xdl+pa*cotet-xgm*sitet**2
      dxdrdt=-sitet-two*xgmp*sitet*cotet
      dzdrdt=xlyv*cotet
      dxdtdt=-pa*cotet-two*xgm*(cotet**2-sitet**2)
      dzdtdt=-xly*pa*sitet
      x0t=dxdt
c--------------------------------------
c components of metric tensor
c--------------------------------------
      g11=dxdr**2+dzdr**2
      g22=dxdt**2+dzdt**2
      g12=dxdr*dxdt+dzdr*dzdt
      g33=x0**2
      xj=(dzdr*dxdt-dxdr*dzdt)**2  !gg=g11*g22-g12*g12
      gg=xj
      g=xj*g33
      g2v1=one/dsqrt(g22)
      g2jq=dsqrt(g22/xj)
      g3v=one/dsqrt(g33)
c--------------------------------------
c  magnetic field
c--------------------------------------
      bt=btor*(r0/rm)/x0
      bp=g2jq*g3v*xmy
      b=dsqrt(bp*bp+bt*bt)
      si=bp/b
      co=bt/b
      if(ivar.eq.1) return
c---------------------------------------
c components of dielectric tensor
c---------------------------------------
!sav2008      pn=fn(pa)
!var      pn=fn1(pa,fnr)
!!      pn=fn1(pa,fnr)
!!      pn=fn2(pa,fnr,fnrr) !sav2008
      if(inew.eq.0) then !vardens
       pn=fn1(pa,fnr)
      else
       pn=fn2(pa,fnr,fnrr)
      end if
      wpq=c0**2*pn
      whe=b*c1
      v=wpq/ww**2
      u1=whe/ww
      u=u1**2
      e1=one-v*(one/xmi-one/u)
      e2=v/u1
      e3=one-v
c-------------------------------------
c dispersion equation
c--------------------------------------
!sav2008      if(ivar.eq.2) yn2=(ynz-yn3*co*g3v)/(si*g2v1)
      ynz=yn2*si*g2v1+yn3*co*g3v
      ynzq=ynz**2
      as=e1
      bs=-(e1**2-e2**2+e1*e3-(e1+e3)*ynzq)
      cs=e3*(e1**2-e2**2-two*e1*ynzq+ynzq**2)
c-----------------------------
!est !sav2009
      pnew=zero
      yny= - (yn2*g2v1*co-yn3*g3v*si)
      if(inew.gt.0) then
       if(inew.eq.1) then
        yny= - (yn2*g2v1*co-yn3*g3v*si)
       else if(inew.eq.2) then
        yny= - g2jq*(yn2*g2v1*co-yn3*g3v*si)
       end if
       gpr=c0**2/ww**2/u1*fnr*xsz
       pnew=yny*gpr
       bs=bs+pnew
       cs=cs+pnew*(ynzq-e3)
      end if
c------------------------------------
      dls=bs*bs-4d0*as*cs

c      write(*,*)'rho=',pa,' teta=',ptet
c      write(*,*)'N2=',yn2,' N3=',yn3
c      write(*,*)'Npar=',ynz,' e1=',e1
c      write(*,*)'v=',v,' u=',u
c      write(*,*)'whe=',whe,' ww=',ww
c      write(*,*)'e2=',e2,' e3=',e3
c      write(*,*)'bs=',bs,' as=',as
c      write(*,*)'cs=',cs,' dls=',dls
c      pause
      
      if(dls.lt.zero) then
        goto (60,20,10) iroot
10        xnr1=1d+10
          xnr2=1d+10
          xnr3=1d+10
          xnr4=1d+10
          return
20      prt=dls
        prm=666d0
        return
      end if
30    continue
      dl1=dfloat(iw)*dsqrt(dls)/two/as
      if(iw.eq.-1) ynpopq=-bs/(two*as)+dl1
      if(iw.eq.1)  ynpopq=two*cs/(-bs-two*as*dl1)
      if(iroot.eq.3) ynpopq1=-bs/(two*as)-dl1

cc      write(*,*)'iw=',iw,' izn=',izn,' Nperp=',dsqrt(ynpopq)
cc      write(*,*)'Nperp2=',ynpopq,' ynpopq1=',-bs/(two*as)-dl1
cc      pause

      if(ynpopq.lt.zero.and.iroot.eq.1) goto70
      al=g22/xj
      bl=-yn2*g12/xj
      cl=g11*yn2**2/xj+yn3**2/g33-ynzq-ynpopq
      if(iroot.eq.3) cl1=g11*yn2**2/xj+yn3**2/g33-ynzq-ynpopq1
      dll=bl*bl-al*cl
      if(iroot.eq.2) then
       prt=dls
       prm=dll
       if(dll.ge.zero) then !sav2008
!!old variant:
cc        dl2=-dfloat(izn)*dsqrt(dll)/al
cc        if(izn.eq.1) xnr=-bl/al+dl2
cc        if(izn.eq.-1) xnr=cl/(-bl-al*dl2)
cc        xnro=xnr
cc       end if
cc       return
cc      end if
!!!!!!!!!!!!!!

!!new variant:
        izn=1
        dl2=-dsqrt(dll)/al
        xnr=-bl/al+dl2
        call dhdomega(pa,ptet,xnr,yn2)
cc        write(*,*)'#1: izn=',izn,' dl2=',dl2,' xnr=',xnr
cc        write(*,*)'znak=',znakstart,' -znak*dhdnr=',-znakstart*dhdnr
        if(-znakstart*dhdnr.gt.zero) then
         izn=-1
         dl2=dsqrt(dll)/al
         xnr=cl/(-bl-al*dl2)
         call dhdomega(pa,ptet,xnr,yn2)
cc         write(*,*)'#2: izn=',izn,' dl2=',dl2,' xnr=',xnr
cc         write(*,*)'znak=',znakstart,' -znak*dhdnr=',-znakstart*dhdnr
         if(-znakstart*dhdnr.gt.zero) then
          write(*,*)'Exception: both modes go outward !!'
          stop
         end if
        end if
        xnro=xnr
cc        pause
       end if
       return
      end if
      if(dll.lt.zero) goto(70,70,50) iroot
40    dl2=-dfloat(izn)*dsqrt(dll)/al
      if(izn.eq.1) xnr=-bl/al+dl2
      if(izn.eq.-1) xnr=cl/(-bl-al*dl2)
      xnro=xnr
      if(ivar.gt.1) then
ccccccc  find Nr of reflected wave
        dnx=two*as*ynpopq+bs
        dhdnr=dnx*(two*g22*xnr-two*g12*yn2)/xj
        if(-znakstart*dhdnr.gt.zero) then
          izn=-izn
          goto40
        end if
        return
      end if
50    if(iroot.eq.3) then
c---------------------------
c  find all roots
c----------------------------
        if(dll.ge.zero) then
          xnr1=xnr
          xnr2=-bl/al-dl2
        else
          xnr1=1d+10
          xnr2=1d+10
        end if
        dll1=bl**2-al*cl1
        if(dll1.lt.zero) then
          xnr3=1d+10
          xnr4=1d+10
        else
          xnr3=-bl/al-izn*dsqrt(dll1)/al
          xnr4=-bl/al+izn*dsqrt(dll1)/al
        end if
      end if
      if(ider.eq.0) then
        prt=0d0
        prm=0d0
        return
      end if
c--------------------------------------
c   calculation of derivatives
c--------------------------------------
      g11t=two*(dxdr*dxdrdt+dzdr*dzdrdt)
      g22t=two*(dxdt*dxdtdt+dzdt*dzdtdt)
      g33t=two*x0*(-pa*sitet-two*xgm*sitet*cotet)
      g12t=dxdrdt*dxdt+dxdr*dxdtdt+dzdrdt*dzdt+dzdr*dzdtdt
      xjt=g11t*g22+g22t*g11-two*g12*g12t
      btt=-btor*(r0/rm)/x0**2*x0t
      g2jqt=(g22t/xj-g22/xj**2*xjt)/(g2jq*two)
      bpt=xmy*(g2jqt*g3v-.5d0*g2jq*g3v/g33*g33t)
      bat=one/b*(bp*bpt+bt*btt)
      sit=bpt/b-bp/b**2*bat
      cot=btt/b-bt/b**2*bat
      u1t=c1/ww*bat
      e1t=-v/u**2*two*u1*u1t
      e2t=-v/u*u1t
      ynzt=yn2*sit*g2v1-yn2*si*g2v1**3/two*g22t+
     *     yn3*cot*g3v-yn3*co*g3v**3/two*g33t
      p1=two*ynz*ynzt
      p2=e1t
      p3=(e2*e2t)/e1-e2**2/(two*e1**2)*e1t

      s1=-p2/(two*e1**2)*e3*(ynzq-e1)+
     *(e3+e1)/(two*e1)*(p1-p2)+p3
      s2=two*e3/e1*(ynzq-e1)*(p1-p2)-
     *p2*e3/(e1**2)*(ynzq-e1)**2-two*e3*p3
      dnm=two*ynz*si*g2v1
      v1=(e3+e1)/(two*e1)*dnm
      v2=two*e3/e1*(ynzq-e1)*dnm
c-----------------------------------
!est !sav2009
      if(inew.gt.0) then
       gprt=-c0**2/ww**2*fnr/u*u1t*xsz
       if(inew.eq.1) then
        ynyt= - (yn2*cot*g2v1-yn2*co*g2v1**3/two*g22t
     &            -(yn3*sit*g3v-yn3*si*g3v**3/two*g33t))
        dnym= - co*g2v1
       else if(inew.eq.2) then
        ynyt= - g2jq*(yn2*cot*g2v1-yn2*co*g2v1**3/two*g22t
     &                 -(yn3*sit*g3v-yn3*si*g3v**3/two*g33t))
        ynyt=ynyt - g2jqt*(yn2*g2v1*co-yn3*g3v*si)
        dnym= - g2jq*co*g2v1
       end if
       pnewt=(ynyt*gpr+yny*gprt)
       s1=s1+pnewt/(two*e1)-pnew/(two*e1**2)*e1t
       s2=s2+(pnewt*(ynzq-e3)+pnew*p1)/e1-pnew*(ynzq-e3)/e1**2*e1t
       v1=v1+dnym*gpr/(two*e1)
       v2=v2+gpr*(dnym*(ynzq-e3)+yny*dnm)/e1
      end if
c---------------------------------------------
      vvt=-s1+(bs/as*s1-s2)/(two*dl1)
      vvm=-v1+(bs/as*v1-v2)/(two*dl1)
      s1=-yn2*(g12t/g22-g12/g22**2*g22t)
      s21=yn2**2*(g11t/g22-g11/g22**2*g22t)
      s22=yn3**2*( xjt/(g33*g22)-
     *xj/(g33*g22)**2*(g33t*g22+g22t*g33) )
      sjg=(xjt*g22-xj*g22t)/g22**2
      s23=two*ynz*ynzt*xj/g22+sjg*ynzq
      s24=vvt*xj/g22+ynpopq*sjg
      s2=s21+s22-s23-s24
      prt=-s1+(two*(bl/al)*s1-s2)/(two*dl2)
      s1=-g12/g22
      s21=two*yn2*g11/g22
      s22=dnm*xj/g22
      s23=vvm*xj/g22
      s2=s21-s22-s23
      prm=-s1+(two*(bl/al)*s1-s2)/(two*dl2)
      if(ipow.gt.0) then
c--------------------------------------
c  calculation of decrements
c--------------------------------------
        dnx=two*as*ynpopq+bs
        dhdnr=dnx*(two*g22*xnr-two*g12*yn2)/xj
        sl1=(ynzq-e1)*(ynzq+ynpopq-e1)-e2**2
        cf3=ynz
        cf4=xnr
        cf5=yn2
        vz=cltn/dabs(ynz)
        if(vz.gt.cltn) vz=cltn !sav2010
        vt=fvt(pa)
        jr=jfoundr
        icf1=iw
        icf2=izn
        call distr(vz,jr,ifound,fder)
        dfdv=fder
        vfound=vz
        cf2=ptet
	cf6=yny
        aimh=wpq/ww**2*pi*sl1*cltn**2/ynzq
        pdecv=dabs(aimh/dhdnr/xsz)
!!        pdec1=-pdecv*dfdv
        pdec1=dabs(pdecv*dfdv)
        pnye=cnye*wpq**2/(pn*vt**3)
        pnyi=cnyi*pnye*zefff(pa)
        pdec2=dabs(pnyi/ww*(wpq/whe**2*ynpopq+wpq/ww**2*ynzq)*ynpopq/
     *  dhdnr/xsz)
        cf1=dsqrt(ynpopq)
        if(itend0.gt.0) then
          tmp=ft(pa)/0.16d-8
          fcoll=.5d-13*pn*zalfa**2*xlog/xmalfa/tmp**1.5d0
cc          ddens=dn1*pn
cc          tdens=dn2*pn
cc          tt=fti(pa)**0.33333d0    ! (ti, kev)^1/3
cc          source=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
           call source_new(pa,source)
           dek1=cnstal*pdecv*(1.d0-e3/ynpopq)**2/cf1
           dek2=source/(fcoll*pn)
           pdecal=dek1*dek2
           pdec3=zero
             if(itend0.gt.0) then
              argum=clt/(cf1*valfa)
              dek3=zatukh(argum,jr,vperp,kv)
              pdec3=pdecal*dek3
             end if
        end if
      end if
      return
c   conversion
60    iconv=1
      if (ivar.ne.0) ivar=-1
      return
c    reflection
70    irefl=1
      if (ivar.gt.1.and.ivar.ne.10) then
        iw=-iw
        ivar=10
        goto30
      end if
      if (ivar.eq.10) ivar=-1
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine disp4(pa,ptet,xnr,yn2)
      implicit real*8 (a-h,o-z)
      common /a0befr/ pi,pi2
      common /a0ef1/ r0,z0,rm,cltn
      common /a0ef2/ btor,ww
      common /a0ef3/ xmi,c0,c1,cnye,cnyi,xsz,vt0
      common /bcef/ ynz,ynpopq
      common /abefo/ yn3
      common /aef2/ icall1,icall2
      common /cefn/ iconv,irefl
      common /df/ pdec14,pdec24,pdec34,idec
      common/a00/ xlog,zalfa,xmalfa,dn1,dn2,factor
      common /a0i5/ vperp(50,100),cnstal,zza,zze,valfa,kv
      common/b0/ itend0
      common /a0k/ cdl(10),cly(10),cgm(10),cmy(10),ncoef
      common /cnew/ inew !est !sav2008
      common/plasma/v,u,e1,e2,e3,dvdr,dudr,dudt
      common/metrika/g11,g12,g22,g33,gg,g,si,co
      common/fj/dhdm,dhdnr,dhdtet,dhdr,ddn,dhdn3,dhdv2v,dhdu2u
      common/fjham/ham
      parameter(zero=0d0, one=1d0, two=2d0)
      irefl=0
      iconv=0
      if(pa.eq.zero) pa=1.d-7
      if(pa.lt.zero) pa=dabs(pa)
!sav2008      if(pa.gt.one) then
!sav2008       dhdm=666d0
!sav2008       dhdtet=-666d0
!sav2008       dhdnr=666d0
!sav2008       dhdr=-666d0
!sav2008       irefl=1
!sav2008       return
!sav2008      end if

      icall2=icall2+1
!!      pn=fn1(pa,fnr)
!!      pn=fn2(pa,fnr,fnrr)
      if(inew.eq.0) then !vardens
       pn=fn1(pa,fnr)
      else
       pn=fn2(pa,fnr,fnrr)
      end if

cc        hstp=1.d-7
cc        pplus=fn2(pa+hstp,fnr2,fnrr2)
cc        pminus=fn2(pa-hstp,fnr1,fnrr1)
cc        fnr=0.5d0*(pplus-pminus)/hstp
cc        fnrr=0.5d0*(fnr2-fnr1)/hstp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      wpq=c0**2*pn
      xdl=fdfddf(pa,cdl,ncoef,xdlp,xdlpp)
      xly=fdfddf(pa,cly,ncoef,xlyp,xlypp)
      xgm=fdfddf(pa,cgm,ncoef,xgmp,xgmpp)
      xmy=fdf(pa,cmy,ncoef,xmyp)
      cotet=dcos(ptet)
      sitet=dsin(ptet)
      xlyv=xly+xlyp*pa
c--------------------------------------
c components of metric tensor
c--------------------------------------
      dxdr=-xdlp+cotet-xgmp*sitet**2
      dxdt=-(pa+two*xgm*cotet)*sitet
      dzdr=xlyv*sitet
      dzdt=xly*pa*cotet
      x0=r0/rm-xdl+pa*cotet-xgm*sitet**2
      dxdrdr=-xdlpp-xgmpp*sitet**2
      dxdtdt=-cotet*(pa+two*xgm*cotet)+sitet**2*two*xgm
      dxdtdr=-sitet*(one+two*xgmp*cotet)
      dxdrdt=dxdtdr
      dzdrdr=(two*xlyp+pa*xlypp)*sitet
      dzdtdt=-xly*pa*sitet
      dzdtdr=xlyv*cotet
      dzdrdt=dzdtdr
      x0t=dxdt
      x0r=dxdr
      g11=dxdr**2+dzdr**2
      g22=dxdt**2+dzdt**2
      g12=dxdr*dxdt+dzdr*dzdt
      g33=x0**2
      xj=(dzdr*dxdt-dxdr*dzdt)**2  !gg=g11*g22-g12*g12
      gg=xj
      g=xj*g33
      g2jq=dsqrt(g22/xj)
      g2gq=dsqrt(g22/g)
      g22q=dsqrt(g22)
      g33q=dsqrt(g33)
c--------------------------------------
c  magnetic field
c--------------------------------------
      bt=btor*(r0/rm)/x0
      bp=g2gq*xmy
      b=dsqrt(bp*bp+bt*bt)
      whe=b*c1
      si=bp/b
      co=bt/b
c---------------------------------------
c components of dielectric tensor
c---------------------------------------
      v=wpq/ww**2
      u1=whe/ww
      u=u1**2
      e1=one-v*(one/xmi-one/u)
      e2=v/u1
      e3=one-v
c-------------------------------------
c dispersion equation
c--------------------------------------
      ynz=yn2*si/g22q+yn3*co/g33q
      ynzq=ynz**2
      vpop=xnr**2*g22-two*xnr*yn2*g12+g11*yn2**2
      ynpopq=vpop/xj+yn3**2/g33-ynzq
      as=e1
      bs=-(e1**2-e2**2+e1*e3-(e1+e3)*ynzq)
      cs=e3*(e1**2-e2**2-two*e1*ynzq+ynzq**2)
c----------------------------------------------------
!sav2009
      dhdv=(1.d0/(u**2*xmi**2))*((2.d0-2.d0*ynpopq-3.d0*v)*v*xmi**2-
     &     u**2*(3.d0*v**2+2.d0*v*(-1.d0+ynpopq*(1.d0+xmi)
     &     +2.d0*xmi*(-1.d0+ynzq))+
     &     xmi*(-2.d0+ynpopq+xmi*(-1.d0+ynzq))*(-1.d0+ynpopq+ynzq))
     &     +u*xmi*(3.d0*v**2*(2.d0+xmi)
     &     +(-2.d0+ynpopq)*xmi*(-1.d0+ynpopq+ynzq)
     &     +v*(-4.d0+4.d0*ynpopq*(1.d0+xmi)+xmi*(-6.d0+4.d0*ynzq))))
      dhdu=-(1.d0/(u**3*xmi))*(v*(-1.d0+ynpopq+v)*(2.d0*u*v-2.d0*v*xmi 
     &     +u*(-2.d0+ynpopq+v)*xmi)+u*v*(-2.d0+ynpopq+2.d0*v)*xmi*ynzq)
      dhdv2v=2.d0*v*dhdv !w*d(-H)/dv
      dhdu2u=2.d0*u*dhdu !w*d(-H)/du
c----------------------------------------------------
!est !sav2009
      if(inew.gt.0) then
       if(inew.eq.1) then
        yny= - (yn2*co/g22q-yn3*si/g33q)
       else if(inew.eq.2) then
        yny= - g2jq*(yn2*co/g22q-yn3*si/g33q)
       end if
       gpr=c0**2/ww**2/u1*fnr*xsz
       gdop=yny*gpr
       bs=bs+gdop
       cs=cs+gdop*(ynzq-e3)
!sav2009:
       dgpr=-gpr        !w*d(gpr)/dw
       wde3dw=2.d0*v    !w*d(e3)/dw
       wdbsdw=yny*dgpr  !w*d(bs)/dw
       wdcsdw=(ynzq-e3)*wdbsdw-gdop*wde3dw   !w*d(cs)/dw
       wdhdw=wdbsdw*ynpopq+wdcsdw   !w*d(H1)/dw
       dhdv2v=dhdv2v-wdhdw !correction to dhdv2v: w*d(-H)/dv+w*d(-H1)/dw
      end if
      ham=as*ynpopq**2+bs*ynpopq+cs !sav2009
c--------------------------------------------------------
!!      dl=bs**2-4d0*as*bs
c--------------------------------------
c   calculation of derivatives
c--------------------------------------
      g11r=two*dxdr*dxdrdr+two*dzdr*dzdrdr
      g22r=two*dxdt*dxdtdr+two*dzdt*dzdtdr
      g11t=two*dxdr*dxdrdt+two*dzdr*dzdrdt
      g22t=two*dxdt*dxdtdt+two*dzdt*dzdtdt
      g12r=dxdrdr*dxdt+dxdr*dxdtdr+dzdrdr*dzdt+dzdr*dzdtdr
      g12t=dxdrdt*dxdt+dxdr*dxdtdt+dzdrdt*dzdt+dzdr*dzdtdt
      g33r=two*x0*x0r
      g33t=two*x0*x0t
      g22qr=g22r/(g22q*two)
      g22qt=g22t/(g22q*two)
      g33qr=g33r/(g33q*two)
      g33qt=g33t/(g33q*two)
      xjr=g11r*g22+g22r*g11-two*g12*g12r
      xjt=g11t*g22+g22t*g11-two*g12*g12t
      g2jqr=(g22r/xj-g22/xj**2*xjr)/(g2jq*two) !sav2009
      g2jqt=(g22t/xj-g22/xj**2*xjt)/(g2jq*two) !sav2009
      gr=xjr*g33+g33r*xj
      gt=xjt*g33+g33t*xj
      g2gqt=(g22t/g-g22/g**2*gt)/(g2gq*two)
      g2gqr=(g22r/g-g22/g**2*gr)/(g2gq*two)
      bpt=xmy*g2gqt
      bpr=g2gqr*xmy+g2gq*xmyp
      btr=-btor*(r0/rm)/x0**2*x0r
      btt=-btor*(r0/rm)/x0**2*x0t
      bat=one/b*(bp*bpt+bt*btt)
      bar=one/b*(bp*bpr+bt*btr)
      sit=bpt/b-bp/b**2*bat
      cot=btt/b-bt/b**2*bat
      sir=bpr/b-bp/b**2*bar
      cor=btr/b-bt/b**2*bar
      dvdr=fnr*c0**2/ww**2
      du1dr=c1*bar/ww
      dudr=two*u1*du1dr
      du1dt=c1*bat/ww
      dudt=two*u1*du1dt
      e1r=-dvdr*(one/xmi-one/u)-v*dudr/u**2
      e1t=-v*dudt/u**2
      e2r=dvdr/u1-v/u1**2*du1dr
      e2t=-v/u1**2*du1dt
      e3r=-dvdr
      ynzr=yn2*(sir/g22q-si/g22q**2*g22qr)+
     *     yn3*(cor/g33q-co/g33q**2*g33qr)
      ynzt=yn2*(sit/g22q-si/g22q**2*g22qt)+
     *     yn3*(cot/g33q-co/g33q**2*g33qt)
      ynzqr=two*ynz*ynzr
      ynzqt=two*ynz*ynzt
      vpopr=(xnr**2*g22r-two*xnr*yn2*g12r+yn2**2*g11r)
      vpopt=(xnr**2*g22t-two*xnr*yn2*g12t+yn2**2*g11t)
      ynpopqr=vpopr/xj-vpop/xj**2*xjr-yn3**2/g33**2*g33r-ynzqr
      ynpopqt=vpopt/xj-vpop/xj**2*xjt-yn3**2/g33**2*g33t-ynzqt
      asr=e1r
      bsr=(e3r+e1r)*(ynzq-e1)+(e3+e1)*(ynzqr-e1r)+two*e2*e2r
      csr=e3r*((ynzq-e1)**2-e2**2)+e3*(two*(ynzq-e1)*(ynzqr-e1r)-
     *    two*e2*e2r)
      ast=e1t
      bst=e1t*(ynzq-e1)+(e3+e1)*(ynzqt-e1t)+two*e2*e2t
      cst=e3*(two*(ynzq-e1)*(ynzqt-e1t)-two*e2*e2t)
c---------------------------------------------------
!est !sav2009
      if (inew.gt.0) then
       if(inew.eq.1) then
        ynyr= - (yn2*(cor/g22q-co/g22q**2*g22qr)
     *            -yn3*(sir/g33q-si/g33q**2*g33qr))
        ynyt= - (yn2*(cot/g22q-co/g22q**2*g22qt)-
     *            -yn3*(sit/g33q-si/g33q**2*g33qt))
       else if(inew.eq.2) then
        ynyr= - g2jq*(yn2*(cor/g22q-co/g22q**2*g22qr)
     *                 -yn3*(sir/g33q-si/g33q**2*g33qr))
        ynyr=ynyr - g2jqr*(yn2*co/g22q-yn3*si/g33q)
        ynyt= - g2jq*(yn2*(cot/g22q-co/g22q**2*g22qt)-
     *                 -yn3*(sit/g33q-si/g33q**2*g33qt))
        ynyt=ynyt - g2jqt*(yn2*co/g22q-yn3*si/g33q)
       end if
       gprr=c0**2/ww**2*(fnrr/u1-fnr/u1**2*du1dr)*xsz
       gprt=-c0**2/ww**2*fnr/u1**2*du1dt*xsz
       gdopr=ynyr*gpr+yny*gprr
       gdopt=ynyt*gpr+yny*gprt
       bsr=bsr+gdopr
       csr=csr+gdopr*(ynzq-e3)+gdop*(ynzqr-e3r)
       bst=bst+gdopt
       cst=cst+gdopt*(ynzq-e3)+gdop*ynzqt
      end if
c---------------------------------------------------

      dhdr=asr*ynpopq**2+bsr*ynpopq+
     *     as*two*ynpopq*ynpopqr+bs*ynpopqr+csr
      dhdtet=ast*ynpopq**2+bst*ynpopq+
     *       as*two*ynpopq*ynpopqt+bs*ynpopqt+cst
      dnx=two*as*ynpopq+bs
      dnz=ynpopq*(e1+e3)+two*(ynzq-e1)*e3
      dhdnr=dnx*two*(g22*xnr-g12*yn2)/xj
      dhdm=dnx*two*(yn2*g11-xnr*g12)/xj+(dnz-dnx)*two*ynz*si/g22q
!sav2009
      dhdn3=two*((yn3-ynz*co*g33q)*dnx+ynz*co*g33q*dnz)/g33 !sav2009
c----------------------------------------------------------------
!est !sav2009
      if (inew.gt.0) then
       if(inew.eq.1) then
        dny= - gpr*(ynpopq+ynzq-e3)
       else if(inew.eq.2) then
        dny= - g2jq*gpr*(ynpopq+ynzq-e3)
       end if
       dhdm=dhdm+dny*co/g22q+two*ynz*yny*gpr*si/g22q
       dhdn3=dhdn3-dny*si/g33q+two*ynz*yny*gpr*co/g33q !sav2009
      end if
c---------------------------------------------------
!sav2009
      ddn2=g11*dhdnr**2+g22*dhdm**2+2.d0*g12*dhdnr*dhdm+g33*dhdn3**2
      ddn=dsqrt(ddn2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(idec.ne.0) then
       vt=fvt(pa)
       sl1=(ynzq-e1)*(ynzq+ynpopq-e1)-e2**2
       aimh=wpq/ww**2*pi*sl1*cltn**2/ynzq
       pdec14=dabs(aimh/xsz/ddn)
       pnye=cnye*wpq**2/(pn*vt**3)
       pnyi=cnyi*pnye*zefff(pa)
       pdec24=dabs(pnyi/ww*(wpq/whe**2*ynpopq+wpq/ww**2*ynzq)*ynpopq/
     * xsz/ddn)
       if(itend0.gt.0) then
        tmp=ft(pa)/0.16d-8
        fcoll=.5d-13*pn*zalfa**2*xlog/xmalfa/tmp**1.5d0
cc        ddens=dn1*pn
cc        tdens=dn2*pn
cc        tt=fti(pa)**0.33333d0    ! (ti, kev)^1/3
cc        source=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
        call source_new(pa,source)
        dek1=cnstal*pdec14*(1.d0-e3/ynpopq)**2/dsqrt(ynpopq)
        dek2=source/(fcoll*pn)
        pdec34=dek1*dek2
       end if
      end if
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine extd2(x,y,dydx)
      implicit real*8 (a-h,o-z)
      dimension y(*),dydx(*)
      tt=y(1)
      xm=y(2)
      call disp2(x,xm,tt,xnr,prt,prm)
      dydx(1)=-prm
      dydx(2)=prt
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine extd4(x,y,dydx)
      implicit real*8 (a-h,o-z)
      dimension y(*),dydx(*)
      common/fj/dhdm,dhdnr,dhdtet,dhdr,ddn,dhdn3,dhdv2v,dhdu2u
      common/direct/znakstart
      znak=znakstart
      xxx=x
      ptet=y(1)
      yn2=y(2)
      pa=y(3)
      yn1=y(4)
      call disp4(pa,ptet,yn1,yn2)

!new variant
      dydx(1)=-znak*dhdm/ddn
      dydx(2)=znak*dhdtet/ddn
      dydx(3)=-znak*dhdnr/ddn
      dydx(4)=znak*dhdr/ddn
      dydx(5)=-znak*dhdn3/ddn

c      dydx(1)=znak*dhdm/ddn
c      dydx(2)=-znak*dhdtet/ddn
c      dydx(3)=znak*dhdnr/ddn
c      dydx(4)=-znak*dhdr/ddn
c      dydx(5)=znak*dhdn3/ddn

!old variant:
!      dydx(1)=dhdm/ddn
!      dydx(2)=-dhdtet/ddn
!      dydx(3)=dhdnr/ddn
!      dydx(4)=-dhdr/ddn
!      dydx(5)=dhdn3/ddn
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dhdomega(rho,theta,yn1,yn2)
      implicit real*8 (a-h,o-z)
      common /a0ef2/ btor,ww
      common /abefo/ yn3
      common/fj/dhdm,dhdnr,dhdtet,dhdr,ddn,dhdn3,dhdv2v,dhdu2u
      common/direct/znakstart
      parameter(zero=0.d0,h=1.d-6)

      call disp4(rho,theta,yn1,yn2)

!!w*dH/dw=wdhdw:
      wdhdw=-(yn1*dhdnr+yn2*dhdm+yn3*dhdn3+dhdv2v+dhdu2u)
      znak=dsign(1.d0,wdhdw)
      znakstart=znak
c      write(*,*)'formula: znak=',znak
c      write(*,*)'wdhdw=',wdhdw,' H=',ham
c      write(*,*)'rho=',rho,' teta=',theta
c      write(*,*)'yn1=',yn1,' yn2=',yn2
c      write(*,*)'dhdnr=',dhdnr,' dhdm=',dhdm
c      write(*,*)'dhdr=',dhdr,' dhdtet=',dhdtet
c      write(*,*)'dhdn3=',dhdn3,' yn3=',yn3
c      write(*,*)'yn1*dhdnr=',yn1*dhdnr,' yn2*dhdm=',yn2*dhdm
c      write(*,*)'yn1*dhdnr+yn2*dhdm=',yn1*dhdnr+yn2*dhdm
cc      pause

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine difeq(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
      implicit none
      external derivs
      integer nv,nmax,kmaxx,imax
      double precision eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv)
     *,safe1,safe2,redmax,redmin,tiny,scalmx
     &,dysav(nv) !sav#
      parameter(nmax=50,kmaxx=8,imax=kmaxx+1,safe1=.25d0,safe2=.7d0
     *,redmax=1.d-5,redmin=.7d0,tiny=1.d-30,scalmx=.1d0)
cu    uses derivs,mmid,pzextr
      integer i,iq,k,kk,km,kmax,kopt,nseq(imax)
      double precision eps1,epsold,errmax,fact,h,red,scale,work,wrkmin
     *,xest,xnew,a(imax),alf(kmaxx,kmaxx),err(kmaxx)
     *,yerr(nmax),ysav(nmax),yseq(nmax)
      logical first,reduct
      save a,alf,epsold,first,kmax,kopt,nseq,xnew
      double precision hmin1,dyd
      integer ii,ind
      common /a0cdm/ hmin1
      common /cmn/ ind
      data first/.true./,epsold/-1.d0/
      data nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=safe1*eps
        a(1)=nseq(1)+1
        do 11 k=1,kmaxx
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,kmaxx
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.d0)*(2*k+
     *1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,kmaxx-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
        dysav(i)=dydx(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if(xnew.eq.x) write(*,*) 'step size underflow in difeq'
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
 !sav#
      if(ind.eq.1) then
        h=h/2d0
        if (dabs(h).lt.hmin1) then
          do ii=1,nv
            y(ii)=ysav(ii)
          end do
          hnext=h
          return
        end if
        do ii=1,nv
          dyd=dabs(dysav(ii))
          yscal(ii)=dabs(ysav(ii))+dabs(h*dyd)+1.d-30/(1d0+dyd)
          y(ii)=ysav(ii)
        end do
        goto 2
      end if
!sav#
        xest=(h/nseq(k))**2
!var        call pzextr(k,xest,yseq,y,yerr,nv)  !polynomial extrapolation
        call rzextr(k,xest,yseq,y,yerr,nv) !rational extrapolation
        if(k.ne.1)then
          errmax=tiny
          do 16 i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/safe1)**(1.d0/(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.d0)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=safe2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1.d0/err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*safe2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,redmin)
      red=max(red,redmax)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.d35
      do 18 kk=1,km
        fact=max(err(kk),scalmx)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),scalmx)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
      implicit none
      external derivs
      integer nstep,nvar,nmax
      double precision htot,xs,dydx(nvar),y(nvar),yout(nvar)
      parameter (nmax=50)
      integer i,n
      double precision h,h2,swap,x,ym(nmax),yn(nmax)
      double precision yz1,yz2
      integer iconv,irefl,ind
      common /cefn/ iconv,irefl
      common /cmn/ ind
      h=htot/nstep
      yz1=y(1) !sav#
      yz2=y(2) !sav#
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout)
      if (iconv+irefl.ne.0) goto 10 !sav#
      h2=2.d0*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout)
        if (iconv+irefl.ne.0) goto 10 !sav#
13    continue
      do 14 i=1,nvar
        yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14    continue
      ind=0 !sav#
      return
10    ind=1 !sav#
      yout(1)=yz1 !sav#
      yout(2)=yz2 !sav#
      return !sav#
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rzextr(iest,xest,yest,yz,dy,nv)
      integer iest,nv,imax,nmax
      double precision xest,dy(nv),yest(nv),yz(nv)
      parameter (imax=13,nmax=50)
      integer j,k
      double precision b,b1,c,ddy,v,yy,d(nmax,imax),fx(imax),x(imax)
      save d,x
      x(iest)=xest
      if(iest.eq.1) then
        do 11 j=1,nv
          yz(j)=yest(j)
          d(j,1)=yest(j)
          dy(j)=yest(j)
11      continue
      else
        do 12 k=1,iest-1
          fx(k+1)=x(iest-k)/xest
12      continue
        do 14 j=1,nv
          yy=yest(j)
          v=d(j,1)
          c=yy
          d(j,1)=yy
          do 13 k=2,iest
            b1=fx(k)*v
            b=b1-c
            if(b.ne.0.d0) then
              b=(c-v)/b
              ddy=c*b
              c=b1*b
            else
              ddy=v
            endif
            if (k.ne.iest) v=d(j,k)
            d(j,k)=ddy
            yy=yy+ddy
13        continue
          dy(j)=ddy
          yz(j)=yy
14      continue
      endif
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pzextr(iest,xest,yest,yz,dy,nv)
      integer iest,nv,imax,nmax
      double precision xest,dy(nv),yest(nv),yz(nv)
      parameter (imax=13,nmax=50)
      integer j,k1
      double precision delta,f1,f2,q,d(nmax),qcol(nmax,imax),x(imax)
      save qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1.d0/(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8 function rini(xm,tet,xnr,yn,hr,ifail) !sav2009
      implicit real*8 (a-h,o-z)
      dimension vgrp(3),vph(3)
      common /bcef/ ynz,ynpopq
      common /abefo/ yn3
      common /beo/ iroot
      common /bdeo/ ivar
      common /cnew/ inew !est !sav2008
      common/metrika/g11,g12,g22,g33,gg,g,si,co !sav2009
      parameter(zero=0.d0,rhostart=1.d0,ntry_max=5)

      ifail=1
      rini=zero
      ntry=0
      pa=rhostart
      do while(ntry.lt.ntry_max.and.pa.ge.2d0*hr)
       pa=rhostart-hr*dble(ntry)-1.d-4
       ntry=ntry+1
       ivar=1
       call disp2(pa,xm,tet,xnr,prt,prm)
       if(inew.gt.0) then !g' in ST and poloidal grill direction
        yn3=zero                 !Nfi=0
        xm=yn*dsqrt(g22)/si      !given Npar at Nfi=0
!!        xm=yn*dsqrt(g22)         !given yn=(N*jpol) at Nfi=0
       else !usual tokamak and toroidal grill direction
        xm=zero               !N2=0
        yn3=yn*dsqrt(g33)/co  !if given Npar at Nteta=0
!!        yn3=yn*dsqrt(g33)       !if given Nfi at Nteta=0
       end if
       ivar=0
       iroot=2
       call disp2(pa,xm,tet,xnr,f1,f2)
       if(f1.ge.zero.and.f2.ge.zero) then
        rini=pa
        ifail=0
        return
       end if
      end do
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine distr(vz,j,ifound,fder)
      implicit none
      real*8 vz,fder
      integer j,ifound,i,klo,khi,ierr,nvpt,nvp
      real*8,dimension(:),allocatable:: vzj,dfdvj
      real*8 vgrid,dfundv,vlf,vrt,dflf,dfrt,dfout
      common/gridv/vgrid(101,100),dfundv(101,100),nvpt
      common /a0ghp/ vlf,vrt,dflf,dfrt
      nvp=nvpt
      allocate(vzj(nvp),dfdvj(nvp))
      do i=1,nvp
       vzj(i)=vgrid(i,j)
       dfdvj(i)=dfundv(i,j)
      end do
      call lock2(vzj,nvp,vz,klo,khi,ierr)
      if(ierr.eq.0) then !vgrid(1,j) <= vz <= vgrid(nvpt,j)
       call linf(vzj,dfdvj,vz,dfout,klo,khi)
       ifound=klo
       vlf=vzj(klo)
       vrt=vzj(khi)
       fder=dfout
       dflf=dfdvj(klo)
       dfrt=dfdvj(khi)
      else if(ierr.eq.1) then !vz < vgrid(1,j)
       write(*,*)'exception: ierr=1 in distr()'
       pause'next key = stop'
       stop
      else if(ierr.eq.2) then !vz > vgrid(nvpt,j)
       write(*,*)'exception: ierr=2 in distr()'
       pause'next key = stop'
       stop
      else if(ierr.eq.3) then
       write(*,*)'exception in distr, klo=khi=',klo,' j=',j,' nvp=',nvp
       write(*,*)'vz=',vz,' v1=',vzj(1),' v2=',vzj(nvp)
       pause'next key = stop'
       stop
      end if
      deallocate(vzj,dfdvj)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      double precision  function obeom(ptet,pa)
      implicit real*8 (a-h,o-z)
      common /a0befr/ pi,pi2
      common /a0ef1/ r0,z0,rm,cltn
      common /a0k/ cdl(10),cly(10),cgm(10),cmy(10),ncoef
      parameter(two=2.d0, pa0=0.d0)
      xdl=fdf(pa,cdl,ncoef,xdlp)
      xly=fdf(pa,cly,ncoef,xlyp)
      xgm=fdf(pa,cgm,ncoef,xgmp)
      xlyv=xlyp*pa+xly
      cotet=dcos(ptet)
      sitet=dsin(ptet)
      dxdr=-xdlp+cotet-xgmp*sitet**2
      dxdt=-(pa+two*xgm*cotet)*sitet
      dzdr=xlyv*sitet
      dzdt=xly*pa*cotet
      x0=r0/rm-xdl+pa*cotet-xgm*sitet**2
      dxdrdt=-sitet-two*xgmp*sitet*cotet
      dzdrdt=xlyv*cotet
      dxdtdt=-pa*cotet-two*xgm*(cotet**2-sitet**2)
      dzdtdt=-xly*pa*sitet
      x0t=dxdt
c--------------------------------------
c components of metric tensor
c--------------------------------------
      g11=dxdr**2+dzdr**2
      g22=dxdt**2+dzdt**2
      g12=dxdr*dxdt+dzdr*dzdt
      g33=x0**2
      xj=(dzdr*dxdt-dxdr*dzdt)**2  !gg=g11*g22-g12*g12
      g=xj*g33
      obeom=dsqrt(g)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision  function ploshad(ptet,pa)
      implicit real*8 (a-h,o-z)
      common /a0befr/ pi,pi2
      common /a0ef1/ r0,z0,rm,cltn
      common /a0k/ cdl(10),cly(10),cgm(10),cmy(10),ncoef
      parameter(two=2.d0, pa0=0.d0)
      xdl=fdf(pa,cdl,ncoef,xdlp)
      xly=fdf(pa,cly,ncoef,xlyp)
      xgm=fdf(pa,cgm,ncoef,xgmp)
      xlyv=xlyp*pa+xly
      cotet=dcos(ptet)
      sitet=dsin(ptet)
      dxdr=-xdlp+cotet-xgmp*sitet**2
      dxdt=-(pa+two*xgm*cotet)*sitet
      dzdr=xlyv*sitet
      dzdt=xly*pa*cotet
      x0=r0/rm-xdl+pa*cotet-xgm*sitet**2
      dxdrdt=-sitet-two*xgmp*sitet*cotet
      dzdrdt=xlyv*cotet
      dxdtdt=-pa*cotet-two*xgm*(cotet**2-sitet**2)
      dzdtdt=-xly*pa*sitet
      x0t=dxdt
c--------------------------------------
c components of metric tensor
c--------------------------------------
      g11=dxdr**2+dzdr**2
      g22=dxdt**2+dzdt**2
      g12=dxdr*dxdt+dzdr*dzdt
      xj=(dzdr*dxdt-dxdr*dzdt)**2  !gg=g11*g22-g12*g12
      ploshad=dsqrt(xj)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision  function fvt(r)
      implicit real*8 (a-h,o-z)
      pt=ft(r)
      fvt=dsqrt(pt/9.11d-28)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision  function ft(x)
! electron temperature, erg
      implicit real*8 (a-h,o-z)
      common /a0l3/ rh(501),y2dn(501),y2tm(501),y2tmi(501)
      common /a0l4/ con(501),tem(501),temi(501),nspl
      parameter(zero=0.d0,alfa=4.d0,dr=.02d0)
      pa=dabs(x) !#@sav
      if(pa.le.rh(nspl)) then
       call splnt(rh,tem,y2tm,nspl,pa,y,dy)
      else
       r=pa-rh(nspl)
       y=tem(nspl)*dexp(-alfa*(r/dr)**2)
      end if
!!      ft=y            ! kev
      ft=y*0.16d-8      ! erg
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision  function fti(x)
! ion temperature, kev
      implicit real*8 (a-h,o-z)
      common /a0l3/ rh(501),y2dn(501),y2tm(501),y2tmi(501)
      common /a0l4/ con(501),tem(501),temi(501),nspl
      parameter(zero=0.d0,alfa=4.d0,dr=.02d0)
      pa=dabs(x) !#@sav
      if(pa.le.rh(nspl)) then
       call splnt(rh,temi,y2tmi,nspl,pa,y,dy)
      else
       r=pa-rh(nspl)
       y=temi(nspl)*dexp(-alfa*(r/dr)**2)
      end if
      fti=y              ! kev
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision  function zefff(x)
! z_effective profile
      implicit real*8 (a-h,o-z)
      common /a0l3/ rh(501),y2dn(501),y2tm(501),y2tmi(501)
      common /a0l4/ con(501),tem(501),temi(501),nspl
      common /a0l5/ zeff(501),y2zeff(501)
      parameter(zero=0.d0,alfa=4.d0,dr=.02d0)
      pa=dabs(x) !#@sav
      if(pa.le.rh(nspl)) then
       call splnt(rh,zeff,y2zeff,nspl,pa,y,dy)
      else
       r=pa-rh(nspl)
       y=zeff(nspl)*dexp(-alfa*(r/dr)**2)
      end if
      zefff=y
      end
c----------------------------------------------------------------
      double precision  function fn(x)
! plasma  density,  cm^-3
      implicit real*8 (a-h,o-z)
      common /a0l3/ rh(501),y2dn(501),y2tm(501),y2tmi(501)
      common /a0l4/ con(501),tem(501),temi(501),nspl
      parameter(zero=0.d0,alfa=4.d0,dr=.02d0)
      pa=dabs(x)
      if(pa.le.rh(nspl)) then
       call splnt(rh,con,y2dn,nspl,pa,y,dy)
      else
       r=pa-rh(nspl)
       y=con(nspl)*dexp(-alfa*(r/dr)**2)
      end if
      fn=y*1.d+13    !cm^-3
      end
c----------------------------------------------------------------
      double precision  function fn1(x,fnp)
! plasma density and its derivative
      implicit real*8 (a-h,o-z)
      common /a0l3/ rh(501),y2dn(501),y2tm(501),y2tmi(501)
      common /a0l4/ con(501),tem(501),temi(501),nspl
      parameter(zero=0.d0,alfa=4.d0,dr=.02d0)
      pa=dabs(x)
      if(pa.le.rh(nspl)) then
       call splnt(rh,con,y2dn,nspl,pa,y,dy)
      else
       call splnt(rh,con,y2dn,nspl,rh(nspl),y1,dy1)
       r=pa-rh(nspl)
       y=rh(nspl)*dexp(-alfa*(r/dr)**2)
       dy=-2.d0*alfa*y*r/dr**2 !corrected
      end if
      fn1=y*1.d+13    !cm^-3
      fnp=dy*1.d+13
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      double precision  function fn2(r,fnp,fnpp)
! plasma density and its first and second derivatives
      implicit real*8 (a-h,o-z)
      common/ne_cheb/chebne(50),chebdne(50),chebddne(50),ncheb
      parameter(zero=0.d0,alfa=4.d0,dr=.02d0)
      x=dabs(r)
      if(x.le.1.d0) then
       y=chebev(zero,1.d0,chebne,ncheb,x)
       dy=chebev(zero,1.d0,chebdne,ncheb,x)
       ddy=chebev(zero,1.d0,chebddne,ncheb,x)
      else
       y1=chebev(zero,1.d0,chebne,ncheb,1.d0)
       s=x-1.d0
       y=y1*dexp(-alfa*(s/dr)**2)
       dy=-2.d0*alfa*y*s/dr**2
       ddy=-2.d0*alfa*y*(1.d0-2.d0*alfa*(s/dr)**2)/dr**2
      end if
      fn2=y    !cm^-3
      fnp=dy
      fnpp=ddy
      end
c----------------------------------------------------------------

      subroutine rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      implicit real*8 (a-h,o-z)
      parameter (nmax=10,fcor=.0666666667d0,
     *    one=1.d0,safety=0.9d0,errcon=6.d-4)
      external derivs
      dimension y(n),dydx(n),yscal(n),ytemp(nmax),ysav(nmax),dysav(nmax)
      pgrow=-0.20d0
      pshrnk=-0.25d0
      xsav=x
      do 11 i=1,n
        ysav(i)=y(i)
        dysav(i)=dydx(i)
11    continue
      h=htry
1     hh=0.5d0*h
      call rk4(ysav,dysav,n,xsav,hh,ytemp,derivs)
      x=xsav+hh
      call derivs(x,ytemp,dydx)
      call rk4(ytemp,dydx,n,x,hh,y,derivs)
      x=xsav+h
      if(x.eq.xsav) then
       write(*,*)' stepsize not significant in rkqc'
       write(*,*)'xsav=',xsav,' h=',h,' htry=',htry
       write(*,88) y,dydx
       pause
      end if
88    format(1x,10(e14.7,1x))

      call rk4(ysav,dysav,n,xsav,h,ytemp,derivs)
      errmax=0.d0
      do 12 i=1,n
        v=ytemp(i)
        ytemp(i)=y(i)-ytemp(i)
        errmax=dmax1(errmax,dabs(ytemp(i)/yscal(i)))
12    continue
      errmax=errmax/eps
      if(errmax.gt.one) then
        h=safety*h*(errmax**pshrnk)
        goto 1
      else
        hdid=h
        if(errmax.gt.errcon)then
          hnext=safety*h*(errmax**pgrow)
        else
          hnext=4.d0*h
        endif
      endif
      do 13 i=1,n
        y(i)=y(i)+ytemp(i)*fcor
13    continue
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rk4(y,dydx,n,x,h,yout,derivs)
      implicit real*8 (a-h,o-z)
      parameter (nmax=10)
      dimension y(n),dydx(n),yout(n),yt(nmax),dyt(nmax),dym(nmax)
      hh=h*0.5d0
      h6=h/6.d0
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs(xh,yt,dyt)
        dv1=dyt(3)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs(xh,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue

      call derivs(x+h,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))
14    continue
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision  function gaussint(f,a,b,r,eps)
      implicit real*8 (a-h,o-z)
      dimension w(12),x(12)
!!      save w,x,const !sav#
      data const /1.0d-12/
      data w
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      data x
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
      delta=const*dabs(a-b)
      gaussint=0d0
      aa=a
    5 y=b-aa
      if(dabs(y) .le. delta) return
    2 bb=aa+y
      c1=0.5d0*(aa+bb)
      c2=c1-aa
      s8=0d0
      s16=0d0
      do 1 i = 1,4
      u=x(i)*c2
    1 s8=s8+w(i)*(f(c1+u,r)+f(c1-u,r))
      do 3 i = 5,12
      u=x(i)*c2
    3 s16=s16+w(i)*(f(c1+u,r)+f(c1-u,r))
      s8=s8*c2
      s16=s16*c2
      if(dabs(s16-s8) .gt. eps*(1d0+dabs(s16))) go to 4
      gaussint=gaussint+s16
      aa=bb
      go to 5
    4 y=0.5d0*y
      if(dabs(y) .gt. delta) go to 2
      write(*,7)
      gaussint=0d0
      return
    7 format(1x,'gaussint ... too high accuracy required')
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine alphas(d,u,j,kmax,g)
      implicit real*8 (a-h,o-z)
      dimension d(50,100),u(50,100),g(50,100)
      common /arr/ dgdu(50,100),kzero(100)
      parameter(zero=0.d0, one=1.d0, tiny=1.d-30)
        km=kzero(j)
        um=u(km,j)
        if(um.ge.one) then
          do k=1,kmax
            if(u(k,j).lt.one) then
             uk=u(k,j)
             uk2=uk**2
             w=dsqrt(one-uk2)
             g(k,j)=w/uk2
             dgdu(k,j)=-one/(w*uk)-2.d0*w/(uk*uk2)
            else
             g(k,j)=zero
             dgdu(k,j)=zero
            end if
          end do
          return
        end if

          do k=1,km
           uk=u(k,j)
           uk2=uk**2
           w=dsqrt(one-uk2)
           g(k,j)=w/uk2
           dgdu(k,j)=-one/(w*uk)-2.d0*w/(uk*uk2)
          end do

       do k=km+1,kmax
         du=u(k,j)-u(k-1,j)
          if(u(k,j).lt.one) then
            beta=u(k,j)*dsqrt(one-u(k,j)**2)
          else
            beta=zero
          end if
         alfa=u(k,j)**3
         g(k,j)=(d(k,j)*g(k-1,j)+beta*du)/(d(k,j)+alfa*du)
           if(d(k,j).ne.zero) then
            dgdu(k,j)=(beta-alfa*g(k,j))/d(k,j)
           else
            dgdu(k,j)=(g(k,j)-g(k-1,j))/du
           end if
       end do
               do k=1,kmax
                if(g(k,j).lt.tiny) g(k,j)=zero
                if(dabs(dgdu(k,j)).lt.tiny) dgdu(k,j)=zero
               end do
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function zatukh(psy,j,u,n)
      implicit real*8 (a-h,o-z)
      dimension u(50,100)
      dimension x(50),y(50),a(50),b(50)
      common /a0befr/ pi,pi2
      common /arr/ dgdu(50,100),kzero(100)
      parameter(zero=0.d0,one=1.d0)
      km=kzero(j)
      um=u(km,j)
      if(um.ge.one) then
       zatukh=zero
       if(psy.lt.one) zatukh=.5d0*pi/psy**3
       return
      end if
      if(psy-um.le.zero.or.u(n,j)-psy.le.zero) then
       zatukh=zero
       return
      end if
      do k=1,n
       x(k)=u(k,j)
       y(k)=dgdu(k,j)
      end do
      i=n-1
      do l=1,n-1
       if(x(l+1)-psy.gt.zero.and.psy-x(l).ge.zero) i=l
      end do
      do k=i,n-1
       b(k)=(y(k+1)-y(k))/(x(k+1)-x(k))
       a(k)=y(k)-b(k)*x(k)
      end do
        s2=dsqrt((x(i+1)-psy)*(x(i+1)+psy))
        ss2=x(i+1)+s2
        sum=a(i)*dlog(psy/ss2)-b(i)*s2
         do k=2,n-i
          s1=dsqrt((x(i+k-1)-psy)*(x(i+k-1)+psy))
          ss1=x(i+k-1)+s1
          s2=dsqrt((x(i+k)-psy)*(x(i+k)+psy))
          ss2=x(i+k)+s2
          sum=sum+a(i+k-1)*dlog(ss1/ss2)+b(i+k-1)*(s1-s2)
         end do
        zatukh=sum
       return
       end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function fdf(x,c,n,df)
      implicit real*8 (a-h,o-z)
      dimension c(n)
      p=c(n)
      dp=0.d0
      do j=n-1,1,-1
        dp=dp*x+p
        p=p*x+c(j)
      end do
      fdf=p
      df=dp
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function fdfddf(x,c,n,df,ddf)
      implicit real*8 (a-h,o-z)
      dimension c(n)
      p=c(n)
      dp=0d0
      ddp=0d0
      do j=n-1,1,-1
        ddp=ddp*x+2d0*dp
        dp=dp*x+p
        p=p*x+c(j)
      end do
      fdfddf=p
      df=dp
      ddf=ddp
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine splne(x,y,n,y2)
      implicit real*8 (a-h,o-z)
      parameter(nn=1001, zero=0d0)
      dimension x(n),y(n),y2(n),u(nn)
      if(n.gt.nn) stop 'n>nn in splne!'
      y2(1)=zero
      u(1)=zero
      do i=2,n-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*y2(i-1)+2.d0
       y2(i)=(sig-1.d0)/p
       u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *     /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
       qn=zero
       un=zero
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
         do k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
         end do
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine splnt(xa,ya,y2a,n,x,y,dy)
      implicit real*8 (a-h,o-z)
      parameter(zero=0.d0)
      dimension xa(n),ya(n),y2a(n)
      klo=1
      khi=n
      do while(khi-klo.gt.1)
       k=(khi+klo)/2
       if(xa(k).gt.x)then
        khi=k
       else
        klo=k
       endif
      end do
      h=xa(khi)-xa(klo)
      if(h.eq.zero) then
       write(*,*)'bad x input in splnt(), x=',x
       write(*,*)'klo=',klo,' kho=',khi
       stop
      end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      aa=a**2
      bb=b**2
      hh=h**2/6d0
      ax=-1d0/h
      bx=-ax
      y=a*ya(klo)+b*ya(khi)+
     *  (a*(aa-1d0)*y2a(klo)+b*(bb-1d0)*y2a(khi))*hh
      dy=ax*ya(klo)+bx*ya(khi)+
     *  ax*((3.d0*aa-1d0)*y2a(klo)-(3.d0*bb-1d0)*y2a(khi))*hh
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dsplnt(xa,ya,y2a,n,x,y,dy,ddy)
      implicit real*8 (a-h,o-z)
      parameter(zero=0.d0)
      dimension xa(n),ya(n),y2a(n)
      klo=1
      khi=n
      do while(khi-klo.gt.1)
       k=(khi+klo)/2
       if(xa(k).gt.x)then
        khi=k
       else
        klo=k
       endif
      end do
      h=xa(khi)-xa(klo)
      if(h.eq.zero) then
       write(*,*)'bad x input in splnt(), x=',x
       write(*,*)'klo=',klo,' kho=',khi
       stop
      end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      aa=a**2
      bb=b**2
      hh=h**2/6d0
      ax=-1d0/h
      bx=-ax
      y=a*ya(klo)+b*ya(khi)+
     *  (a*(aa-1d0)*y2a(klo)+b*(bb-1d0)*y2a(khi))*hh
      dy=ax*ya(klo)+bx*ya(khi)+
     *  ax*((3.d0*aa-1d0)*y2a(klo)-(3.d0*bb-1d0)*y2a(khi))*hh
      ddy=6.d0*ax*ax*(a*y2a(klo)+b*y2a(khi))*hh
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine diff(x,y,n,dy)
      implicit real*8 (a-h,o-z)
      dimension y(*),x(*),dy(*)
        dy(1)=(y(2)-y(1))/(x(2)-x(1))
        do k=2,n-1
         dy(k)=(y(k+1)-y(k-1))/(x(k+1)-x(k-1))
        end do
        dy(n)=(y(n)-y(n-1))/(x(n)-x(n-1))
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function polin(k,x)
      implicit real*8 (a-h,o-z)
       polin=1d0
       if(k.gt.1) polin=x**(k-1)
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function polin1(k,x)
      implicit real*8 (a-h,o-z)
       polin1=x**k
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function polin2(k,x)
      implicit real*8 (a-h,o-z)
       polin2=x**(k+1)
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine approx(x,y,n,f,m,b)
c
c     y(i)=y(x(i))  the data to be approximated
c     n             number of points in the input data
c     m             number of coefficients of decomposition
c                   over base functions "f(k,x)" :
c                          y(x)=sum_1^m [b(k)*f(k,x)]
c     b(i)          found decomposition coefficients
c
      implicit real*8 (a-h,o-z)
      parameter(zero=0.d0, np=20)
      dimension a(np,np),indx(np)
      dimension y(n),x(n),b(*)

       if(m.gt.np) then
         write(*,*)'index error subroutine "approx"'
         return
       end if

         do j=1,m
          do k=1,j
           a(k,j)=zero
            do i=1,n
             a(k,j)=a(k,j)+f(j,x(i))*f(k,x(i))
            end do
          end do
         end do
                 do k=2,m
                  do j=1,k-1
                    a(k,j)=a(j,k)
                  end do
                 end do

          do k=1,m
           b(k)=zero
            do i=1,n
             b(k)=b(k)+y(i)*f(k,x(i))
            end do
          end do

        call ludcmp(a,m,np,indx,d)
        call lubksb(a,m,np,indx,b)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ludcmp(a,n,np,indx,d)
      implicit real*8 (a-h,o-z)
      parameter (nmax=501, tiny=1.d-20, zero=0.d0)
      dimension a(np,np),indx(n),vv(nmax)
      d=1.d0
      do 12 i=1,n
        aamax=zero
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.zero) pause 'singular matrix.'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=zero
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n) then
          if(a(j,j).eq.zero) a(j,j)=tiny
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.zero) a(n,n)=tiny
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lubksb(a,n,np,indx,b)
      implicit real*8 (a-h,o-z)
      parameter(zero=0.d0)
      dimension a(np,np),indx(n),b(n)
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.zero) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
13        continue
        endif
        b(i)=sum/a(i,i)
14    continue
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine source_new(r,out)
      implicit real*8 (a-h,o-z)
      common /asou/ rsou(102),sou(102),npta
      call lock2(rsou,npta,r,klo,khi,ierr)
      if(ierr.ne.0) then
       write(*,*)'lock2 error in source_new'
       write(*,*)'ierr=',ierr,' rho=',r
       stop
      else
       call linf(rsou,sou,r,fout,klo,khi)
       out=dabs(fout)
      end if
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine linf(x,y,t,fout,klo,khi)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
       dout=(y(khi)-y(klo))/(x(khi)-x(klo))
       fout=y(klo)+dout*(t-x(klo))
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lock(xa,n,x,klo,khi,ierr)
      implicit real*8 (a-h,o-z)
      dimension xa(*)
      parameter(tiny=1.d-14)
      klo=0
      khi=0
      dx1=x-xa(1)
      dx2=x-xa(n)
      if(dx1*dx2.ge.tiny) then
       ierr=1
       return
      end if
      ierr=0
      klo=1
      khi=n
      do while(khi-klo.gt.1)
       k=(khi+klo)/2
       if(xa(k).gt.x)then
         khi=k
       else
         klo=k
       endif
      end do
      if(khi.eq.klo) ierr=1
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lock2(xa,n,x,klo,khi,ierr)
      implicit real*8 (a-h,o-z)
      dimension xa(*)
      parameter(zero=0.d0,tiny=1.d-7)
      ierr=0
      klo=0
      khi=0
      dx1=x-xa(1)
      if(dabs(dx1).lt.tiny) then
       klo=1
       khi=2
       return
      else if(dx1.lt.zero) then
       ierr=1
       return
      end if
      dx2=x-xa(n)
      if(dabs(dx2).lt.tiny) then
       klo=n-1
       khi=n
       return
      else if(dx2.gt.zero) then
       ierr=2
       return
      end if
      klo=1
      khi=n
      do while(khi-klo.gt.1)
       k=(khi+klo)/2
       if(xa(k).gt.x)then
         khi=k
       else
         klo=k
       endif
      end do
      if(khi.eq.klo) ierr=3
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sav2008: below this line there are new subroutins and functions
      SUBROUTINE chebft1(a,b,c,n,func)
! Chebyshev fit: Given a function func, lower and upper limits
! of the interval [a,b], and a maximum degree n, this routine 
! computes the n coefficients c(k) such that func(x) approximately =
! SUMM_(k=1)^(k=n)[c(k)*T(k-1)(y)]-c(1)/2, where y and x are related by
! (5.8.10). This routine is to be used with moderately large n 
! (e.g., 30 or 50), the array of cs subsequently to be truncated
! at the smaller value m such that c(m+1) and subsequent elements 
! are negligible. Parameters: Maximum expected value of n, and ð. 
      implicit none
      INTEGER n,NMAX
      DOUBLE PRECISION a,b,c(n),func,PI
      EXTERNAL func
      PARAMETER (NMAX=50, PI=3.141592653589793d0)
      INTEGER j,k
      DOUBLE PRECISION bma,bpa,fac,y,f(NMAX)
      DOUBLE PRECISION sum
      bma=0.5d0*(b-a)
      bpa=0.5d0*(b+a)
      do 11 k=1,n
        y=cos(PI*(k-0.5d0)/n)
        f(k)=func(y*bma+bpa)
11    continue
      fac=2.d0/n
      do 13 j=1,n
        sum=0.d0
        do 12 k=1,n
          sum=sum+f(k)*cos((PI*(j-1))*((k-0.5d0)/n))
12      continue
        c(j)=fac*sum
13    continue
      return
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION chebev(a,b,c,m,x)
! Chebyshev evaluation: All arguments are input. 
! c(1:m) is an array of Chebyshev coefficients, the first m elements 
! of c output from chebft (which must have been called with
! the same a and b). The Chebyshev polynomial evaluated
! and the result is returned as the function value.
      implicit none
      INTEGER m
      DOUBLE PRECISION chebev,a,b,x,c(m)
      INTEGER j
      DOUBLE PRECISION d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.d0) pause 'x not in range in chebev'
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5d0*c(1)
      return
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE chder(a,b,c,cder,n)
! Given a,b,c(1:n), as output from routine chebft(), and given n, 
! the desired degree of approximation (length of c to be used), 
! this routine returns the array cder(1:n), the Chebyshev 
! coefficients of the derivative of the function whose coefficients 
! are c(1:n).
      implicit none
      INTEGER n
      DOUBLE PRECISION a,b,c(n),cder(n)
      INTEGER j
      DOUBLE PRECISION con
      cder(n)=0.d0
      cder(n-1)=2*(n-1)*c(n)
      do 11 j=n-2,1,-1
        cder(j)=cder(j+2)+2*j*c(j+1)
11    continue
      con=2.d0/(b-a)
      do 12 j=1,n
        cder(j)=cder(j)*con
12    continue
      return
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      integer n,nmax
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      external derivs
      parameter (nmax=50)
cu    uses derivs,rkck
      integer i
      double precision errmax,h,htemp,xnew,yerr(nmax),ytemp(nmax)
     *,safety,pgrow,pshrnk,errcon
      parameter (safety=0.9d0,pgrow=-.2d0,pshrnk=-.25d0,errcon=1.89d-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.d0
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.d0)then
        htemp=safety*h*(errmax**pshrnk)
        h=sign(max(abs(htemp),0.1d0*abs(h)),h)
        xnew=x+h
        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.errcon)then
          hnext=safety*h*(errmax**pgrow)
        else
          hnext=5.d0*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rkck(y,dydx,n,x,h,yout,yerr,derivs)
      integer n,nmax
      double precision h,x,dydx(n),y(n),yerr(n),yout(n)
      external derivs
      parameter (nmax=50)
cu    uses derivs
      integer i
      double precision ak2(nmax),ak3(nmax),ak4(nmax),ak5(nmax),ak6(nmax)
     *,ytemp(nmax),a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,b52,b53,
     *b54,b61,b62,b63,b64,b65,c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
      parameter (a2=.2d0,a3=.3d0,a4=.6d0,a5=1.d0,a6=.875d0,b21=.2d0,b31
     *=3.d0/40.d0,
     *b32=9.d0/40.d0,b41=.3d0,b42=-.9d0,b43=1.2d0,b51=-11.d0/54.d0,b52
     *=2.5d0,
     *b53=-70.d0/27.d0,b54=35.d0/27.d0,b61=1631.d0/55296.d0,b62=175.d0
     */512.d0,
     *b63=575.d0/13824.d0,b64=44275.d0/110592.d0,b65=253.d0/4096.d0,c1
     *=37.d0/378.d0,
     *c3=250.d0/621.d0,c4=125.d0/594.d0,c6=512.d0/1771.d0,dc1=c1-2825.d0
     */27648.d0,
     *dc3=c3-18575.d0/48384.d0,dc4=c4-13525.d0/55296.d0,dc5=-277.d0
     */14336.d0,
     *dc6=c6-.25d0)
      do 11 i=1,n
        ytemp(i)=y(i)+b21*h*dydx(i)
11    continue
      call derivs(x+a2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(b31*dydx(i)+b32*ak2(i))
12    continue
      call derivs(x+a3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(b41*dydx(i)+b42*ak2(i)+b43*ak3(i))
13    continue
      call derivs(x+a4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(b51*dydx(i)+b52*ak2(i)+b53*ak3(i)+b54*ak4(i))
14    continue
      call derivs(x+a5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(b61*dydx(i)+b62*ak2(i)+b63*ak3(i)+b64*ak4(i)+
     *b65*ak5(i))
15    continue
      call derivs(x+a6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(c1*dydx(i)+c3*ak3(i)+c4*ak4(i)+c6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(dc1*dydx(i)+dc3*ak3(i)+dc4*ak4(i)+dc5*ak5(i)+dc6*
     *ak6(i))
17    continue
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_unit(iunit)
!! GET_UNIT returns a free FORTRAN unit number.
!! usage:
!!!      call get_unit(iunit)
!!!      if(iunit.ne.0) then
!!!       write(*,*)'iunit=',iunit
!!!      else
!!!       write(*,*)'no free units up to 299'
!!!       pause
!!!       stop
!!!      end if
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 299 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 299, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
      implicit none
      integer ( kind = 4 ) i
      integer ( kind = 4 ) ios
      integer ( kind = 4 ) iunit
      logical lopen

      iunit = 0
      do i = 1, 299
       if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
        inquire ( unit = i, opened = lopen, iostat = ios )
        if ( ios == 0 ) then
         if ( .not. lopen ) then
          iunit = i
          return
         end if
        end if
       end if
      end do
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine integral(ibeg,iend,x,y,fout)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      fout=0.d0
      if(ibeg.eq.iend) return
      znak=1.d0
      n1=ibeg
      n2=iend
      if(n2.lt.n1) then
       znak=-1.d0
       ie=n1
       n1=n2
       n2=ie
      end if
      sum=0.d0
      do i=n1+1,n2
       dx=x(i)-x(i-1)
       dsum=y(i)+y(i-1)
       sum=sum+.5d0*dsum*dx
      end do
      fout=znak*sum
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine gridvel(v1,v2,vmax,cdel,ni1,ni2,ipt1,kpt3,vrj)
      implicit none
      integer ni1,ni2,ipt1,kpt1,kpt2,kpt3,k
      double precision vrj(*),v1,v2,v12,vmax,cdel
      kpt1=ipt1-1
      kpt2=ni1+ni2+1
      do k=1,kpt1  !0<=v<v1
       vrj(k)=dble(k-1)*v1/dble(kpt1)
      end do
      v12=v1+(v2-v1)*cdel
      do k=1,ni1+1 !v1<=v<=v12
       vrj(k+kpt1)=v1+dble(k-1)*(v12-v1)/dble(ni1)
      end do
      do k=2,ni2+1 !!v12<v<=v2
       vrj(k+kpt1+ni1)=v12+dble(k-1)*(v2-v12)/dble(ni2)
      end do     
      do k=1,kpt3  !v2<v<=vmax
       vrj(k+kpt1+kpt2)=v2+dble(k)*(vmax-v2)/dble(kpt3)
      end do
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fsmoth4(x,y,n,ys)
      implicit real*8 (a-h,o-z)
      external polin2
      parameter(zero=0.d0,np=10,imax=601)
      parameter(m0=1,ndp=1)
      ! m0,ndp - parameters of smoothing procedure
      dimension y(n),x(n),ys(n)
      dimension yy(imax),xx(imax)
      dimension coeffs(np),cffs(np)
      dimension dys(imax)
      if(n.gt.imax) stop 'small imax in subroutine fsmoth4()'
        call diff(x,y,n,dys)
         do k=1,n
          ys(k)=y(k)
         end do
       m=m0
       m2=m+2
       id=m+ndp
       nmax=n-id
       xs=x(1)
       do j=1,nmax
         do i=1,id
          xx(i)=x(j+i)-xs
          yy(i)=y(j+i)-ys(j)-dys(j)*xx(i)
         end do
        call approx(xx,yy,id,polin2,m,coeffs)
         cffs(1)=ys(j)
         cffs(2)=dys(j)
          do k=1,m
           cffs(k+2)=coeffs(k)
          end do
        xs=x(j+1)
        ys(j+1)=fdf(xx(1),cffs,m2,dys(j+1))
       end do

      j=nmax+1
1     continue
         jlast=j
         id=n-jlast
         m=id-1
           if(m.eq.0) then
            j=j+1
            xs=x(j)
            ys(j)=fdf(xx(2),cffs,m2,dys(j))
            return
           end if
         m2=m+2
         do i=1,id
          xx(i)=x(jlast+i)-xs
          yy(i)=y(jlast+i)-ys(jlast)-dys(jlast)*xx(i)
         end do
         call approx(xx,yy,id,polin2,m,coeffs)
         cffs(1)=ys(jlast)
         cffs(2)=dys(jlast)
          do k=1,m
           cffs(k+2)=coeffs(k)
          end do
        j=j+1
        xs=x(j)
        ys(j)=fdf(xx(1),cffs,m2,dys(j))
      go to 1
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine interpol(fj)
      implicit none
      integer p,k,i0,i
      real*8 p2,zero
      double precision fj(1002)
      parameter(zero=0.d0,i0=1002)
      p=zero
      do i=1,i0
      if (fj(i).gt.zero) then
      if (p.gt.zero) then
      p2=(fj(i)-fj(p))/(i-p)
      do k=p+1,i-1
      fj(k)=fj(p)+p2*(k-p)
      end do
      end if
      p=i
      end if
      end do
      end
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine view(tview,iview,nnz,ntet) !sav2008
      implicit real*8 (a-h,o-z)
      integer iview  !sav#
      parameter(length=5000000, mpnt=10000)
      dimension dland(length),dcoll(length),perpn(length),dalf(length)
      dimension vel(length),jrad(length),iww(length),tetai(length)
      dimension xnpar(length),izz(length)
      common/agh/xnpar,vel,dland,dcoll,dalf,perpn,tetai,jrad,iww,izz
cc      dimension an1(length),an2(length)
cc      common /xn1xn2/ an1,an2
      dimension mbeg(mpnt),mend(mpnt),mbad(mpnt),rbeg(mpnt) !sav2008
      dimension tetbeg(mpnt),xnrbeg(mpnt),xmbeg(mpnt),yn3beg(mpnt)
      common/viewdat/mbeg,mend,mbad,rbeg,tetbeg,xnrbeg,xmbeg,yn3beg
      common /a0k/ cdl(10),cly(10),cgm(10),cmy(10),ncoef
      common /a0i5/ vperp(50,100),cnstal,zza,zze,valfa,kv
      common /a0ef1/ r0,z0,rm,cltn
      common /bcef/ ynz,ynpopq
      common /a0befr/ pi,pi2
      common /a0a1/ ynzm(1001),pm(1001),nmaxm(4)
      common /a0a2/ tet1,tet2
      common /a0ab/ nr
      common /vth/ vthc(length),poloidn(length)
      common/b0/ itend0
      real*8 vthcg,npoli
      common /a0ghp/ vlf,vrt,dflf,dfrt
      integer unit_bias
      parameter (unit_bias = 10)
      parameter(m=7,pleft=1.d-10) !m may be chaged together with name(m)
      dimension ptt(m),pll(m),pcc(m),paa(m)
      dimension pt_c(m),pl_c(m),pc_c(m),pa_c(m)
      character(40) fname
      character*40 name(m)
      save name !sav#
      data name/'lhcd/out/1.dat','lhcd/out/2.dat','lhcd/out/3.dat'
     &,'lhcd/out/4.dat','lhcd/out/5.dat'
     &,'lhcd/out/rest.dat','lhcd/out/traj.dat'/
      parameter(zero=0.d0, clt=3.d10)
      if(iview.eq.0) return
      print *, 'view_time=',tview
      print *, name(m)
      write(fname,'("lhcd/out/traj", f17.16,".dat")') tview
      print *, fname
      name(m) = fname
      print *, name(m)

      htet=zero
      h=1d0/dble(nr+1)
      if(ntet.ne.1) htet=(tet2-tet1)/(ntet-1)
      open(1,file='lhcd/out/lcms.dat')
      write(1,*)'     R(m)            Z(m)'
      write(1,*)
      xr=1.d0
      xdl=fdf(xr,cdl,ncoef,xdlp)
      xly=fdf(xr,cly,ncoef,xlyp)
      xgm=fdf(xr,cgm,ncoef,xgmp)
      do i=1,101
       th=dble(i-1)*pi2/dble(100)
       cotet=dcos(th)
       sitet=dsin(th)
       xx=-xdl+xr*cotet-xgm*sitet**2
       zz=xr*xly*sitet
       x=(r0+rm*xx)/1d2
       z=(z0+rm*zz)/1d2
       write(1,5) x,z
      end do
      close(1)
      open(1,file='lhcd/out/npar_crit.dat')
      write(1,*)'  Npar_crit=sqrt(50/Te(keV))'
      write(1,*)
      write(1,*)'   rho         Npar_strong absorption'
      write(1,*)
      do i=1,101
       xr=dble(i-1)/dble(100)
       tmp=ft(xr)/0.16d-8  !Te,  KeV
       parn_c=dsqrt(50d0/tmp)
       write(1,5) xr,parn_c
      end do
      close(1)

      do n=1,m
       open(n+unit_bias,file=name(n))
       write(n+unit_bias,*)
       close(n+unit_bias)
       open(n+unit_bias,file=name(n))
       write(n+unit_bias,3)
       write(n+unit_bias,*)
       ptt(n)=zero
       pll(n)=zero
       pcc(n)=zero
       paa(n)=zero
      end do
      ntraj=0 !sav2008
      do itr=1,nnz*ntet !sav2008
       pow=1.d0
       pl=zero
       pc=zero
       pa=zero
       pdec1=zero
       pdec1z=zero
       pdec3=zero
       pdec3z=zero
       pdecv=zero
       pintld=zero
       pintal=zero
       jrc=nr+1
       jznak=-1
       nturn=1
       if(mbad(itr).eq.0) then 
        ntraj=ntraj+1
        ib=mbeg(itr)
        ie=mend(itr)
10      continue
        do i=ib,ie
         v=vel(i)
         jr=jrad(i)
         refr=perpn(i)
         npoli=poloidn(i)
         ifast=iww(i)
         vthcg=vthc(i)
         idir=izz(i)
         dek3=zero
         th=tetai(i)
         parn=xnpar(i)
cc         xn1=an1(i)
cc         xn2=an2(i)
         if(itend0.gt.0) then
          argum=clt/(refr*valfa)
          dek3=zatukh(argum,abs(jr),vperp,kv)
         end if
!!!old variant
!!!         call raspr(v,abs(jr),iv,df)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call distr(v,abs(jr),iv,df)
         if(jr.lt.0) then    !case of turn
          jr=-jr
!variant          pintld=-dland(i)*df
!!          pintld=-dland(i)*(dflf+dfrt)/2d0
          pintld=dabs(dland(i)*(dflf+dfrt)/2d0)
          pdec2=dexp(-2d0*dcoll(i))
          pintal=dabs(dalf(i)*dek3)
         else
          pdec2=dcoll(i)
          pdecv=dland(i)
!!          pdec1=-pdecv*df
          pdec1=dabs(pdecv*df)
          pdec3=dabs(dalf(i)*dek3)
          pintld=(pdec1+pdec1z)/2d0*h
          pintal=(pdec3+pdec3z)/2d0*h
          pdec1z=pdec1
          pdec3z=pdec3
         end if
         powpr=pow
         powd=pow*dexp(-2d0*pintld)
         powcol=powd*pdec2
         powal=powcol*dexp(-2d0*pintal)
         pow=powal
         pil=pintld
         pic=.5d0*dabs(dlog(pdec2))
         pia=pintal
         pt=1.d0-pow  !total absorbed power
         denom=pil+pic+pia
         powdamped=1.d0-dexp(-2.d0*denom)
         domin=powpr*powdamped
         if(denom.ne.zero) then
          fff=domin/denom
          pl=pl+dabs(pil*fff)  !el. Landau absorbed power
          pc=pc+dabs(pic*fff)  !el. collisions absorbed power
          pa=pa+dabs(pia*fff)  !alpha Landau absorbed power
         end if
         xr=h*dble(jr)
         cotet=dcos(th)
         sitet=dsin(th)
         xdl=fdf(xr,cdl,ncoef,xdlp)
         xly=fdf(xr,cly,ncoef,xlyp)
         xgm=fdf(xr,cgm,ncoef,xgmp)
         xx=-xdl+xr*cotet-xgm*sitet**2
         zz=xr*xly*sitet
         x=(r0+rm*xx)/1d2
         z=(z0+rm*zz)/1d2 
         jdlt=jr-jrc
         jrc=jr
         if(jdlt*jznak.lt.0.and.nturn.lt.m-1) then
          nturn=nturn+1
          jznak=-jznak
         end if
         mn = nturn + unit_bias
         write(mn, 7) x,z,xr,th,parn,npoli,pt,pl,pc,pa,ifast,idir,itr
         mm = m + unit_bias
         write(mm, 7) x,z,xr,th,parn,npoli,pt,pl,pc,vthcg,ifast,idir,itr
         do n=m,nturn,-1
          pt_c(n)=pt
          pl_c(n)=pl
          pc_c(n)=pc
          pa_c(n)=pa
         end do
         if(pt.ge.1d0-pleft) go to 11 !maximal absorbed power along a ray
        end do
        jchek=jrad(ie+1)
        if(jchek.ne.0) then  !continue this trajectory
         ib=idnint(dland(ie+1))
         ie=idnint(dcoll(ie+1))
         goto10
        end if
11      continue
        do n=1,m
         ptt(n)=ptt(n)+pt_c(n)
         pll(n)=pll(n)+pl_c(n)
         pcc(n)=pcc(n)+pc_c(n)
         paa(n)=paa(n)+pa_c(n)
        end do
        if(itr.lt.nnz*ntet) then
         write(m+unit_bias,*)
         write(nturn+unit_bias,*)
        end if
       end if
      end do
      do n=1,m
       close(n+unit_bias)
      end do
      do n=1,m
       ptt(n)=ptt(n)/dble(ntraj)
       pll(n)=pll(n)/dble(ntraj)
       pcc(n)=pcc(n)/dble(ntraj)
       paa(n)=paa(n)/dble(ntraj)
      end do
      open(1,file='lhcd/out/info_traj.dat')
      write(1,20) tview
      write(1,*)
      do n=1,m-1
       if(n.lt.m-1) then
        write(1,8) n,ptt(n),pll(n),pcc(n),paa(n)
       else
        write(1,9) ptt(n),pll(n),pcc(n),paa(n)
       end if
       write(1,*)
      end do
      write(1,*)
      write(1,1)
      write(1,*)
      itr=0
      do itet=1,ntet
       tetin=tet1+htet*(itet-1)
       do inz=1,nnz
        itr=itr+1
        write(1,6) itr,mbad(itr),tetin,ynzm(inz),rbeg(itr)
       end do
      end do
      close(1)
      open(1,file='lhcd/out/absorp.dat')
      write(1,2)
      write(1,*)
      do n=1,m-1
       if(n.eq.1) then
        dpt=ptt(n)
        dpl=pll(n)
        dpc=pcc(n)
        dpa=paa(n)
       else
        dpt=ptt(n)-ptt(n-1)
        dpl=pll(n)-pll(n-1)
        dpc=pcc(n)-pcc(n-1)
        dpa=paa(n)-paa(n-1)
       end if
       write(1,4) n,ptt(n),pll(n),pcc(n),paa(n)
     &             ,dpt,dpl,dpc,dpa
      end do
      close(1)

1     format(2x,'N_traj',3x,'mbad',6x,'theta',9x,'Npar',9x,'rho_start')
2     format('R_pass',4x,'Ptot',6x,'Pland',6x,'Pcoll',8x,'Pa'
     &,7x,'dPtot',6x,'dPland',5x,'dPcoll',6x,'dPa')
3     format(5x,'R',10x,'Z',11x,'rho',8x,'theta',7x,'N_par',7x,'N_pol'
     &,6x,'P_tot',7x,'P_land',6x,'P_coll',6x,'vth'
     &,4x,'slow=1',4x,'out=1',2x,'N_traj',6x)
4     format(i3,5x,8(f6.3,5x))
5     format(6(e13.6,3x))
6     format(2(i6,2x),4(e13.6,1x))
7     format(10(e11.4,1x),i5,2x,i5,2x,i5)
8     format('after radial pass=',i3,2x,' P_tot=',f6.3,2x
     &    ,' P_land=',f5.3,2x,' P_coll=',f6.3,2x,' P_a=',f6.3)
9     format('Total passes:           P_tot=',f6.3,2x
     & ,' P_land=',f5.3,2x,' P_coll=',f6.3,2x,' P_a=',f6.3)
20    format('written time slice (seconds) =',f9.3)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

