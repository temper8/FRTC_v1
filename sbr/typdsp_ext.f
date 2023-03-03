C Modules >>> TYPDSP, PCOVA, UF1DWA, UF2DWA, OPENAP, OPENRD, OPENWT <<<
C >>> SETFRM, AFRAME, DNSTR, UPSTR,  TIMEDT, PUTXY,  NEGA,   SETFNA <<< 
C >>> SKIPM  <<< 
C======================================================================|
	subroutine	TYPDSP_EXT(NCH,YN,ITIMES,TTOUT,TOUT)
C NCH= 5 - terminal, 0 - file (old format), 1 - file (new format)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
	integer	ITIMES
	double precision	YN,YQ,TTOUT(ITIMES),TOUT(ITIMES,NRW)
	integer NCH,NP1,ITBE,ITEND,ITEN,IERR,MODEX,NLINSC,NVAR,NNJ,JN0
	integer JBE,JEND,J,JEN,JJ,J1,JLR,KILLBL,length,WarningColor
	character*6 CH6,STRMN*40,CONN(6)
	character(132) FNAME, STRI, tmp
	logical	EXI
	save	STRMN,CONN,JN0,NLINSC,WarningColor
	data	STRMN/' R=     a=     B=     I=     q=     <n>='/
     +	CONN/' CF   ',' CV   ',' CH   ',' CCD  ',' CBND ',' CRAD '/
     +	JN0/0/ NLINSC /50/ WarningColor/30/
C NLINSC - maximum line number 
	MODEX = XOUT+.49
C MODEX = 0	[0,AB]  against "a"
C MODEX = 1	[0,ABC] against "a"
C MODEX = 2	[0,ROC] against "rho"
C MODEX = 3	[FP(1),FP(NA1)] against "psi"
C otherwise	Unknown option => MODEX=0
	NP1 = NAB
	if (MODEX .ge. 1 .and. MODEX .le. 3 .or. MOD10 .eq. 3)
     +	NP1 = NA1

	if(NCH.eq.5)	goto	300
	if(NCH.ne.0 .and. NCH.ne.1)	return
C Write output to a file:
	J = length(RDNAME)
	J1 = length(EQNAME)
	!FNAME=AWD(1:length(AWD))//'dat/XData'
	!JEND = KILLBL(FNAME,132)
	!call	SETFNA_EXT(FNAME,JEND)
	!FNAME=REPEAT('*',132)
	write(FNAME,'("dat/XData/", f9.7,".dat")')  time
C	write(*,*)'Returned:  "',FNAME(1:JEND),'"',JEND
	call colovm(WarningColor)
	JEND = KILLBL(FNAME,132)
	print '(">>>  XData are written into the file: ", A)', FNAME(1:JEND)
	JLR = XWH-125
C	call textvm(JN0,JLR,STRI,37+JEND)
	call	OPENWT(7,FNAME,0,IERR)
	if(IERR.gt.0)	then
	   write(*,*)'>>> TYPDSP: Output file error'
	   stop
	endif
C Creating UPSTRI
	STRI=XLINE1(2:17)
	STRI(17:)=STRMN
	call FMTF4(STRI(20:23),RTOR)
	call FMTF4(STRI(27:30),ABC)
	call FMTF4(STRI(34:37),BTOR)
	call FMTF4(STRI(41:44),IPL)
C Triangularity corrected MHD q (accoding to ITER guidelines)
C	YQ	=ELON(NA)**2
C	YD	=TRIA(NA)
C	YQ=(1.+YQ*(1.+YD**2*(2.-1.2*YD)))/(MU(NA)*(1.+YQ))
	YQ=1./MU(NA)
	call FMTF4(STRI(48:51),YQ)
	call FMTF4(STRI(57:60),YN)
	write(STRI(62:78),103)TIME
	call FMTF4(STRI(79:82),1000.*TAU)
 102	format('   Time',16(3X,1A4))
 103	format('Time=',1F6.5,' dt=')
 104	format(1X,1A120)
 	 write(7,*) 'Ext format'
	 write(7,104)STRI
	 if(NCH.eq.1)	goto	400

C Writing radial data
	write(7,*)
	write(7,*) 'Writing radial data'
	if (MOD10.le.5)			then
	    JBE=1
	    JEND=16
 3	    JEN=MIN0(NTOUT,JEND)
 		write(7,'(100A10)'), "Time", (NAMET(J),J=JBE,JEN)
 		write(7,'(100f10.1)') TIME, (TOUT(LTOUT,j), j=JBE,JEN )

	    if(JEN.eq.NTOUT)	go to 4
	    JBE=JEN+1
	    JEND=JEN+16
		go to 3

 4	    write(7,*)
 		JBE=1
	    JEND=16
 1	    JEN=MIN0(NROUT,JEND)

 		 write(7,'(99A20))') 'a',(NAMER(J),J=1,NROUT)

		do  j=1, NP1 
			SELECT CASE (MODEX) 
			CASE (0:1)		
				write(7,'(100ES22.14)') AMETR(j), (ROUT(j,jj), jj=1,NROUT )
			CASE (2)
				write(7,'(100ES22.14)') RHO(j), (ROUT(j,jj), jj=1,NROUT )				
			CASE (3)
				write(7,'(100ES22.14)') FP(j), (ROUT(j,jj), jj=1,NROUT )					
			CASE DEFAULT
				write(7,'(100ES22.14)') AMETR(j), (ROUT(j,jj), jj=1,NROUT )
		    END SELECT
			
		enddo
		write(7,*)
	endif

C Writing time data
	if(MOD10.eq.6)	then
		JBE=1
		JEND=16
 6		JEN=MIN0(NTOUT,JEND)
		write(7,102)	(NAMET(J),J=JBE,JEN)
		do	15	J1=1,LTOUT-1
		STRI=' '
		call FMTXF5(STRI(1:5),TTOUT(J1))
		do	14	J=JBE,JEN
		JJ=7*(J-JBE)+8
 14		call FMTXF5(STRI(JJ:JJ+5),TOUT(J1,J))
 15		write(7,104)STRI
		if(JEN.eq.NTOUT)	go to 2
		JBE=JEN+1
		JEND=JEN+16
					go to 6
					endif
 400	continue
C Writing constants
 2	write(7,'(10X,1A80)')RUNID
	J1=0
	do	93 JEN=1,100
	   STRI=' '
	   do	J=1,16
	      J1=J1+1
	      if (J1.gt.NCFNAM)	go to 94
	      call FMTXF5(CH6,CONSTF(J1))
	      JJ=7*(J-1)+1
	      STRI(JJ:JJ+5)=CH6
	   enddo
 93	write(7,101) CONN(JEN),STRI
 94	write(7,101) CONN(JEN),STRI
 101	format(1X,1A6,1A111)
	if(NCH.eq.1)	goto	401
	close(7)
	return

 401	continue
C Writing radial data
	if(MOD10.le.5)			then
	JBE=1
	JEND=16
 402	JEN=MIN0(NTOUT,JEND)
	write(7,'(3X,"Time",16(3X,1A4))')(NAMET(J),J=JBE,JEN)
	STRI=' '
	call FMTXF4(STRI(1:5),TIME)
	do	J=JBE,JEN
	   JJ=7*(J-JBE)+8
	   call FMTXF5(STRI(JJ:JJ+5),TOUT(LTOUT,J))
	enddo
	write(7,104)STRI
	if(JEN.EQ.NTOUT)	go to 403
	JBE=JEN+1
	JEND=JEN+16
				go to 402
 403	JBE=1
	JEND=NRW
	JEN=MIN0(NROUT,JEND)
C Different options for radial variable 
	if (MODEX .eq. 2)	then
	   write(7,'(8X,"rho ",64(8X,1A4))')(NAMER(J),J=JBE,JEN)
	   do	j=1,NP1
		write(7,408)RHO(j),(ROUT(J,JJ),JJ=JBE,JEN)
	   enddo
	elseif (MODEX .eq. 3 .or. MOD10 .eq. 3)	then
	   write(7,'(8X,"psi ",64(8X,1A4))')(NAMER(J),J=JBE,JEN)
	   do	j=1,NP1
		write(7,408)FP(j),(ROUT(J,JJ),JJ=JBE,JEN)
	   enddo
	else
	   write(7,'(8X,"a   ",64(8X,1A4))')(NAMER(J),J=JBE,JEN)
	   do	j=1,NP1
		write(7,408)AMETR(j),(ROUT(J,JJ),JJ=JBE,JEN)
	   enddo
	endif
	endif

 408    format(1PE12.3,64(1PE12.3))

C Writing time data
	if(MOD10.ne.6)	goto	409
	STRI=' '
	write(7,104)STRI
	write(7,104)STRI
	write(7,104)STRI
	write(7,104)STRI
	JBE = 1
	JEND=MIN(NTOUT,NRW)
 404	continue
	JEN = JBE+7
	write(7,'(8X,"Time",64(8X,1A4))')(NAMET(J),J=JBE,JEN)
	do	J1=1,LTOUT-1
	   STRI=' '
	   call FMTXF5(STRI(1:5),TTOUT(J1))
	   write(7,408)TTOUT(J1),(TOUT(J1,J),J=JBE,min(JEN,JEND))
	enddo
	JBE = JBE+8
	if (JBE.lt.JEND)	goto	404
 409	close(7)
	return

 300	continue
C Output to the terminal
	if(MOD10.gt.5)	goto	305
	JBE=1
	JEND=16
 301	JEN=MIN0(NROUT,JEND)
	write(STRI,302)(NAMER(J),J=JBE,JEN)
 302	format(16(1X,1A4))
	call prntxt(STRI)
	do 304 J=1,NP1
	STRI=' '
	do	303	JJ=JBE,JEN
	J1=5*(JJ-JBE+1)-4
 303	call FMTXF4(STRI(J1:J1+4),ROUT(J,JJ))
	NNJ=9*J+3
	call prntxt(STRI)
 304	continue
	if(JEN.eq.NROUT)	return
	JBE=JEN+1
	JEND=JEN+16
	GO TO 301
 305	if(MOD10.gt.7)	return
	JBE=1
	JEND=15
 306	JEN=MIN(NTOUT,JEND)
	ITBE=1
	ITEND=NLINSC
 307	ITEN=MIN(LTOUT-1,ITEND)
	write(STRI,308)	(NAMET(J),J=JBE,JEN)
 308	format(1X,'Time',15(1X,1A4))
	call	prntxt(STRI)
	do	315	J1=ITBE,ITEN
	STRI=' '
	call FMTXF4(STRI(1:5),TTOUT(J1))
	do	314	J=JBE,JEN
	JJ=5*(J-JBE)+6
 314	call FMTXF4(STRI(JJ:JJ+4),TOUT(J1,J))
	NNJ=9*(J1-ITBE)+12
	call prntxt(STRI)
 315	continue
	if(ITEN.eq.LTOUT-1)	go to 316
	ITBE=ITEN
	ITEND=ITEN+NLINSC-1
				go to 307
 316	if(JEN.eq.NTOUT)	return
	JBE=JEN+1
	JEND=JEN+15
				GO TO 306
	END

C======================================================================|
	subroutine	SETFNA_EXT(FNAME,NLEN)
C FNAME - input name (without blanks) is appended with an extension.
C	  The extension is the ordinal number of the file
C NLEN  - returns length of the new name
C
		implicit none
		integer	killbl,nvar,nlen,length,j
		character FNAME*(*)
		logical	EXI
		NLEN = length(FNAME)
		NVAR = 0
1		NVAR = NVAR+1
		write(FNAME(NLEN+1:),'(1A1,I0.3)')'.',NVAR
		j = killbl(FNAME,NLEN+4)
C	write(*,*)'Inquired:  "',FNAME(1:j),'"',j
		inquire(FILE=FNAME,EXIST=EXI)
		if(EXI)		goto 1
		NLEN = j
		return
		end
C   ======================================================================|