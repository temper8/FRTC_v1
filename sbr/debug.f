      subroutine debug
      implicit none
      double precision LINEAV
      double precision CHORDN
      integer NTRUN
      include	'for/parameter.inc'
      include	'for/const.inc'
      include	'for/status.inc'
      include	'for/outcmn.inc'
      include	'for/timeoutput.inc'      
      !print *, "time=", time
      CHORDN = LINEAV(1)
      NTRUN = NTIMES
      call TYPDSP_EXT(0,CHORDN,NTRUN,TTOUT,TOUT)
      end subroutine