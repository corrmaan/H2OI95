      SUBROUTINE GETLU(nlu,nerr)

c     This subroutine finds a currently unused unit number.

c     This subroutine is called by:

c       OPENIN.f
c       OPENOU.f

c-----------------------------------------------------------------------

c     Input:

c       None

c     Output:

c       nlu    = first currently unused unit number
c       nerr   = error flag:
c                  = 0   Okay
c                  = 1   Error

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      INTEGER nlu,nerr

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER iumax,iumin

      LOGICAL qopen

c-----------------------------------------------------------------------

      data iumax,iumin /40,0/

c-----------------------------------------------------------------------

      nerr = 0

c     Loop through all valid file numbers, beginning with the largest.

      DO nlu = iumax,iumin,-1
        INQUIRE(unit=nlu,opened=qopen)
        IF (.NOT.qopen) GO TO 999
      ENDDO

      nerr = 1

  999 CONTINUE

      END