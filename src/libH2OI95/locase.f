      SUBROUTINE LOCASE(ustr)

c     This subroutine converts a string from upper case to lower case.

c     This subroutine is called by:

c       Any

c-----------------------------------------------------------------------

c     Input:

c       ustr   = the input string variable

c     Output:

c       ustr   = the output string variable

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      CHARACTER(LEN=*) ustr

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER idel,j,nchars

      CHARACTER(LEN=1) u1

c-----------------------------------------------------------------------

      idel = ichar('A') - ichar('a')
      IF (idel .NE. 0) THEN
        nchars = len(ustr)
        DO j = 1,nchars
          u1 = ustr(j:j)
          IF (u1.GE.'A' .AND. u1.LE.'Z')
     $    ustr(j:j) = char(ichar(u1) - idel)
        ENDDO
      ENDIF

      END