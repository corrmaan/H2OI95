      SUBROUTINE OPENIN(noutpt,nttyo,ufiln,uform,ilu)

c     This subroutine opens an input type file. The file must already
c     exist. An unused logical unit number is obtained.

c     This subroutine is called by:

c       H2OI95

c-----------------------------------------------------------------------

c     Input:

c       noutpt = the unit number of the output file
c       nttyo  = the unit number of the screen file
c       ufiln  = the name of the file to open
c       uform  = the file format, 'formatted' or 'unformatted'

c     Output:

c       ilu     = the logical unit number of opened file

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

c     Note- the character length for the variables ufiln and uform
c     can not be exactly specified.

      INTEGER ilu,noutpt,nttyo

      CHARACTER(LEN=*) ufiln,uform

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER j2,nerr

      INTEGER ILNOBL

      LOGICAL qex

      CHARACTER(LEN=8) uformo

c-----------------------------------------------------------------------

      j2 = ILNOBL(ufiln)

c     Check to make sure the file exists.

      INQUIRE(file=ufiln,exist=qex,formatted=uformo)
      IF (.NOT.qex) THEN
        IF (noutpt .GT. 0) WRITE (noutpt,1000) ufiln(1:j2)

        WRITE (nttyo,1000) ufiln(1:j2)

 1000   FORMAT(/6x,'ERROR (OPENIN): The file "',a,'"',
     $  /9x,'does not exist.')

        STOP
      ENDIF

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

c     Check the file format.

      IF (uform.EQ.'formatted' .AND. uformo.EQ.'no') THEN

        IF (noutpt .GT. 0) WRITE (noutpt,1010) ufiln(1:j2)

        WRITE (nttyo,1010) ufiln(1:j2)

 1010   FORMAT(/6x,'ERROR (OPENIN): The file "',a,'"',
     $  /9x,'should be formatted, but it is not.')

        STOP
      ENDIF

      IF (uform.EQ.'unformatted' .AND .uformo.EQ.'yes') THEN

        IF (noutpt .GT. 0) WRITE (noutpt,1020) ufiln(1:j2)

        WRITE (nttyo,1020) ufiln(1:j2)

 1020   FORMAT(/6x,'ERROR (OPENIN): The file "',a,'"',
     $  /9x,'should be unformatted, but it is not.')

        STOP
      ENDIF

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

c     Get the next available logical unit number.

      CALL GETLU(ilu,nerr)
      IF (nerr .NE. 0) THEN
        IF (noutpt .GT. 0) WRITE (noutpt,1050) ufiln(1:j2)

        WRITE (nttyo,1050) ufiln(1:j2)

 1050   FORMAT(/6x,'ERROR (OPENIN): No logical unit number',
     $  /9x,'is available for the file "',a,'".')

        STOP
      ENDIF

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      OPEN(ilu,file=ufiln,form=uform,status='old',err=10)
      GO TO 999

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

   10 CONTINUE
      IF (noutpt .GT. 0) WRITE (noutpt,1060) ufiln(1:j2)

      WRITE (nttyo,1060) ufiln(1:j2)

 1060 FORMAT(/6x,"ERROR (OPENIN): Can't open the file",
     $ /9x,'"',a,'".')

      STOP

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  999 CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

