      SUBROUTINE OPENOU(noutpt,nttyo,ufiln,uform,nrecl,ilu)

c     This subroutine opens an output type file. If a file of the
c     same name already exists, it is first destroyed. An unused
c     logical unit number is obtained.

c     This subroutine is called by:

c       H2OI95

c-----------------------------------------------------------------------

c     Input:

c       noutpt = the unit number of the output file
c       nttyo  = the unit number of the screen file
c       ufiln  = the name of the file to open
c       uform  = the file format, 'formatted' or 'unformatted'
c       nrecl  = the record length (number of characters per line)

c     Output:

c       ilu     = the logical unit number of opened file

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      INTEGER noutpt,nttyo

      INTEGER ilu,nrecl

c     Note- the character length for the variables ufiln and uform
c     can not be exactly specified.

      CHARACTER(LEN=*) ufiln,uform

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER j2,nerr

      INTEGER ILNOBL

      LOGICAL qex

      CHARACTER(LEN=8) ustat

c-----------------------------------------------------------------------

      j2 = ILNOBL(ufiln)

c     See if a file of the same name already exists. If so,
c     destroy it. This makes the logical unit number available.
c     If a file of the same name does not exist, get the next
c     available logical unit number.

      INQUIRE(file=ufiln,exist=qex)
      IF (qex) THEN
        ustat = 'old'
        CALL GETLU(ilu,nerr)
        IF (nerr .NE. 0) THEN

          IF (noutpt .GT. 0) WRITE (noutpt,1000) ustat,ufiln(1:j2)

          WRITE (nttyo,1000) ustat,ufiln(1:j2)

 1000     FORMAT(/6x,'ERROR (OPENOU): No logical unit number',
     $    /9x,'is available to open the ',a3,' file "',a,'".')

          STOP
        ENDIF
        OPEN(ilu,file=ufiln,status=ustat,err=10)
        CLOSE(ilu,status='delete',err=15)
      ELSE
        ustat = 'new'
        CALL GETLU(ilu,nerr)
        IF (nerr .NE. 0) THEN
          IF (noutpt .GT. 0) WRITE (noutpt,1000) ustat,ufiln(1:j2)
          WRITE (nttyo,1000) ustat,ufiln(1:j2)
          STOP
        ENDIF
      ENDIF

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

c     Open the new file.

      IF (nrecl .GT. 0) THEN

c       Use the specified record length.

        OPEN(ilu,file=ufiln,form=uform,status='new',recl=nrecl,err=10)
      ELSE

c       The record length is not specified. Open the file at the
c       default record length.

        OPEN(ilu,file=ufiln,form=uform,status='new',err=10)
      ENDIF
      GO TO 999

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

   10 IF (noutpt .GT. 0) WRITE (noutpt,1010) ustat,ufiln(1:j2)

      WRITE (nttyo ,1010) ustat,ufiln(1:j2)

 1010 FORMAT(/6x,"ERROR (OPENOU): Can't open the ",a3,' copy',
     $ /9x,'of the file "',a,'".')

      STOP

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

   15 IF (noutpt .GT. 0) WRITE (noutpt,1020) ufiln(1:j2)

      WRITE (nttyo,1020) ufiln(1:j2)

 1020 FORMAT(/6x,"ERROR (OPENOU): Can't delete the old copy",
     $ /9x,'of the file "',a,'".')

      STOP

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  999 CONTINUE

      END