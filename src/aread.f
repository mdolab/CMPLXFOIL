
      SUBROUTINE AREAD(LU,FNAME,NMAX,X,Y,N,NAME,ISPARS,ITYPE,INFO)
      DIMENSION X(NMAX), Y(NMAX)
      CHARACTER*(*) FNAME
      CHARACTER*(*) NAME
      CHARACTER*(*) ISPARS
C--------------------------------------------------------
C     Reads in several types of airfoil coordinate file.
C
C  Input:
C       LU      logical unit to use for reading
C       FNAME   name of coordinate file to be read,
C               if FNAME(1:1).eq.' ', unit LU is assumed
C               to be already open
C       INFO   0 keep quiet
C              1 print info on airfoil
C  Output:
C       X,Y     coordinates
C       N       number of X,Y coordinates
C       NAME    character name string        (if ITYPE > 1)
C       ISPARS  ISES/MSES domain-size string (if ITYPE > 2)
C       ITYPE returns type of file:
C           0  None.  Read error occurred.
C           1  Generic.
C           2  Labeled generic.
C           3  MSES single element.
C           4  MSES multi-element.
C--------------------------------------------------------
C*****************************************************************
CModified 2.24.06 by D. Berkenstock to suppress text output
C*****************************************************************
C
      CHARACTER*80 LINE1,LINE2,LINE
      LOGICAL LOPEN, ERROR
      DIMENSION A(10)
C
      IEL = 0
      NEL = 0
C
C---- assume read error will occur
      ITYPE = 0
C
      LOPEN = FNAME(1:1) .NE. ' '
      IF(LOPEN) OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=98)
C
 11   READ(LU,1000,END=99,ERR=98) LINE1
      IF(INDEX('#!',LINE1(1:1)) .NE. 0) GO TO 11
C
 12   READ(LU,1000,END=99) LINE2
      IF(INDEX('#!',LINE2(1:1)) .NE. 0) GO TO 12
C
      I = 1
C
C---- try to read numbers from first line
      NA = 10
      CALL GETFLT(LINE1,A,NA,ERROR)
      IF(ERROR .OR. NA.LT.2) THEN
        NAME = LINE1
      ELSE
        IF(INFO.GT.0) THEN
C         WRITE(*,*)
C         WRITE(*,*) 'Plain airfoil file'
        ENDIF
        ITYPE = 1
        REWIND(LU)
        GO TO 50
      ENDIF
C
C---- try to read numbers from second line
      NA = 10
      CALL GETFLT(LINE2,A,NA,ERROR)
      IF(ERROR .OR. NA.GE.4) THEN
        ISPARS = LINE2
      ELSE
        NAME = LINE1
        IF(INFO.GT.0) THEN
C         WRITE(*,*)
C         WRITE(*,*) 'Labeled airfoil file.  Name:  ', NAME
        ENDIF
        ITYPE = 2
        REWIND(LU)
        READ(LU,1000,END=99) LINE1
        GO TO 50
      ENDIF
C
      IF(INFO.GT.0) THEN
C       WRITE(*,*)
C       WRITE(*,*) 'MSES airfoil file.  Name:  ', NAME
      ENDIF
      ITYPE = 3
C
C---- read each element until 999.0 or end of file is encountered
   50 NEL = NEL + 1
      DO 55 I=1, NMAX
 51     READ(LU,1000,END=60) LINE
C
C------ skip comment line
        IF(INDEX('#!',LINE(1:1)) .NE. 0) GO TO 51
C
        NA = 2
        CALL GETFLT(LINE,A,NA,ERROR)
        IF(ERROR) GO TO 99
C
C------ skip line without at least two numbers
        IF(NA.LT.2) GO TO 51
C
        X(I) = A(1)
        Y(I) = A(2)
C
        IF (X(I) .EQ. 999.0) THEN
C-------- if this is the element we want, just exit
          IF(IEL .EQ. NEL) GO TO 60
C
          IF(IEL.EQ.0) THEN
           CALL ASKI('Enter element number^',IEL)
           ITYPE = 4
          ENDIF
C
C-------- if this is the specified element, exit.
          IF(IEL .EQ. NEL) GO TO 60
          GO TO 50
        ENDIF
   55 CONTINUE
C      WRITE(*,5030) NMAX
C      WRITE(*,5900)
      IF(LOPEN) CLOSE(LU)
      ITYPE = 0
      RETURN
C
   60 N = I-1
      IF(LOPEN) CLOSE(LU)
      RETURN
C
   98 CONTINUE
C      WRITE(*,5050)
C      WRITE(*,5900)
      ITYPE = 0
      RETURN
C
   99 CONTINUE
      IF(LOPEN) CLOSE(LU)
C      WRITE(*,5100)
C      WRITE(*,5900)
      ITYPE = 0
      RETURN
C...............................................................
 1000 FORMAT(A)
 5030 FORMAT(/' Buffer array size exceeded'
     &       /' Maximum number of points: ', I4 )
 5050 FORMAT(/' File OPEN error.  Nonexistent file')
 5100 FORMAT(/' File READ error.  Unrecognizable file format')
 5900 FORMAT( ' *** LOAD NOT COMPLETED ***' )
      END ! AREAD
