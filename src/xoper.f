C***********************************************************************
C    Module:  xoper.f
C
C    Copyright (C) 2000 Mark Drela
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************
C
      SUBROUTINE OPER
      INCLUDE 'XFOIL.INC'



      CHARACTER*1 ANS
      CHARACTER*4 COMAND, COMOLD
      LOGICAL LRECALC, LCPX, LCONV
C
      CHARACTER*128 COMARG, ARGOLD, LINE
C
      PARAMETER (NPRX = 101)
      DIMENSION XPR(NPRX), YPR(NPRX)
C
      DIMENSION NBLP(NPX)
      DIMENSION IPPAI(NPX), NAPOLT(NPX)
C
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
C
C---- retain last-command info if OPER is exited and then re-entered
      SAVE COMOLD, ARGOLD
C
C---- logical units for  polar save file,  polar dump file
      LUPLR = 9
      LUPLX = 11
C
      COMAND = '****'
      COMARG = ' '
      LRECALC = .FALSE.
      LCPX = .FALSE.
      LPLOT = .FALSE.
      LEXITFLAG = .FALSE.
C
      IF(N.EQ.0) THEN
       WRITE(*,*)
       WRITE(*,*) '***  No airfoil available  ***'
       RETURN
      ENDIF
C
      IF(IPACT.NE.0) THEN
       WRITE(*,5000) IPACT
 5000  FORMAT(/'  Polar', I3,'  is active')
      ENDIF
C
ccc 500  CONTINUE
      COMOLD = COMAND
      ARGOLD = COMARG
C
C====================================================
C******************************************************
C	START OF DB COMMAND STRUCTURE
C
C	Main purpose is to calculate Cl/Cd either
C	through inviscid or viscous calculation.
C
c       LVISC = .NOT. LVISC
       LVISC = .TRUE.
C
c       IF(LVISC) THEN
C         IF(NINPUT.GE.1) THEN
C           REINF1 = RINPUT(1)
C         ELSE IF(REINF1 .EQ. 0.0) THEN
C           CALL ASKR('Enter Reynolds number^',REINF1)
C         ENDIF
C
C******************************************************
C	REINF1 = 100000
C	WRITE(*,*) REINF1
C	STOP

C*****CHANGE REYNOLDS NUMBER***************************
C       ITMAX = 75
C       CALL MRSHOW(.TRUE.,.TRUE.)
c       ENDIF
       LCONV = .FALSE.
       IF(.NOT.LRECALC) THEN
C------- set inviscid solution only if point is not being recalculated
C         IF(NINPUT.GE.1) THEN
C          ADEG = RINPUT(1)
C         ELSE
C          ADEG = ALFA/DTOR
C          CALL ASKR('Enter angle of attack (deg)^',ADEG)
C         ENDIF

C	 ADEG = 0.0	SET IN XFOIL.F NOW
         LALFA = .TRUE.
         ALFA = DTOR*ADEG
         QINF = 1.0

         CALL SPECAL

         IF(ABS(ALFA-AWAKE) .GT. 1.0E-5) LWAKE  = .FALSE.
         IF(ABS(ALFA-AVISC) .GT. 1.0E-5) LVCONV = .FALSE.
         IF(ABS(MINF-MVISC) .GT. 1.0E-5) LVCONV = .FALSE.
       ENDIF
C
       IF(LVISC) CALL VISCAL(ITMAX)
C       CALL CPX
       CALL FCPMIN

C
C       IF(LVISC .AND. LPACC .AND. LVCONV) THEN
C        CALL PLRADD(LUPLR,IPACT)
C        CALL PLXADD(LUPLX,IPACT)
C       ENDIF
C
C       IF(LVISC .AND. .NOT.LPACC .AND. .NOT.LVCONV) THEN
C        WRITE(*,*) 'Type "!" to continue iterating'
C       ENDIF


C      call cpcalc(N+NW,QVIS,QINF,MINF,CPV)
       call cdcalc
C      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF,XCMREF,YCMREF,
C     &            CL,CM,CDP, CL_ALF,CL_MSQ)

	RETURN
C
	END ! OPER

      SUBROUTINE FCPMIN
C------------------------------------------------
C     Finds minimum Cp on dist for cavitation work
C------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      XCPMNI = X(1)
      XCPMNV = X(1)
      CPMNI = CPI(1)
      CPMNV = CPV(1)
C
      DO I = 2, N + NW
        IF(CPI(I) .LT. CPMNI) THEN
         XCPMNI = X(I)
         CPMNI = CPI(I)
        ENDIF
        IF(CPV(I) .LT. CPMNV) THEN
         XCPMNV = X(I)
         CPMNV = CPV(I)
        ENDIF
      ENDDO
C

      IF(LVISC)THEN
        CPMN = CPMNV
      ELSE
        CPMN = CPMNI
C
        CPMNV = CPMNI
        XCPMNV = XCPMNI
      ENDIF
C
      RETURN
      END ! FCPMIN



      SUBROUTINE MRSHOW(LM,LR)
      INCLUDE 'XFOIL.INC'
      LOGICAL LM, LR
C
      IF(LM .OR. LR) WRITE(*,*)
C
      IF(LM) THEN
       IF(MATYP.EQ.1) WRITE(*,1100) MINF1
       IF(MATYP.EQ.2) WRITE(*,1100) MINF1, ' / sqrt(CL)'
       IF(MATYP.EQ.3) WRITE(*,1100) MINF1, ' / CL'
      ENDIF
C
      IF(LR) THEN
       IF(RETYP.EQ.1) WRITE(*,1200) INT(REINF1)
       IF(RETYP.EQ.2) WRITE(*,1200) INT(REINF1), ' / sqrt(CL)'
       IF(RETYP.EQ.3) WRITE(*,1200) INT(REINF1), ' / CL'
      ENDIF
C
      RETURN
C
 1100 FORMAT(1X,'M  =' , F10.4, A)
 1200 FORMAT(1X,'Re =' , I10  , A)
      END ! MRSHOW



      SUBROUTINE NAMMOD(NAME,KDEL,KMOD0)
      CHARACTER*(*) NAME
C-------------------------------------------
C     Requests new modified NAME with
C     version number in brackets, e.g.
C            NACA 0012  [5]
C
C     If bracketed index exists in NAME,
C        it is incremented by KDEL.
C     If no bracketed index exists, it
C        is added with initial value KMOD0,
C        unless KMOD0 is negative in which
C        case nothing is added.
C-------------------------------------------
      CHARACTER*48 NAMDEF
C
      CALL STRIP(NAME,NNAME)
      KBRACK1 = INDEX(NAME,'[')
      KBRACK2 = INDEX(NAME,']')
C
      NAMDEF = NAME(1:NNAME)
C
      IF(KBRACK1.NE.0 .AND.
     &   KBRACK2.NE.0 .AND. KBRACK2-KBRACK1.GT.1) THEN
C----- brackets exist... get number, (go get user's input on READ error)
       READ(NAME(KBRACK1+1:KBRACK2-1),*,ERR=40) KMOD
       KMOD = IABS(KMOD)
       KMODP = MOD( KMOD+KDEL , 100 )
       IF(KBRACK1.GE.2) THEN
        NAME = NAME(1:KBRACK1-1)
       ELSE
        NAME = ' '
       ENDIF
       CALL STRIP(NAME,NNAME)
      ELSEIF(KMOD0.GT.0) THEN
       KMODP = MOD( KMOD0 , 100 )
      ELSE
       KMODP = 0
      ENDIF
C
      IF    (KMODP.GE.10) THEN
       NAMDEF = NAME(1:NNAME) // ' [  ]'
       WRITE(NAMDEF(NNAME+3:NNAME+4),1020) KMODP
 1020  FORMAT(I2)
      ELSEIF(KMODP.GE. 1) THEN
       NAMDEF = NAME(1:NNAME) // ' [ ]'
       WRITE(NAMDEF(NNAME+3:NNAME+3),1025) KMODP
 1025  FORMAT(I1)
      ENDIF
C
 40   WRITE(*,1040) NAMDEF
 1040 FORMAT(/' Enter airfoil name or <return> for default:  ',A)
      READ(*,1000) NAME
 1000 FORMAT(A)
      IF(NAME .EQ. ' ') NAME = NAMDEF
C
      RETURN
      END ! NAMMOD



      SUBROUTINE BLDUMP(FNAME1)
      INCLUDE 'XFOIL.INC'
      CHARACTER*(*) FNAME1
C
      CHARACTER*80 FILDEF
C
 1000 FORMAT(A)
C
      IF(FNAME1(1:1).NE.' ') THEN
       FNAME = FNAME1
      ELSE
C----- no argument... get it somehow
       IF(NPREFIX.GT.0) THEN
C------ offer default using existing prefix
        FILDEF = PREFIX(1:NPREFIX) // '.bl'
        WRITE(*,1100) FILDEF
 1100   FORMAT(/' Enter filename:  ', A)
        READ(*,1000) FNAME
        CALL STRIP(FNAME,NFN)
        IF(NFN.EQ.0) FNAME = FILDEF
       ELSE
C------ nothing available... just ask for filename
        CALL ASKS('Enter filename^',FNAME)
       ENDIF
      ENDIF
C
      LU = 19
      OPEN(LU,FILE=FNAME,STATUS='UNKNOWN')
      REWIND(LU)
C
      WRITE(LU,1000)
     & '#    s        x        y     Ue/Vinf    Dstar     Theta      Cf'
C         1.23456  0.23451  0.23451  0.23451  0.012345  0.001234  0.004123
C
      CALL COMSET
      DO 10 I=1, N
        IS = 1
        IF(GAM(I) .LT. 0.0) IS = 2
C
        IF(LIPAN .AND. LVISC) THEN
          IF(IS.EQ.1) THEN
            IBL = IBLTE(IS) - I + 1
          ELSE
            IBL = IBLTE(IS) + I - N
          ENDIF
          DS = DSTR(IBL,IS)
          TH = THET(IBL,IS)
          CF =  TAU(IBL,IS)/(0.5*QINF**2)
        ELSE
          DS = 0.
          TH = 0.
          CF = 0.
        ENDIF
        UE = (GAM(I)/QINF)*(1.0-TKLAM) / (1.0 - TKLAM*(GAM(I)/QINF)**2)
C
        WRITE(LU,8500) S(I), X(I), Y(I), UE, DS, TH, CF
 8500   FORMAT(1X,4F9.5,3F10.6)
  10  CONTINUE
C
      IF(LWAKE) THEN
        IS = 2
        DO 20 I=N+1, N+NW
          IBL = IBLTE(IS) + I - N
          DS = DSTR(IBL,IS)
          TH = THET(IBL,IS)
          CF = 0.
          UI = UEDG(IBL,IS)
          UE = (UI/QINF)*(1.0-TKLAM) / (1.0 - TKLAM*(UI/QINF)**2)
C
          WRITE(LU,8500) S(I), X(I), Y(I), UE, DS, TH, CF
 20     CONTINUE
      ENDIF
C
      CLOSE(LU)
      RETURN
      END ! BLDUMP



      SUBROUTINE CPDUMP(FNAME1)
      INCLUDE 'XFOIL.INC'
      CHARACTER*(*) FNAME1
C
      CHARACTER*80 FILDEF
C
 1000 FORMAT(A)
C
      IF(FNAME1(1:1).NE.' ') THEN
       FNAME = FNAME1
      ELSE
C----- no argument... get it somehow
       IF(NPREFIX.GT.0) THEN
C------ offer default using existing prefix
        FILDEF = PREFIX(1:NPREFIX) // '.cp'
        WRITE(*,1100) FILDEF
 1100   FORMAT(/' Enter filename:  ', A)
        READ(*,1000) FNAME
        CALL STRIP(FNAME,NFN)
        IF(NFN.EQ.0) FNAME = FILDEF
       ELSE
C------ nothing available... just ask for filename
        CALL ASKS('Enter filename^',FNAME)
       ENDIF
      ENDIF
C
C
      LU = 19
      OPEN(LU,FILE=FNAME,STATUS='UNKNOWN')
      REWIND(LU)
C
      WRITE(LU,1000)
     & '#    x        Cp  '
C         0.23451  0.23451
C
      CALL COMSET
C
      BETA = SQRT(1.0 - MINF**2)
      BFAC = 0.5*MINF**2 / (1.0 + BETA)
C
      DO 10 I=1, N
        CPINC = 1.0 - (GAM(I)/QINF)**2
        DEN = BETA + BFAC*CPINC
        CPCOM = CPINC / DEN
C
        WRITE(LU,8500) X(I), CPCOM
 8500   FORMAT(1X,2F9.5)
  10  CONTINUE
C
      CLOSE(LU)
      RETURN
      END ! CPDUMP



      SUBROUTINE MHINGE
C----------------------------------------------------
C     Calculates the hinge moment of the flap about
C     (XOF,YOF) by integrating surface pressures.
C----------------------------------------------------
      INCLUDE 'XFOIL.INC'
C
      IF(.NOT.LFLAP) THEN
C
        CALL GETXYF(X,XP,Y,YP,S,N, TOPS,BOTS,XOF,YOF)
        LFLAP = .TRUE.
C
      ELSE
C
C------ find top and bottom y at hinge x location
        TOPS = XOF
        BOTS = S(N) - XOF
        CALL SINVRT(TOPS,XOF,X,XP,S,N)
        CALL SINVRT(BOTS,XOF,X,XP,S,N)
C
      ENDIF
C
      TOPX = SEVAL(TOPS,X,XP,S,N)
      TOPY = SEVAL(TOPS,Y,YP,S,N)
      BOTX = SEVAL(BOTS,X,XP,S,N)
      BOTY = SEVAL(BOTS,Y,YP,S,N)
C
C
      HMOM = 0.
      HFX  = 0.
      HFY  = 0.
C
C---- integrate pressures on top and bottom sides of flap
      DO 20 I=2, N
        IF(S(I-1).GE.TOPS .AND. S(I).LE.BOTS) GO TO 20
C
         DX = X(I) - X(I-1)
         DY = Y(I) - Y(I-1)
         XMID = 0.5*(X(I)+X(I-1)) - XOF
         YMID = 0.5*(Y(I)+Y(I-1)) - YOF
         IF(LVISC) THEN
          PMID = 0.5*(CPV(I) + CPV(I-1))
         ELSE
          PMID = 0.5*(CPI(I) + CPI(I-1))
         ENDIF
         HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
         HFX  = HFX  - PMID* DY
         HFY  = HFY  + PMID* DX
   20 CONTINUE
C
C---- find S(I)..S(I-1) interval containing s=TOPS
      DO I=2, N
        IF(S(I).GT.TOPS) GO TO 31
      ENDDO
C
   31 CONTINUE
C---- add on top surface chunk TOPS..S(I-1),  missed in the DO 20 loop.
      DX = TOPX - X(I-1)
      DY = TOPY - Y(I-1)
      XMID = 0.5*(TOPX+X(I-1)) - XOF
      YMID = 0.5*(TOPY+Y(I-1)) - YOF
      IF(S(I) .NE. S(I-1)) THEN
       FRAC = (TOPS-S(I-1))/(S(I)-S(I-1))
      ELSE
       FRAC = 0.
      ENDIF
      IF(LVISC) THEN
       TOPP = CPV(I)*FRAC + CPV(I-1)*(1.0-FRAC)
       PMID = 0.5*(TOPP+CPV(I-1))
      ELSE
       TOPP = CPI(I)*FRAC + CPI(I-1)*(1.0-FRAC)
       PMID = 0.5*(TOPP+CPI(I-1))
      ENDIF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
C
C---- add on inside flap surface contribution from hinge to top surface
      DX = XOF - TOPX
      DY = YOF - TOPY
      XMID = 0.5*(TOPX+XOF) - XOF
      YMID = 0.5*(TOPY+YOF) - YOF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
C
C---- find S(I)..S(I-1) interval containing s=BOTS
      DO I=N, 2, -1
        IF(S(I-1).LT.BOTS) GO TO 41
      ENDDO
C
   41 CONTINUE
C---- add on bottom surface chunk BOTS..S(I),  missed in the DO 20 loop.
      DX = X(I) - BOTX
      DY = Y(I) - BOTY
      XMID = 0.5*(BOTX+X(I)) - XOF
      YMID = 0.5*(BOTY+Y(I)) - YOF
      IF(S(I) .NE. S(I-1)) THEN
       FRAC = (BOTS-S(I-1))/(S(I)-S(I-1))
      ELSE
       FRAC = 0.
      ENDIF
      IF(LVISC) THEN
       BOTP = CPV(I)*FRAC + CPV(I-1)*(1.0-FRAC)
       PMID = 0.5*(BOTP+CPV(I))
      ELSE
       BOTP = CPI(I)*FRAC + CPI(I-1)*(1.0-FRAC)
       PMID = 0.5*(BOTP+CPI(I))
      ENDIF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
C
C---- add on inside flap surface contribution from hinge to bottom surface
      DX = BOTX - XOF
      DY = BOTY - YOF
      XMID = 0.5*(BOTX+XOF) - XOF
      YMID = 0.5*(BOTY+YOF) - YOF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
C
C---- add on TE base thickness contribution
      DX = X(1) - X(N)
      DY = Y(1) - Y(N)
      XMID = 0.5*(X(1)+X(N)) - XOF
      YMID = 0.5*(Y(1)+Y(N)) - YOF
      IF(LVISC) THEN
       PMID = 0.5*(CPV(1)+CPV(N))
      ELSE
       PMID = 0.5*(CPI(1)+CPI(N))
      ENDIF
      HMOM = HMOM + PMID*(XMID*DX + YMID*DY)
      HFX  = HFX  - PMID* DY
      HFY  = HFY  + PMID* DX
C
      RETURN
      END ! MHINGE


      SUBROUTINE VPAR
C---------------------------------------------
C     Viscous parameter change menu routine.
C---------------------------------------------
      INCLUDE 'XFOIL.INC'
      CHARACTER*4 COMAND
      CHARACTER*128 COMARG
C
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
C
C
      TURB = 100.0 * EXP( -(ACRIT + 8.43)/2.4 )
      WRITE(*,1200) XSTRIP(1), XSTRIP(2), ACRIT, TURB, VACCEL
C
  500 CONTINUE
      CALL ASKC('..VPAR^',COMAND,COMARG)
C
      DO I=1, 20
        IINPUT(I) = 0
        RINPUT(I) = 0.0
      ENDDO
      NINPUT = 20
      CALL GETINT(COMARG,IINPUT,NINPUT,ERROR)
      NINPUT = 20
      CALL GETFLT(COMARG,RINPUT,NINPUT,ERROR)
C
      IF(COMAND.EQ.'    ') RETURN
      IF(COMAND.EQ.'?   ') GO TO 5
      IF(COMAND.EQ.'SHOW') GO TO 10
      IF(COMAND.EQ.'XTR ') GO TO 40
      IF(COMAND.EQ.'N   ') GO TO 50
      IF(COMAND.EQ.'VACC') GO TO 70
      IF(COMAND.EQ.'INIT') GO TO 80
C
      WRITE(*,1000) COMAND
      GO TO 500
C
    5 WRITE(*,1050)
      GO TO 500
C
   10 TURB = 100.0 * EXP( -(ACRIT + 8.43)/2.4 )
      WRITE(*,1200) XSTRIP(1), XSTRIP(2), ACRIT, TURB, VACCEL
      GO TO 500
C
   40 IF(NINPUT.GE.2) THEN
       XSTRIP(1) = RINPUT(1)
       XSTRIP(2) = RINPUT(2)
      ELSE
       CALL ASKR('Enter top    side Xtrip/c^',XSTRIP(1))
       CALL ASKR('Enter bottom side Xtrip/c^',XSTRIP(2))
      ENDIF
      LVCONV = .FALSE.
      GO TO 500
C
   50 IF(NINPUT.GE.1) THEN
       ACRIT = RINPUT(1)
      ELSE
       CALL ASKR('Enter critical amplification ratio^',ACRIT)
      ENDIF
      LVCONV = .FALSE.
      GO TO 500
C
   70 IF(NINPUT.GE.1) THEN
       VACCEL = RINPUT(1)
      ELSE
       CALL ASKR('Enter viscous acceleration parameter^',VACCEL)
      ENDIF
      GO TO 500
C
   80 LBLINI = .NOT.LBLINI
      IF(.NOT.LBLINI) WRITE(*,*) 'BLs will be initialized on next point'
      IF(     LBLINI) WRITE(*,*) 'BLs are assumed to be initialized'
      IF(.NOT.LBLINI) LIPAN = .FALSE.
      GO TO 500
C
C...................................................................
C
 1000 FORMAT(1X,A4,' command not recognized.  Type a "?" for list')
 1050 FORMAT(
     & /'   <cr>    Return to OPER menu'
     & /'   SHOW    Display viscous parameters'
     & /'   XTR  rr Change trip positions Xtr/c'
     & /'   N    r  Change critical amplification exponent Ncrit'
     & /'   VACC r  Change Newton solution acceleration parameter'
     & /'   INIT    BL initialization flag toggle')
 1200 FORMAT(/' Xtr/c     =', F8.4, '    top    side'
     &       /' Xtr/c     =', F8.4, '    bottom side'
     &       /' Ncrit     =', F8.2, '   (', F6.3, ' % turb. level )'
     &       /' Vacc      =', F8.4  )
      END ! VPAR




      SUBROUTINE SPECAL
C-----------------------------------
C     Converges to specified alpha.
C-----------------------------------
      INCLUDE 'XFOIL.INC'
      REAL MINF_CLM, MSQ_CLM
C
C---- calculate surface vorticity distributions for alpha = 0, 90 degrees
      IF(.NOT.LGAMU .OR. .NOT.LQAIJ) CALL GGCALC
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
C---- superimpose suitably weighted  alpha = 0, 90  distributions
      DO 50 I=1, N
        GAM(I)   =  COSA*GAMU(I,1) + SINA*GAMU(I,2)
        GAM_A(I) = -SINA*GAMU(I,1) + COSA*GAMU(I,2)
   50 CONTINUE
      PSIO = COSA*GAMU(N+1,1) + SINA*GAMU(N+1,2)
C
      CALL TECALC
      CALL QISET
C
C---- set initial guess for the Newton variable CLM
      CLM = 1.0
C
C---- set corresponding  M(CLM), Re(CLM)
      CALL MRCL(CLM,MINF_CLM,REINF_CLM)
      CALL COMSET
C
C---- set corresponding CL(M)
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
C
C---- iterate on CLM
      DO 100 ITCL=1, 20
C
        MSQ_CLM = 2.0*MINF*MINF_CLM
        DCLM = (CL - CLM)/(1.0 - CL_MSQ*MSQ_CLM)
C
        CLM1 = CLM
        RLX = 1.0
C
C------ under-relaxation loop to avoid driving M(CL) above 1
        DO 90 IRLX=1, 12
C
          CLM = CLM1 + RLX*DCLM
C
C-------- set new freestream Mach M(CLM)
          CALL MRCL(CLM,MINF_CLM,REINF_CLM)
C
C-------- if Mach is OK, go do next Newton iteration
          IF(MATYP.EQ.1 .OR. MINF.EQ.0.0 .OR. MINF_CLM.NE.0.0) GO TO 91
C
          RLX = 0.5*RLX
   90   CONTINUE
   91   CONTINUE
C
C------ set new CL(M)
        CALL COMSET
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &              CL,CM,CDP,CL_ALF,CL_MSQ)
C
        IF(ABS(DCLM).LE.1.0E-6) GO TO 110
C
  100 CONTINUE
      WRITE(*,*) 'SPECAL:  Minf convergence failed'
  110 CONTINUE
C
C---- set final Mach, CL, Cp distributions, and hinge moment
      CALL MRCL(CL,MINF_CL,REINF_CL)
      CALL COMSET
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
      CALL CPCALC(N,QINV,QINF,MINF,CPI)
      IF(LVISC) THEN
       CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
       CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
      ELSE
       CALL CPCALC(N,QINV,QINF,MINF,CPI)
      ENDIF
      IF(LFLAP) CALL MHINGE
C
      RETURN
      END ! SPECAL


      SUBROUTINE SPECCL
C-----------------------------------------
C     Converges to specified inviscid CL.
C-----------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- calculate surface vorticity distributions for alpha = 0, 90 degrees
      IF(.NOT.LGAMU .OR. .NOT.LQAIJ) CALL GGCALC
C
C---- set freestream Mach from specified CL -- Mach will be held fixed
      CALL MRCL(CLSPEC,MINF_CL,REINF_CL)
      CALL COMSET
C
C---- current alpha is the initial guess for Newton variable ALFA
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
      DO 10 I=1, N
        GAM(I)   =  COSA*GAMU(I,1) + SINA*GAMU(I,2)
        GAM_A(I) = -SINA*GAMU(I,1) + COSA*GAMU(I,2)
   10 CONTINUE
      PSIO = COSA*GAMU(N+1,1) + SINA*GAMU(N+1,2)
C
C---- get corresponding CL, CL_alpha, CL_Mach
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
C
C---- Newton loop for alpha to get specified inviscid CL
      DO 100 ITAL=1, 20
C
        DALFA = (CLSPEC - CL) / CL_ALF
        RLX = 1.0
C
        ALFA = ALFA + RLX*DALFA
C
C------ set new surface speed distribution
        COSA = COS(ALFA)
        SINA = SIN(ALFA)
        DO 40 I=1, N
          GAM(I)   =  COSA*GAMU(I,1) + SINA*GAMU(I,2)
          GAM_A(I) = -SINA*GAMU(I,1) + COSA*GAMU(I,2)
   40   CONTINUE
        PSIO = COSA*GAMU(N+1,1) + SINA*GAMU(N+1,2)
C
C------ set new CL(alpha)
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &              CL,CM,CDP,CL_ALF,CL_MSQ)
C
        IF(ABS(DALFA).LE.1.0E-6) GO TO 110
  100 CONTINUE
      WRITE(*,*) 'SPECCL:  CL convergence failed'
  110 CONTINUE
C
C---- set final surface speed and Cp distributions
      CALL TECALC
      CALL QISET
      IF(LVISC) THEN
       CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
       CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
      ELSE
       CALL CPCALC(N,QINV,QINF,MINF,CPI)
      ENDIF
      IF(LFLAP) CALL MHINGE
C
      RETURN
      END ! SPECCL


      SUBROUTINE VISCAL(NITER1)
C----------------------------------------
C     Converges viscous operating point
C----------------------------------------
      INCLUDE 'XFOIL.INC'
C
C---- convergence tolerance
      DATA EPS1 / 1.0E-10 /
C
      NITER = NITER1
C
C---- calculate wake trajectory from current inviscid solution if necessary
      IF(.NOT.LWAKE) THEN
       CALL XYWAKE
      ENDIF
C
C---- set velocities on wake from airfoil vorticity for alpha=0, 90
      CALL QWCALC
C
C---- set velocities on airfoil and wake for initial alpha
      CALL QISET
C
      IF(.NOT.LIPAN) THEN
C
       IF(LBLINI) CALL GAMQV
C
C----- locate stagnation point arc length position and panel index
       CALL STFIND
C
C----- set  BL position -> panel position  pointers
       CALL IBLPAN
C
C----- calculate surface arc length array for current stagnation point location
       CALL XICALC
C
C----- set  BL position -> system line  pointers
       CALL IBLSYS
C
      ENDIF
C
C---- set inviscid BL edge velocity UINV from QINV
      CALL UICALC
C
      IF(.NOT.LBLINI) THEN
C
C----- set initial Ue from inviscid Ue
       DO IBL=1, NBL(1)
         UEDG(IBL,1) = UINV(IBL,1)
       ENDDO
C
       DO IBL=1, NBL(2)
         UEDG(IBL,2) = UINV(IBL,2)
       ENDDO
C
      ENDIF
C
      IF(LVCONV) THEN
C----- set correct CL if converged point exists
       CALL QVFUE
       IF(LVISC) THEN
        CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
        CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
       ELSE
        CALL CPCALC(N,QINV,QINF,MINF,CPI)
       ENDIF
       CALL GAMQV
       CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &             CL,CM,CDP, CL_ALF,CL_MSQ)
       CALL CDCALC
      ENDIF
C
C---- set up source influence matrix if it doesn't exist
      IF(.NOT.LWDIJ .OR. .NOT.LADIJ) CALL QDCALC
C
C---- Newton iteration for entire BL solution
      IF(NITER.EQ.0) CALL ASKI('Enter number of iterations^',NITER)
C      WRITE(*,*)
C      WRITE(*,*) 'Solving BL system ...'
      DO 1000 ITER=1, NITER
C
C------ fill Newton system for BL variables
c        WRITE (*,*) 'Calling SETBL...'
	CALL SETBL

c------ DB 040306
c	WRITE (*,*) 'Exit flag = ', LEXITFLAG
	IF (LEXITFLAG) THEN
		WRITE(*,*) 'LEXITFLAG TRUE, GOING TO 90...'
		GOTO 90
c		WRITE(*,*) 'Exiting on account of test...'
c		STOP
	ENDIF

C
C------ solve Newton system with custom solver
c	WRITE(*,*) 'CALLING BLSOLV...'
        CALL BLSOLV
C
C------ update BL variables
        CALL UPDATE
C
        IF(LALFA) THEN
C------- set new freestream Mach, Re from new CL
         CALL MRCL(CL,MINF_CL,REINF_CL)
         CALL COMSET
        ELSE
C------- set new inviscid speeds QINV and UINV for new alpha
         CALL QISET
         CALL UICALC
        ENDIF
C
C------ calculate edge velocities QVIS(.) from UEDG(..)
        CALL QVFUE
C
C------ set GAM distribution from QVIS
        CALL GAMQV
C
C------ relocate stagnation point
        CALL STMOVE
C
C------ set updated CL,CD
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &              CL,CM,CDP,CL_ALF,CL_MSQ)
        CALL CDCALC
C
C------ display changes and test for convergence
c        IF(RLX.LT.1.0)
c     &   WRITE(*,2000) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL,RLX
c        IF(RLX.EQ.1.0)
c     &   WRITE(*,2010) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL
c        CDP = CD - CDF
c         WRITE(*,2020) ALFA/DTOR, CL, CM, CD, CDF, CDP
C
        IF(RMSBL .LT. EPS1) THEN
         LVCONV = .TRUE.
         AVISC = ALFA
         MVISC = MINF
         GO TO 90
        ENDIF
C
 1000 CONTINUE
C      WRITE(*,*) 'VISCAL:  Convergence failed'
C
   90 CONTINUE
      CALL CPCALC(N+NW,QINV,QINF,MINF,CPI)
      CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV)
      IF(LFLAP) CALL MHINGE
      RETURN
C....................................................................
 2000   FORMAT
     &   (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3,
     &     '   RLX:',F6.3)
 2010   FORMAT
     &   (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3)
 2020   FORMAT
     &   ( 1X,3X,'   a =', F7.3,'      CL =',F8.4  /
     &     1X,3X,'  Cm =', F8.4, '     CD =',F9.5,
     &           '   =>   CDf =',F9.5,'    CDp =',F9.5)
      END ! VISCAL


      subroutine dcpout
      include 'XFOIL.INC'
c
c     Computes and writes upper and lower-surface
c     Cp values at two specified x locations
c
c
      x1 = 0.05
      x2 = 0.15
c
      lu = 60
      open(lu,file='dcp.out',status='old',access='append',err=10)
      go to 20
c
 10   continue
      open(lu,file='dcp.out',status='new')
      write(lu,*) '#  ', name
      write(lu,*) '# alpha   CL       ',
     &            ' Cpl05     Cpu05     dCp05    ',
     &            ' Cpl15     Cpu15     dCp15    '
 20   continue
c
      call spline(cpv,w1,s,n)
c
      su1 = sle + x1*(s(1)-sle)
      sl1 = sle + x1*(s(n)-sle)
      su2 = sle + x2*(s(1)-sle)
      sl2 = sle + x2*(s(n)-sle)
c
      call sinvrt(sl1,x1,x,xp,s,n)
      call sinvrt(su1,x1,x,xp,s,n)
      call sinvrt(sl2,x2,x,xp,s,n)
      call sinvrt(su2,x2,x,xp,s,n)
c
      cpl1 = seval(sl1,cpv,w1,s,n)
      cpu1 = seval(su1,cpv,w1,s,n)
      cpl2 = seval(sl2,cpv,w1,s,n)
      cpu2 = seval(su2,cpv,w1,s,n)
c
      write(lu,1200) alfa/dtor, cl,
     &               cpl1, cpu1, cpl1-cpu1,
     &               cpl2, cpu2, cpl2-cpu2

 1200 format(1x, f7.3, f9.4, 8f10.5)
c
      close(lu)
c
      return
      end
