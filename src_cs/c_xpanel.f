C***********************************************************************
C    Module:  xpanel.f
C
C    Copyright (C) 2000 Mark Drela
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


      SUBROUTINE APCALC
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
C---- set angles of airfoil panels
      DO 10 I=1, N-1
        SX = X(I+1) - X(I)
        SY = Y(I+1) - Y(I)
        IF((SX.ceq.0.0) .AND. (SY.ceq.0.0)) THEN
          APANEL(I) = ATAN2( -NY(I) , -NX(I) )
        ELSE
          APANEL(I) = ATAN2( SX , -SY )
        ENDIF
   10 CONTINUE
C
C---- TE panel
      I = N
      IP = 1
      IF(SHARP) THEN
       APANEL(I) = PI
      ELSE
       SX = X(IP) - X(I)
       SY = Y(IP) - Y(I)
       APANEL(I) = ATAN2( -SX , SY ) + PI
      ENDIF
C
      RETURN
      END


      SUBROUTINE NCALC(X,Y,S,N,XN,YN)
C---------------------------------------
C     Calculates normal unit vector
C     components at airfoil panel nodes
C---------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      DIMENSION X(N), Y(N), S(N), XN(N), YN(N)
C
      IF(N.LE.1) RETURN
C
      CALL SEGSPL(X,XN,S,N)
      CALL SEGSPL(Y,YN,S,N)
      DO 10 I=1, N
        SX =  YN(I)
        SY = -XN(I)
        SMOD = SQRT(SX*SX + SY*SY)
        XN(I) = SX/SMOD
        YN(I) = SY/SMOD
   10 CONTINUE
C
C---- average normal vectors at corner points
      DO 20 I=1, N-1
        IF(S(I) .ceq. S(I+1)) THEN
          SX = 0.5*(XN(I) + XN(I+1))
          SY = 0.5*(YN(I) + YN(I+1))
          SMOD = SQRT(SX*SX + SY*SY)
          XN(I)   = SX/SMOD
          YN(I)   = SY/SMOD
          XN(I+1) = SX/SMOD
          YN(I+1) = SY/SMOD
        ENDIF
 20   CONTINUE
C
      RETURN
      END


      SUBROUTINE PSILIN(I,XI,YI,NXI,NYI,PSI,PSI_NI,GEOLIN,SIGLIN)
C-----------------------------------------------------------------------
C     Calculates current streamfunction Psi at panel node or wake node
C     I due to freestream and all bound vorticity Gam on the airfoil.
C     Sensitivities of Psi with respect to alpha (Z_ALFA) and inverse
C     Qspec DOFs (Z_QDOF0,Z_QDOF1) which influence Gam in inverse cases.
C     Also calculates the sensitivity vector dPsi/dGam (DZDG).
C
C     If SIGLIN=True, then Psi includes the effects of the viscous
C     source distribution Sig and the sensitivity vector dPsi/dSig
C     (DZDM) is calculated.
C
C     If GEOLIN=True, then the geometric sensitivity vector dPsi/dn
C     is calculated, where n is the normal motion of the jth node.
C
C          Airfoil:  1   < I < N
C          Wake:     N+1 < I < N+NW
C-----------------------------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
      complex NXO, NYO, NXP, NYP, NXI, NYI
      LOGICAL GEOLIN,SIGLIN
C
C---- distance tolerance for determining if two points are the same
      SEPS = (S(N)-S(1)) * 1.0E-5

      IO = I
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
      DO 3 JO=1, N
        DZDG(JO) = 0.0
        DZDN(JO) = 0.0
        DQDG(JO) = 0.0
    3 CONTINUE
C
      DO 4 JO=1, N
        DZDM(JO) = 0.0
        DQDM(JO) = 0.0
    4 CONTINUE
C
      Z_QINF = 0.
      Z_ALFA = 0.
      Z_QDOF0 = 0.
      Z_QDOF1 = 0.
      Z_QDOF2 = 0.
      Z_QDOF3 = 0.
C
      PSI    = 0.
      PSI_NI = 0.
C


      QTAN1 = 0.
      QTAN2 = 0.
      QTANM = 0.
C
      IF(SHARP) THEN
       SCS = 1.0
       SDS = 0.0
      ELSE
       SCS = ANTE/DSTE
       SDS = ASTE/DSTE
      ENDIF
C
      DO 10 JO=1, N
        JP = JO+1
C
        JM = JO-1
        JQ = JP+1
C
        IF(JO.ceq.1) THEN
         JM = JO
        ELSE IF(JO.ceq.N-1) THEN
         JQ = JP
        ELSE IF(JO.ceq.N) THEN
         JP = 1
         IF((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2 .LT. SEPS**2) GO TO 12
        ENDIF
C
        DSO = SQRT((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2)
C
C------ skip null panel
        IF(DSO .ceq. 0.0) GO TO 10
C
        DSIO = 1.0 / DSO
C
        APAN = APANEL(JO)
C
        RX1 = XI - X(JO)
        RY1 = YI - Y(JO)
        RX2 = XI - X(JP)
        RY2 = YI - Y(JP)
C
        SX = (X(JP) - X(JO)) * DSIO
        SY = (Y(JP) - Y(JO)) * DSIO
C
        X1 = SX*RX1 + SY*RY1
        X2 = SX*RX2 + SY*RY2
        YY = SX*RY1 - SY*RX1
C
        RS1 = RX1*RX1 + RY1*RY1
        RS2 = RX2*RX2 + RY2*RY2


C------ set reflection flag SGN to avoid branch problems with arctan
        IF(IO.GE.1 .AND. IO.LE.N) THEN
C------- no problem on airfoil surface
         SGN = 1.0
        ELSE
C------- make sure arctan falls between  -/+  Pi/2
         SGN = SIGN(1.0,YY)
        ENDIF
C
C------ set log(r^2) and arctan(x/y), correcting for reflection if any
        IF((IO.cne.JO) .AND. RS1.GT.0.0) THEN
         G1 = LOG(RS1)
         T1 = ATAN2(SGN*X1,SGN*YY) + (0.5 - 0.5*SGN)*PI
        ELSE
         G1 = 0.0
         T1 = 0.0
        ENDIF
C
        IF((IO.cne.JP) .AND. RS2.GT.0.0) THEN
         G2 = LOG(RS2)
         T2 = ATAN2(SGN*X2,SGN*YY) + (0.5 - 0.5*SGN)*PI
        ELSE
         G2 = 0.0
         T2 = 0.0
        ENDIF
c        if (aimag(G1) .ne.0)  write(*,*) 'G1',G1
c        if (aimag(T1) .ne.0)  write(*,*) 'T1',T1
c        if (aimag(T2) .ne.0)  write(*,*) 'T2',T2

        X1I = SX*NXI + SY*NYI
c        if (aimag(X1I) .ne. 0) then
c           write(*,*) 'SX',SX
c           write(*,*) 'NXI',NXI
c           write(*,*) 'SY',SY
c           write(*,*) 'NYI',NYI
c        endif

        X2I = SX*NXI + SY*NYI
        YYI = SX*NYI - SY*NXI
C
        IF(GEOLIN) THEN
         NXO = NX(JO)
         NYO = NY(JO)
         NXP = NX(JP)
         NYP = NY(JP)
C
         X1O =-((RX1-X1*SX)*NXO + (RY1-X1*SY)*NYO)*DSIO-(SX*NXO+SY*NYO)
         X1P = ((RX1-X1*SX)*NXP + (RY1-X1*SY)*NYP)*DSIO
         X2O =-((RX2-X2*SX)*NXO + (RY2-X2*SY)*NYO)*DSIO
         X2P = ((RX2-X2*SX)*NXP + (RY2-X2*SY)*NYP)*DSIO-(SX*NXP+SY*NYP)
         YYO = ((RX1+X1*SY)*NYO - (RY1-X1*SX)*NXO)*DSIO-(SX*NYO-SY*NXO)
         YYP =-((RX1-X1*SY)*NYP - (RY1+X1*SX)*NXP)*DSIO
        ENDIF
C
        IF(JO.ceq.N) GO TO 11
C
        IF(SIGLIN) THEN
C
C------- set up midpoint quantities
         X0 = 0.5*(X1+X2)
         RS0 = X0*X0 + YY*YY
         G0 = LOG(RS0)
         T0 = ATAN2(SGN*X0,SGN*YY) + (0.5 - 0.5*SGN)*PI
C
C------- calculate source contribution to Psi  for  1-0  half-panel
         DXINV = 1.0/(X1-X0)
         PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0)
         PDIF = ((X1+X0)*PSUM + RS1*(T1-APAN) - RS0*(T0-APAN)
     &        + (X0-X1)*YY) * DXINV
C
         PSX1 =  -(T1-APAN)
         PSX0 =    T0-APAN
         PSYY =  0.5*(G1-G0)
C
         PDX1 = ((X1+X0)*PSX1 + PSUM + 2.0*X1*(T1-APAN) - PDIF) * DXINV
         PDX0 = ((X1+X0)*PSX0 + PSUM - 2.0*X0*(T0-APAN) + PDIF) * DXINV
         PDYY = ((X1+X0)*PSYY + 2.0*(X0-X1 + YY*(T1-T0))      ) * DXINV
C
         DSM = SQRT((X(JP)-X(JM))**2 + (Y(JP)-Y(JM))**2)
         DSIM = 1.0/DSM
C
CCC      SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC      SIG1 = (SIG(JP) - SIG(JM))*DSIM
CCC      SSUM = SIG0 + SIG1
CCC      SDIF = SIG0 - SIG1
C
         SSUM = (SIG(JP) - SIG(JO))*DSIO + (SIG(JP) - SIG(JM))*DSIM
         SDIF = (SIG(JP) - SIG(JO))*DSIO - (SIG(JP) - SIG(JM))*DSIM
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JM) = DZDM(JM) + QOPI*(-PSUM*DSIM + PDIF*DSIM)
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*DSIO - PDIF*DSIO)
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*(DSIO+DSIM)
     &                                          + PDIF*(DSIO-DSIM))
C
C------- dPsi/dni
         PSNI = PSX1*X1I + PSX0*(X1I+X2I)*0.5 + PSYY*YYI
         PDNI = PDX1*X1I + PDX0*(X1I+X2I)*0.5 + PDYY*YYI

         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)

c         if (aimag(PSI_NI) .ne. 0) write(*,*) 'position 1',PSI_NI
C
         QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JM) = DQDM(JM) + QOPI*(-PSNI*DSIM + PDNI*DSIM)
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*DSIO - PDNI*DSIO)
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*(DSIO+DSIM)
     &                                          + PDNI*(DSIO-DSIM))
C
C
C------- calculate source contribution to Psi  for  0-2  half-panel
         DXINV = 1.0/(X0-X2)
         PSUM = X2*(T2-APAN) - X0*(T0-APAN) + 0.5*YY*(G0-G2)
         PDIF = ((X0+X2)*PSUM + RS0*(T0-APAN) - RS2*(T2-APAN)
     &        + (X2-X0)*YY) * DXINV
C
         PSX0 =  -(T0-APAN)
         PSX2 =    T2-APAN
         PSYY =  0.5*(G0-G2)
C
         PDX0 = ((X0+X2)*PSX0 + PSUM + 2.0*X0*(T0-APAN) - PDIF) * DXINV
         PDX2 = ((X0+X2)*PSX2 + PSUM - 2.0*X2*(T2-APAN) + PDIF) * DXINV
         PDYY = ((X0+X2)*PSYY + 2.0*(X2-X0 + YY*(T0-T2))      ) * DXINV
C
         DSP = SQRT((X(JQ)-X(JO))**2 + (Y(JQ)-Y(JO))**2)
         DSIP = 1.0/DSP
C
CCC         SIG2 = (SIG(JQ) - SIG(JO))*DSIP
CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC         SSUM = SIG2 + SIG0
CCC         SDIF = SIG2 - SIG0
C
         SSUM = (SIG(JQ) - SIG(JO))*DSIP + (SIG(JP) - SIG(JO))*DSIO
         SDIF = (SIG(JQ) - SIG(JO))*DSIP - (SIG(JP) - SIG(JO))*DSIO
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*(DSIP+DSIO)
     &                                          - PDIF*(DSIP-DSIO))
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*DSIO - PDIF*DSIO)
         DZDM(JQ) = DZDM(JQ) + QOPI*( PSUM*DSIP + PDIF*DSIP)
C
C------- dPsi/dni
         PSNI = PSX0*(X1I+X2I)*0.5 + PSX2*X2I + PSYY*YYI
         PDNI = PDX0*(X1I+X2I)*0.5 + PDX2*X2I + PDYY*YYI

         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
c         if (aimag(PSI_NI) .ne. 0) write(*,*) 'position 2',PSI_NI
C
         QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*(DSIP+DSIO)
     &                                          - PDNI*(DSIP-DSIO))
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*DSIO - PDNI*DSIO)
         DQDM(JQ) = DQDM(JQ) + QOPI*( PSNI*DSIP + PDNI*DSIP)
C
        ENDIF
C
C------ calculate vortex panel contribution to Psi
        DXINV = 1.0/(X1-X2)
        PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
        PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2-RS1*G1 + X1*X1-X2*X2))*DXINV
C
        PSX1 = 0.5*G1
        PSX2 = -.5*G2
        PSYY = T1-T2
C
        PDX1 = ((X1+X2)*PSX1 + PSIS - X1*G1 - PSID)*DXINV
        PDX2 = ((X1+X2)*PSX2 + PSIS + X2*G2 + PSID)*DXINV
        PDYY = ((X1+X2)*PSYY - YY*(G1-G2)         )*DXINV
C
        GSUM1 = GAMU(JP,1) + GAMU(JO,1)
        GSUM2 = GAMU(JP,2) + GAMU(JO,2)
        GDIF1 = GAMU(JP,1) - GAMU(JO,1)
        GDIF2 = GAMU(JP,2) - GAMU(JO,2)
C
        GSUM = GAM(JP) + GAM(JO)
        GDIF = GAM(JP) - GAM(JO)
C
        PSI = PSI + QOPI*(PSIS*GSUM + PSID*GDIF)
C
C------ dPsi/dGam
        DZDG(JO) = DZDG(JO) + QOPI*(PSIS-PSID)
        DZDG(JP) = DZDG(JP) + QOPI*(PSIS+PSID)
C
C------ dPsi/dni
        PSNI = PSX1*X1I + PSX2*X2I + PSYY*YYI
        PDNI = PDX1*X1I + PDX2*X2I + PDYY*YYI
        PSI_NI = PSI_NI + QOPI*(GSUM*PSNI + GDIF*PDNI)
c        if (aimag(PSI_NI) .ne. 0) then
c           write(*,*) 'position 3',PSI_NI
c           write(*,*) 'QOPI',QOPI
c           write(*,*) 'GSUM', GSUM
c           write(*,*) 'PSNI',PSNI
c           write(*,*) 'GDIF',GDIF
c           write(*,*) 'PDNI',PDNI
c           write(*,*) 'PSX1',PSX1
c           write(*,*) 'X1I',X1I
c           write(*,*) 'PSX2',PSX2
c           write(*,*) 'X2I',X2I
c           write(*,*) 'PSYY',PSYY
c           write(*,*) 'YYI',YYI
c        endif

C
        QTAN1 = QTAN1 + QOPI*(GSUM1*PSNI + GDIF1*PDNI)
        QTAN2 = QTAN2 + QOPI*(GSUM2*PSNI + GDIF2*PDNI)
C
        DQDG(JO) = DQDG(JO) + QOPI*(PSNI - PDNI)
        DQDG(JP) = DQDG(JP) + QOPI*(PSNI + PDNI)
C
        IF(GEOLIN) THEN
C
C------- dPsi/dn
         DZDN(JO) = DZDN(JO)+ QOPI*GSUM*(PSX1*X1O + PSX2*X2O + PSYY*YYO)
     &                      + QOPI*GDIF*(PDX1*X1O + PDX2*X2O + PDYY*YYO)
         DZDN(JP) = DZDN(JP)+ QOPI*GSUM*(PSX1*X1P + PSX2*X2P + PSYY*YYP)
     &                      + QOPI*GDIF*(PDX1*X1P + PDX2*X2P + PDYY*YYP)
C------- dPsi/dP
         Z_QDOF0 = Z_QDOF0
     &           + QOPI*((PSIS-PSID)*QF0(JO) + (PSIS+PSID)*QF0(JP))
         Z_QDOF1 = Z_QDOF1
     &           + QOPI*((PSIS-PSID)*QF1(JO) + (PSIS+PSID)*QF1(JP))
         Z_QDOF2 = Z_QDOF2
     &           + QOPI*((PSIS-PSID)*QF2(JO) + (PSIS+PSID)*QF2(JP))
         Z_QDOF3 = Z_QDOF3
     &           + QOPI*((PSIS-PSID)*QF3(JO) + (PSIS+PSID)*QF3(JP))
        ENDIF
C
C
   10 CONTINUE
C
   11 CONTINUE
      PSIG = 0.5*YY*(G1-G2) + X2*(T2-APAN) - X1*(T1-APAN)
      PGAM = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
C
      PSIGX1 = -(T1-APAN)
      PSIGX2 =   T2-APAN
      PSIGYY = 0.5*(G1-G2)
      PGAMX1 = 0.5*G1
      PGAMX2 = -.5*G2
      PGAMYY = T1-T2
C
      PSIGNI = PSIGX1*X1I + PSIGX2*X2I + PSIGYY*YYI
      PGAMNI = PGAMX1*X1I + PGAMX2*X2I + PGAMYY*YYI
C
C---- TE panel source and vortex strengths
      SIGTE1 = 0.5*SCS*(GAMU(JP,1) - GAMU(JO,1))
      SIGTE2 = 0.5*SCS*(GAMU(JP,2) - GAMU(JO,2))
      GAMTE1 = -.5*SDS*(GAMU(JP,1) - GAMU(JO,1))
      GAMTE2 = -.5*SDS*(GAMU(JP,2) - GAMU(JO,2))
C
      SIGTE = 0.5*SCS*(GAM(JP) - GAM(JO))
      GAMTE = -.5*SDS*(GAM(JP) - GAM(JO))
C
C---- TE panel contribution to Psi
      PSI = PSI + HOPI*(PSIG*SIGTE + PGAM*GAMTE)
C
C---- dPsi/dGam
      DZDG(JO) = DZDG(JO) - HOPI*PSIG*SCS*0.5
      DZDG(JP) = DZDG(JP) + HOPI*PSIG*SCS*0.5
C
      DZDG(JO) = DZDG(JO) + HOPI*PGAM*SDS*0.5
      DZDG(JP) = DZDG(JP) - HOPI*PGAM*SDS*0.5
C
C---- dPsi/dni
      PSI_NI = PSI_NI + HOPI*(PSIGNI*SIGTE + PGAMNI*GAMTE)
C
      QTAN1 = QTAN1 + HOPI*(PSIGNI*SIGTE1 + PGAMNI*GAMTE1)
      QTAN2 = QTAN2 + HOPI*(PSIGNI*SIGTE2 + PGAMNI*GAMTE2)
C
      DQDG(JO) = DQDG(JO) - HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS)
      DQDG(JP) = DQDG(JP) + HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS)
C
      IF(GEOLIN) THEN
C
C----- dPsi/dn
       DZDN(JO) = DZDN(JO)
     &          + HOPI*(PSIGX1*X1O + PSIGX2*X2O + PSIGYY*YYO)*SIGTE
     &          + HOPI*(PGAMX1*X1O + PGAMX2*X2O + PGAMYY*YYO)*GAMTE
       DZDN(JP) = DZDN(JP)
     &          + HOPI*(PSIGX1*X1P + PSIGX2*X2P + PSIGYY*YYP)*SIGTE
     &          + HOPI*(PGAMX1*X1P + PGAMX2*X2P + PGAMYY*YYP)*GAMTE
C
C----- dPsi/dP
       Z_QDOF0 = Z_QDOF0 + HOPI*PSIG*0.5*(QF0(JP)-QF0(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF0(JP)-QF0(JO))*SDS
       Z_QDOF1 = Z_QDOF1 + HOPI*PSIG*0.5*(QF1(JP)-QF1(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF1(JP)-QF1(JO))*SDS
       Z_QDOF2 = Z_QDOF2 + HOPI*PSIG*0.5*(QF2(JP)-QF2(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF2(JP)-QF2(JO))*SDS
       Z_QDOF3 = Z_QDOF3 + HOPI*PSIG*0.5*(QF3(JP)-QF3(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF3(JP)-QF3(JO))*SDS
C
      ENDIF
C
   12 CONTINUE
C
C**** Freestream terms
      PSI = PSI + QINF*(COSA*YI - SINA*XI)
C
C---- dPsi/dn
      PSI_NI = PSI_NI + QINF*(COSA*NYI - SINA*NXI)
C
      QTAN1 = QTAN1 + QINF*NYI
      QTAN2 = QTAN2 - QINF*NXI
C
C---- dPsi/dQinf
      Z_QINF = Z_QINF + (COSA*YI - SINA*XI)
C
C---- dPsi/dalfa
      Z_ALFA = Z_ALFA - QINF*(SINA*YI + COSA*XI)
C
      IF(.NOT.LIMAGE) RETURN
C
C
C
      DO 20 JO=1, N
        JP = JO+1
C
        JM = JO-1
        JQ = JP+1
C
        IF(JO.ceq.1) THEN
         JM = JO
        ELSE IF(JO.ceq.N-1) THEN
         JQ = JP
        ELSE IF(JO.ceq.N) THEN
         JP = 1
         IF((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2 .LT. SEPS**2) GO TO 22
        ENDIF
C
        DSO = SQRT((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2)
C
C------ skip null panel
        IF(DSO .ceq. 0.0) GO TO 20
C
        DSIO = 1.0 / DSO
C
ccc     APAN = APANEL(JO)
        APAN = PI - APANEL(JO) + 2.0*ALFA
C
        XJO = X(JO) + 2.0*(YIMAGE+Y(JO))*SINA
        YJO = Y(JO) - 2.0*(YIMAGE+Y(JO))*COSA
        XJP = X(JP) + 2.0*(YIMAGE+Y(JP))*SINA
        YJP = Y(JP) - 2.0*(YIMAGE+Y(JP))*COSA
C
        RX1 = XI - XJO
        RY1 = YI - YJO
        RX2 = XI - XJP
        RY2 = YI - YJP
C
        SX = (XJP - XJO) * DSIO
        SY = (YJP - YJO) * DSIO
C
        X1 = SX*RX1 + SY*RY1
        X2 = SX*RX2 + SY*RY2
        YY = SX*RY1 - SY*RX1
C
        RS1 = RX1*RX1 + RY1*RY1
        RS2 = RX2*RX2 + RY2*RY2
C
C------ set reflection flag SGN to avoid branch problems with arctan
        IF(IO.GE.1 .AND. IO.LE.N) THEN
C------- no problem on airfoil surface
         SGN = 1.0
        ELSE
C------- make sure arctan falls between  -/+  Pi/2
         SGN = SIGN(1.0,YY)
        ENDIF
C
C------ set log(r^2) and arctan(x/y), correcting for reflection if any
        G1 = LOG(RS1)
        T1 = ATAN2(SGN*X1,SGN*YY) + (0.5 - 0.5*SGN)*PI
C
        G2 = LOG(RS2)
        T2 = ATAN2(SGN*X2,SGN*YY) + (0.5 - 0.5*SGN)*PI
C
        X1I = SX*NXI + SY*NYI
        X2I = SX*NXI + SY*NYI
        YYI = SX*NYI - SY*NXI
C
        IF(GEOLIN) THEN
         NXO = NX(JO)
         NYO = NY(JO)
         NXP = NX(JP)
         NYP = NY(JP)
C
         X1O =-((RX1-X1*SX)*NXO + (RY1-X1*SY)*NYO)*DSIO-(SX*NXO+SY*NYO)
         X1P = ((RX1-X1*SX)*NXP + (RY1-X1*SY)*NYP)*DSIO
         X2O =-((RX2-X2*SX)*NXO + (RY2-X2*SY)*NYO)*DSIO
         X2P = ((RX2-X2*SX)*NXP + (RY2-X2*SY)*NYP)*DSIO-(SX*NXP+SY*NYP)
         YYO = ((RX1+X1*SY)*NYO - (RY1-X1*SX)*NXO)*DSIO-(SX*NYO-SY*NXO)
         YYP =-((RX1-X1*SY)*NYP - (RY1+X1*SX)*NXP)*DSIO
        ENDIF
C
        IF(JO.ceq.N) GO TO 21
C
        IF(SIGLIN) THEN
C
C------- set up midpoint quantities
         X0 = 0.5*(X1+X2)
         RS0 = X0*X0 + YY*YY
         G0 = LOG(RS0)
         T0 = ATAN2(SGN*X0,SGN*YY) + (0.5 - 0.5*SGN)*PI
C
C------- calculate source contribution to Psi  for  1-0  half-panel
         DXINV = 1.0/(X1-X0)
         PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0)
         PDIF = ((X1+X0)*PSUM + RS1*(T1-APAN) - RS0*(T0-APAN)
     &        + (X0-X1)*YY) * DXINV
C
         PSX1 =  -(T1-APAN)
         PSX0 =    T0-APAN
         PSYY =  0.5*(G1-G0)
C
         PDX1 = ((X1+X0)*PSX1 + PSUM + 2.0*X1*(T1-APAN) - PDIF) * DXINV
         PDX0 = ((X1+X0)*PSX0 + PSUM - 2.0*X0*(T0-APAN) + PDIF) * DXINV
         PDYY = ((X1+X0)*PSYY + 2.0*(X0-X1 + YY*(T1-T0))      ) * DXINV
C
         DSM = SQRT((X(JP)-X(JM))**2 + (Y(JP)-Y(JM))**2)
         DSIM = 1.0/DSM
C
CCC      SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC      SIG1 = (SIG(JP) - SIG(JM))*DSIM
CCC      SSUM = SIG0 + SIG1
CCC      SDIF = SIG0 - SIG1
C
         SSUM = (SIG(JP) - SIG(JO))*DSIO + (SIG(JP) - SIG(JM))*DSIM
         SDIF = (SIG(JP) - SIG(JO))*DSIO - (SIG(JP) - SIG(JM))*DSIM
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JM) = DZDM(JM) + QOPI*(-PSUM*DSIM + PDIF*DSIM)
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*DSIO - PDIF*DSIO)
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*(DSIO+DSIM)
     &                                          + PDIF*(DSIO-DSIM))
C
C------- dPsi/dni
         PSNI = PSX1*X1I + PSX0*(X1I+X2I)*0.5 + PSYY*YYI
         PDNI = PDX1*X1I + PDX0*(X1I+X2I)*0.5 + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JM) = DQDM(JM) + QOPI*(-PSNI*DSIM + PDNI*DSIM)
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*DSIO - PDNI*DSIO)
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*(DSIO+DSIM)
     &                                          + PDNI*(DSIO-DSIM))
C
C
C------- calculate source contribution to Psi  for  0-2  half-panel
         DXINV = 1.0/(X0-X2)
         PSUM = X2*(T2-APAN) - X0*(T0-APAN) + 0.5*YY*(G0-G2)
         PDIF = ((X0+X2)*PSUM + RS0*(T0-APAN) - RS2*(T2-APAN)
     &        + (X2-X0)*YY) * DXINV
C
         PSX0 =  -(T0-APAN)
         PSX2 =    T2-APAN
         PSYY =  0.5*(G0-G2)
C
         PDX0 = ((X0+X2)*PSX0 + PSUM + 2.0*X0*(T0-APAN) - PDIF) * DXINV
         PDX2 = ((X0+X2)*PSX2 + PSUM - 2.0*X2*(T2-APAN) + PDIF) * DXINV
         PDYY = ((X0+X2)*PSYY + 2.0*(X2-X0 + YY*(T0-T2))      ) * DXINV
C
         DSP = SQRT((X(JQ)-X(JO))**2 + (Y(JQ)-Y(JO))**2)
         DSIP = 1.0/DSP
C
CCC         SIG2 = (SIG(JQ) - SIG(JO))*DSIP
CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC         SSUM = SIG2 + SIG0
CCC         SDIF = SIG2 - SIG0
C
         SSUM = (SIG(JQ) - SIG(JO))*DSIP + (SIG(JP) - SIG(JO))*DSIO
         SDIF = (SIG(JQ) - SIG(JO))*DSIP - (SIG(JP) - SIG(JO))*DSIO
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*(DSIP+DSIO)
     &                                          - PDIF*(DSIP-DSIO))
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*DSIO - PDIF*DSIO)
         DZDM(JQ) = DZDM(JQ) + QOPI*( PSUM*DSIP + PDIF*DSIP)
C
C------- dPsi/dni
         PSNI = PSX0*(X1I+X2I)*0.5 + PSX2*X2I + PSYY*YYI
         PDNI = PDX0*(X1I+X2I)*0.5 + PDX2*X2I + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*(DSIP+DSIO)
     &                                          - PDNI*(DSIP-DSIO))
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*DSIO - PDNI*DSIO)
         DQDM(JQ) = DQDM(JQ) + QOPI*( PSNI*DSIP + PDNI*DSIP)
C
        ENDIF
C
C------ calculate vortex panel contribution to Psi
        DXINV = 1.0/(X1-X2)
        PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
        PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2-RS1*G1 + X1*X1-X2*X2))*DXINV
C
        PSX1 = 0.5*G1
        PSX2 = -.5*G2
        PSYY = T1-T2
C
        PDX1 = ((X1+X2)*PSX1 + PSIS - X1*G1 - PSID)*DXINV
        PDX2 = ((X1+X2)*PSX2 + PSIS + X2*G2 + PSID)*DXINV
        PDYY = ((X1+X2)*PSYY - YY*(G1-G2)         )*DXINV
C
        GSUM1 = GAMU(JP,1) + GAMU(JO,1)
        GSUM2 = GAMU(JP,2) + GAMU(JO,2)
        GDIF1 = GAMU(JP,1) - GAMU(JO,1)
        GDIF2 = GAMU(JP,2) - GAMU(JO,2)
C
        GSUM = GAM(JP) + GAM(JO)
        GDIF = GAM(JP) - GAM(JO)
C
        PSI = PSI - QOPI*(PSIS*GSUM + PSID*GDIF)
C
C------ dPsi/dGam
        DZDG(JO) = DZDG(JO) - QOPI*(PSIS-PSID)
        DZDG(JP) = DZDG(JP) - QOPI*(PSIS+PSID)
C
C------ dPsi/dni
        PSNI = PSX1*X1I + PSX2*X2I + PSYY*YYI
        PDNI = PDX1*X1I + PDX2*X2I + PDYY*YYI
        PSI_NI = PSI_NI - QOPI*(GSUM*PSNI + GDIF*PDNI)
C
        QTAN1 = QTAN1 - QOPI*(GSUM1*PSNI + GDIF1*PDNI)
        QTAN2 = QTAN2 - QOPI*(GSUM2*PSNI + GDIF2*PDNI)
C
        DQDG(JO) = DQDG(JO) - QOPI*(PSNI - PDNI)
        DQDG(JP) = DQDG(JP) - QOPI*(PSNI + PDNI)
C
        IF(GEOLIN) THEN
C
C------- dPsi/dn
         DZDN(JO) = DZDN(JO)- QOPI*GSUM*(PSX1*X1O + PSX2*X2O + PSYY*YYO)
     &                      - QOPI*GDIF*(PDX1*X1O + PDX2*X2O + PDYY*YYO)
         DZDN(JP) = DZDN(JP)- QOPI*GSUM*(PSX1*X1P + PSX2*X2P + PSYY*YYP)
     &                      - QOPI*GDIF*(PDX1*X1P + PDX2*X2P + PDYY*YYP)
C------- dPsi/dP
         Z_QDOF0 = Z_QDOF0
     &           - QOPI*((PSIS-PSID)*QF0(JO) + (PSIS+PSID)*QF0(JP))
         Z_QDOF1 = Z_QDOF1
     &           - QOPI*((PSIS-PSID)*QF1(JO) + (PSIS+PSID)*QF1(JP))
         Z_QDOF2 = Z_QDOF2
     &           - QOPI*((PSIS-PSID)*QF2(JO) + (PSIS+PSID)*QF2(JP))
         Z_QDOF3 = Z_QDOF3
     &           - QOPI*((PSIS-PSID)*QF3(JO) + (PSIS+PSID)*QF3(JP))
        ENDIF
C
C
   20 CONTINUE
C
   21 CONTINUE
      PSIG = 0.5*YY*(G1-G2) + X2*(T2-APAN) - X1*(T1-APAN)
      PGAM = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
C
      PSIGX1 = -(T1-APAN)
      PSIGX2 =   T2-APAN
      PSIGYY = 0.5*(G1-G2)
      PGAMX1 = 0.5*G1
      PGAMX2 = -.5*G2
      PGAMYY = T1-T2
C
      PSIGNI = PSIGX1*X1I + PSIGX2*X2I + PSIGYY*YYI
      PGAMNI = PGAMX1*X1I + PGAMX2*X2I + PGAMYY*YYI
C
C---- TE panel source and vortex strengths
      SIGTE1 = 0.5*SCS*(GAMU(JP,1) - GAMU(JO,1))
      SIGTE2 = 0.5*SCS*(GAMU(JP,2) - GAMU(JO,2))
      GAMTE1 = -.5*SDS*(GAMU(JP,1) - GAMU(JO,1))
      GAMTE2 = -.5*SDS*(GAMU(JP,2) - GAMU(JO,2))
C
      SIGTE = 0.5*SCS*(GAM(JP) - GAM(JO))
      GAMTE = -.5*SDS*(GAM(JP) - GAM(JO))
C
C---- TE panel contribution to Psi
      PSI = PSI + HOPI*(PSIG*SIGTE - PGAM*GAMTE)
C
C---- dPsi/dGam
      DZDG(JO) = DZDG(JO) - HOPI*PSIG*SCS*0.5
      DZDG(JP) = DZDG(JP) + HOPI*PSIG*SCS*0.5
C
      DZDG(JO) = DZDG(JO) - HOPI*PGAM*SDS*0.5
      DZDG(JP) = DZDG(JP) + HOPI*PGAM*SDS*0.5
C
C---- dPsi/dni
      PSI_NI = PSI_NI + HOPI*(PSIGNI*SIGTE - PGAMNI*GAMTE)
C
      QTAN1 = QTAN1 + HOPI*(PSIGNI*SIGTE1 - PGAMNI*GAMTE1)
      QTAN2 = QTAN2 + HOPI*(PSIGNI*SIGTE2 - PGAMNI*GAMTE2)
C
      DQDG(JO) = DQDG(JO) - HOPI*(PSIGNI*0.5*SCS + PGAMNI*0.5*SDS)
      DQDG(JP) = DQDG(JP) + HOPI*(PSIGNI*0.5*SCS + PGAMNI*0.5*SDS)
C
      IF(GEOLIN) THEN
C
C----- dPsi/dn
       DZDN(JO) = DZDN(JO)
     &          + HOPI*(PSIGX1*X1O + PSIGX2*X2O + PSIGYY*YYO)*SIGTE
     &          - HOPI*(PGAMX1*X1O + PGAMX2*X2O + PGAMYY*YYO)*GAMTE
       DZDN(JP) = DZDN(JP)
     &          + HOPI*(PSIGX1*X1P + PSIGX2*X2P + PSIGYY*YYP)*SIGTE
     &          - HOPI*(PGAMX1*X1P + PGAMX2*X2P + PGAMYY*YYP)*GAMTE
C
C----- dPsi/dP
       Z_QDOF0 = Z_QDOF0 + HOPI*PSIG*0.5*(QF0(JP)-QF0(JO))*SCS
     &                   + HOPI*PGAM*0.5*(QF0(JP)-QF0(JO))*SDS
       Z_QDOF1 = Z_QDOF1 + HOPI*PSIG*0.5*(QF1(JP)-QF1(JO))*SCS
     &                   + HOPI*PGAM*0.5*(QF1(JP)-QF1(JO))*SDS
       Z_QDOF2 = Z_QDOF2 + HOPI*PSIG*0.5*(QF2(JP)-QF2(JO))*SCS
     &                   + HOPI*PGAM*0.5*(QF2(JP)-QF2(JO))*SDS
       Z_QDOF3 = Z_QDOF3 + HOPI*PSIG*0.5*(QF3(JP)-QF3(JO))*SCS
     &                   + HOPI*PGAM*0.5*(QF3(JP)-QF3(JO))*SDS
C
      ENDIF
C
   22 CONTINUE
C
      RETURN
      END


      SUBROUTINE PSWLIN(I,XI,YI,NXI,NYI,PSI,PSI_NI)
C--------------------------------------------------------------------
C     Calculates current streamfunction Psi and tangential velocity
C     Qtan at panel node or wake node I due to freestream and wake
C     sources Sig.  Also calculates sensitivity vectors dPsi/dSig
C     (DZDM) and dQtan/dSig (DQDM).
C
C          Airfoil:  1   < I < N
C          Wake:     N+1 < I < N+NW
C--------------------------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
      complex NXI, NYI
C
      IO = I
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
      DO 4 JO=N+1, N+NW
        DZDM(JO) = 0.0
        DQDM(JO) = 0.0
    4 CONTINUE
C
      PSI    = 0.
      PSI_NI = 0.
C
      DO 20 JO=N+1, N+NW-1
C
        JP = JO+1
C
        JM = JO-1
        JQ = JP+1
        IF(JO.ceq.N+1) THEN
         JM = JO
        ELSE IF(JO.ceq.N+NW-1) THEN
         JQ = JP
        ENDIF
C
        DSO = SQRT((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2)
        DSIO = 1.0 / DSO
C
        APAN = APANEL(JO)
C
        RX1 = XI - X(JO)
        RY1 = YI - Y(JO)
        RX2 = XI - X(JP)
        RY2 = YI - Y(JP)
C
        SX = (X(JP) - X(JO)) * DSIO
        SY = (Y(JP) - Y(JO)) * DSIO
C
        X1 = SX*RX1 + SY*RY1
        X2 = SX*RX2 + SY*RY2
        YY = SX*RY1 - SY*RX1
C
        RS1 = RX1*RX1 + RY1*RY1
        RS2 = RX2*RX2 + RY2*RY2
C
        IF(IO.GE.N+1 .AND. IO.LE.N+NW) THEN
         SGN = 1.0
        ELSE
         SGN = SIGN(1.0,YY)
        ENDIF
C
        IF((IO.cne.JO) .AND. RS1.GT.0.0) THEN
         G1 = LOG(RS1)
         T1 = ATAN2(SGN*X1,SGN*YY) - (0.5 - 0.5*SGN)*PI
        ELSE
         G1 = 0.0
         T1 = 0.0
        ENDIF
C
        IF((IO.cne.JP) .AND. RS2.GT.0.0) THEN
         G2 = LOG(RS2)
         T2 = ATAN2(SGN*X2,SGN*YY) - (0.5 - 0.5*SGN)*PI
        ELSE
         G2 = 0.0
         T2 = 0.0
        ENDIF
C
        X1I = SX*NXI + SY*NYI
        X2I = SX*NXI + SY*NYI
        YYI = SX*NYI - SY*NXI
C
C------- set up midpoint quantities
         X0 = 0.5*(X1+X2)
         RS0 = X0*X0 + YY*YY
         G0 = LOG(RS0)
         T0 = ATAN2(SGN*X0,SGN*YY) - (0.5 - 0.5*SGN)*PI
C
C------- calculate source contribution to Psi  for  1-0  half-panel
         DXINV = 1.0/(X1-X0)
         PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0)
         PDIF = ((X1+X0)*PSUM + RS1*(T1-APAN) - RS0*(T0-APAN)
     &        + (X0-X1)*YY) * DXINV
C
         PSX1 =  -(T1-APAN)
         PSX0 =    T0-APAN
         PSYY =  0.5*(G1-G0)
C
         PDX1 = ((X1+X0)*PSX1 + PSUM + 2.0*X1*(T1-APAN) - PDIF) * DXINV
         PDX0 = ((X1+X0)*PSX0 + PSUM - 2.0*X0*(T0-APAN) + PDIF) * DXINV
         PDYY = ((X1+X0)*PSYY + 2.0*(X0-X1 + YY*(T1-T0))      ) * DXINV
C
         DSM = SQRT((X(JP)-X(JM))**2 + (Y(JP)-Y(JM))**2)
         DSIM = 1.0/DSM
C
CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC         SIG1 = (SIG(JP) - SIG(JM))*DSIM
CCC         SSUM = SIG0 + SIG1
CCC         SDIF = SIG0 - SIG1
C
         SSUM = (SIG(JP) - SIG(JO))*DSIO + (SIG(JP) - SIG(JM))*DSIM
         SDIF = (SIG(JP) - SIG(JO))*DSIO - (SIG(JP) - SIG(JM))*DSIM
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JM) = DZDM(JM) + QOPI*(-PSUM*DSIM + PDIF*DSIM)
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*DSIO - PDIF*DSIO)
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*(DSIO+DSIM)
     &                                          + PDIF*(DSIO-DSIM))
C
C------- dPsi/dni
         PSNI = PSX1*X1I + PSX0*(X1I+X2I)*0.5 + PSYY*YYI
         PDNI = PDX1*X1I + PDX0*(X1I+X2I)*0.5 + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JM) = DQDM(JM) + QOPI*(-PSNI*DSIM + PDNI*DSIM)
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*DSIO - PDNI*DSIO)
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*(DSIO+DSIM)
     &                                          + PDNI*(DSIO-DSIM))
C
C
C------- calculate source contribution to Psi  for  0-2  half-panel
         DXINV = 1.0/(X0-X2)
         PSUM = X2*(T2-APAN) - X0*(T0-APAN) + 0.5*YY*(G0-G2)
         PDIF = ((X0+X2)*PSUM + RS0*(T0-APAN) - RS2*(T2-APAN)
     &        + (X2-X0)*YY) * DXINV
C
         PSX0 =  -(T0-APAN)
         PSX2 =    T2-APAN
         PSYY =  0.5*(G0-G2)
C
         PDX0 = ((X0+X2)*PSX0 + PSUM + 2.0*X0*(T0-APAN) - PDIF) * DXINV
         PDX2 = ((X0+X2)*PSX2 + PSUM - 2.0*X2*(T2-APAN) + PDIF) * DXINV
         PDYY = ((X0+X2)*PSYY + 2.0*(X2-X0 + YY*(T0-T2))      ) * DXINV
C
         DSP = SQRT((X(JQ)-X(JO))**2 + (Y(JQ)-Y(JO))**2)
         DSIP = 1.0/DSP
C
CCC         SIG2 = (SIG(JQ) - SIG(JO))*DSIP
CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC         SSUM = SIG2 + SIG0
CCC         SDIF = SIG2 - SIG0
C
         SSUM = (SIG(JQ) - SIG(JO))*DSIP + (SIG(JP) - SIG(JO))*DSIO
         SDIF = (SIG(JQ) - SIG(JO))*DSIP - (SIG(JP) - SIG(JO))*DSIO
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*(DSIP+DSIO)
     &                                          - PDIF*(DSIP-DSIO))
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*DSIO - PDIF*DSIO)
         DZDM(JQ) = DZDM(JQ) + QOPI*( PSUM*DSIP + PDIF*DSIP)
C
C------- dPsi/dni
         PSNI = PSX0*(X1I+X2I)*0.5 + PSX2*X2I + PSYY*YYI
         PDNI = PDX0*(X1I+X2I)*0.5 + PDX2*X2I + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*(DSIP+DSIO)
     &                                          - PDNI*(DSIP-DSIO))
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*DSIO - PDNI*DSIO)
         DQDM(JQ) = DQDM(JQ) + QOPI*( PSNI*DSIP + PDNI*DSIP)
C
   20 CONTINUE
C
      RETURN
      END



C**************************************************
C
C	CALLED FROM SPECAL/SPECCL
C
C**************************************************
      SUBROUTINE GGCALC
C--------------------------------------------------------------
C     Calculates two surface vorticity (gamma) distributions
C     for alpha = 0, 90  degrees.  These are superimposed
C     in SPECAL or SPECCL for specified alpha or CL.
C--------------------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
C---- distance of internal control point ahead of sharp TE
C-    (fraction of smaller panel length adjacent to TE)
      BWT = 0.1
C
C      WRITE(*,*) 'Calculating unit vorticity distributions ...'
C
      DO 10 I=1, N
        GAM(I) = 0.
        GAMU(I,1) = 0.
        GAMU(I,2) = 0.
   10 CONTINUE
      PSIO = 0.
C
C---- Set up matrix system for  Psi = Psio  on airfoil surface.
C-    The unknowns are (dGamma)i and dPsio.
      DO 20 I=1, N
C
C------ calculate Psi and dPsi/dGamma array for current node
        CALL PSILIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N,.FALSE.,.TRUE.)
C
        PSIINF = QINF*(COS(ALFA)*Y(I) - SIN(ALFA)*X(I))
C
C------ RES1 = PSI( 0) - PSIO
C------ RES2 = PSI(90) - PSIO
        RES1 =  QINF*Y(I)
        RES2 = -QINF*X(I)
C
C------ dRes/dGamma
        DO 201 J=1, N
          AIJ(I,J) = DZDG(J)
  201   CONTINUE
C
        DO 202 J=1, N
          BIJ(I,J) = -DZDM(J)
  202   CONTINUE
C
C------ dRes/dPsio
        AIJ(I,N+1) = -1.0
C
        GAMU(I,1) = -RES1
        GAMU(I,2) = -RES2
C
   20 CONTINUE
C
C---- set Kutta condition
C-    RES = GAM(1) + GAM(N)
      RES = 0.
C
      DO 30 J=1, N+1
        AIJ(N+1,J) = 0.0
   30 CONTINUE
C
      AIJ(N+1,1) = 1.0
      AIJ(N+1,N) = 1.0
C
      GAMU(N+1,1) = -RES
      GAMU(N+1,2) = -RES
C
C---- set up Kutta condition (no direct source influence)
      DO 32 J=1, N
        BIJ(N+1,J) = 0.
   32 CONTINUE
C
      IF(SHARP) THEN
C----- set zero internal velocity in TE corner
C
C----- set TE bisector angle
       AG1 = ATAN2(-YP(1),-XP(1)    )
       AG2 = ATANC( YP(N), XP(N),AG1)
       ABIS = 0.5*(AG1+AG2)
       CBIS = COS(ABIS)
       SBIS = SIN(ABIS)
C
C----- minimum panel length adjacent to TE
       DS1 = SQRT( (X(1)-X(2)  )**2 + (Y(1)-Y(2)  )**2 )
       DS2 = SQRT( (X(N)-X(N-1))**2 + (Y(N)-Y(N-1))**2 )
       DSMIN = MIN( DS1 , DS2 )
C
C----- control point on bisector just ahead of TE point
       XBIS = XTE - BWT*DSMIN*CBIS
       YBIS = YTE - BWT*DSMIN*SBIS
ccc       write(*,*) xbis, ybis
C
C----- set velocity component along bisector line
       CALL PSILIN(0,XBIS,YBIS,-SBIS,CBIS,PSI,QBIS,.FALSE.,.TRUE.)
C
CCC--- RES = DQDGj*Gammaj + DQDMj*Massj + QINF*(COSA*CBIS + SINA*SBIS)
       RES = QBIS
C
C----- dRes/dGamma
       DO J=1, N
         AIJ(N,J) = DQDG(J)
       ENDDO
C
C----- -dRes/dMass
       DO J=1, N
         BIJ(N,J) = -DQDM(J)
       ENDDO
C
C----- dRes/dPsio
       AIJ(N,N+1) = 0.
C
C----- -dRes/dUinf
       GAMU(N,1) = -CBIS
C
C----- -dRes/dVinf
       GAMU(N,2) = -SBIS
C
      ENDIF
C
C---- LU-factor coefficient matrix AIJ
      CALL LUDCMP(IQX,N+1,AIJ,AIJPIV)
      LQAIJ = .TRUE.
C
C---- solve system for the two vorticity distributions
      CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,GAMU(1,1))
      CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,GAMU(1,2))
C
C---- set inviscid alpha=0,90 surface speeds for this geometry
      DO 50 I=1, N
        QINVU(I,1) = GAMU(I,1)
        QINVU(I,2) = GAMU(I,2)
   50 CONTINUE
C
      LGAMU = .TRUE.
C
      RETURN
      END



      SUBROUTINE QWCALC
C---------------------------------------------------------------
C     Sets inviscid tangential velocity for alpha = 0, 90
C     on wake due to freestream and airfoil surface vorticity.
C---------------------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
C---- first wake point (same as TE)
      QINVU(N+1,1) = QINVU(N,1)
      QINVU(N+1,2) = QINVU(N,2)
C
C---- rest of wake
      DO 10 I=N+2, N+NW
        CALL PSILIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_NI,.FALSE.,.FALSE.)
        QINVU(I,1) = QTAN1
        QINVU(I,2) = QTAN2
   10 CONTINUE
C
      RETURN
      END


      SUBROUTINE QDCALC
C-----------------------------------------------------
C     Calculates source panel influence coefficient
C     matrix for current airfoil and wake geometry.
C-----------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
C      WRITE(*,*) 'Calculating source influence matrix ...'
C
      IF(.NOT.LADIJ) THEN
C
C----- calculate source influence matrix for airfoil surface if it doesn't exist
       DO 10 J=1, N
C
C------- multiply each dPsi/Sig vector by inverse of factored dPsi/dGam matrix
         CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,BIJ(1,J))
C
C------- store resulting dGam/dSig = dQtan/dSig vector
         DO 105 I=1, N
           DIJ(I,J) = BIJ(I,J)
  105    CONTINUE
C
   10  CONTINUE
       LADIJ = .TRUE.
C
      ENDIF
C
C---- set up coefficient matrix of dPsi/dm on airfoil surface
      DO 20 I=1, N
        CALL PSWLIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N)
        DO 202 J=N+1, N+NW
          BIJ(I,J) = -DZDM(J)
  202   CONTINUE
   20 CONTINUE
C
C---- set up Kutta condition (no direct source influence)
      DO 32 J=N+1, N+NW
        BIJ(N+1,J) = 0.
   32 CONTINUE
C
C---- sharp TE gamma extrapolation also has no source influence
      IF(SHARP) THEN
       DO 34 J=N+1, N+NW
         BIJ(N,J) = 0.
   34  CONTINUE
      ENDIF
C
C---- multiply by inverse of factored dPsi/dGam matrix
      DO 40 J=N+1, N+NW
        CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,BIJ(1,J))
   40 CONTINUE
C
C---- set the source influence matrix for the wake sources
      DO 50 I=1, N
        DO 510 J=N+1, N+NW
          DIJ(I,J) = BIJ(I,J)
  510   CONTINUE
   50 CONTINUE
C
C**** Now we need to calculate the influence of sources on the wake velocities
C
C---- calculcate dQtan/dGam and dQtan/dSig at the wake points
      DO 70 I=N+1, N+NW
C
        IW = I-N
C
C------ airfoil contribution at wake panel node
        CALL PSILIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N,.FALSE.,.TRUE.)
C
        DO 710 J=1, N
          CIJ(IW,J) = DQDG(J)
  710   CONTINUE
C
        DO 720 J=1, N
          DIJ(I,J) = DQDM(J)
  720   CONTINUE
C
C------ wake contribution
        CALL PSWLIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N)
C
        DO 730 J=N+1, N+NW
          DIJ(I,J) = DQDM(J)
  730   CONTINUE
C
   70 CONTINUE
C
C---- add on effect of all sources on airfoil vorticity which effects wake Qtan
      DO 80 I=N+1, N+NW
        IW = I-N
C
C------ airfoil surface source contribution first
        DO 810 J=1, N
          SUM = 0.
          DO 8100 K=1, N
            SUM = SUM + CIJ(IW,K)*DIJ(K,J)
 8100     CONTINUE
          DIJ(I,J) = DIJ(I,J) + SUM
  810   CONTINUE
C
C------ wake source contribution next
        DO 820 J=N+1, N+NW
          SUM = 0.
          DO 8200 K=1, N
            SUM = SUM + CIJ(IW,K)*BIJ(K,J)
 8200     CONTINUE
          DIJ(I,J) = DIJ(I,J) + SUM
  820   CONTINUE
C
   80 CONTINUE
C
C---- make sure first wake point has same velocity as trailing edge
      DO 90 J=1, N+NW
        DIJ(N+1,J) = DIJ(N,J)
   90 CONTINUE
C
      LWDIJ = .TRUE.
C
      RETURN
      END


      SUBROUTINE XYWAKE
C-----------------------------------------------------
C     Sets wake coordinate array for current surface
C     vorticity and/or mass source distributions.
C-----------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
C      WRITE(*,*) 'Calculating wake trajectory ...'
C
C---- number of wake points
      NW = N/8 + 2
      IF(NW.GT.IWX) THEN
C       WRITE(*,*)
C     &  'Array size (IWX) too small.  Last wake point index reduced.'
       NW = IWX
      ENDIF
C
      DS1 = 0.5*(S(2) - S(1) + S(N) - S(N-1))
      CALL SETEXP(SNEW(N+1),DS1,WAKLEN*CHORD,NW)

      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
C
C---- set first wake point a tiny distance behind TE
      I = N+1
      SX = 0.5*(YP(N) - YP(1))
      SY = 0.5*(XP(1) - XP(N))
      SMOD = SQRT(SX**2 + SY**2)
      NX(I) = SX / SMOD
      NY(I) = SY / SMOD
      X(I) = XTE - 0.0001*NY(I)
      Y(I) = YTE + 0.0001*NX(I)
      S(I) = S(N)
C
C---- calculate streamfunction gradient components at first point

      xhat = cmplx(1,0)
      yhat = cmplx(0,0)
      CALL PSILIN(I,X(I),Y(I),xhat,yhat,PSI,PSI_X,.FALSE.,.FALSE.)
      xhat = cmplx(0,0)
      yhat = cmplx(1,0)

      CALL PSILIN(I,X(I),Y(I),xhat,yhat,PSI,PSI_Y,.FALSE.,.FALSE.)
c      write(*,*) 'PSI',PSI
c      write(*,*) 'PSI_X',PSI_X

c      write(*,*) 'PSI_Y',PSI_Y

C---- set unit vector normal to wake at first point
      NX(I+1) = -PSI_X / SQRT(PSI_X**2 + PSI_Y**2)
      NY(I+1) = -PSI_Y / SQRT(PSI_X**2 + PSI_Y**2)
C
C---- set angle of wake panel normal
      APANEL(I) = ATAN2( PSI_Y , PSI_X )
C
C---- set rest of wake points
      DO 10 I=N+2, N+NW
        DS = SNEW(I) - SNEW(I-1)
C
C------ set new point DS downstream of last point
        X(I) = X(I-1) - DS*NY(I)
        Y(I) = Y(I-1) + DS*NX(I)
        S(I) = S(I-1) + DS
C
        IF(I.ceq.N+NW) GO TO 10
C
C------- calculate normal vector for next point

         xhat = cmplx(1,0)
         yhat = cmplx(0,0)

         CALL PSILIN(I,X(I),Y(I),xhat,yhat,PSI,PSI_X,.FALSE.,.FALSE.)
         xhat = cmplx(0,0)
         yhat = cmplx(1,0)

         CALL PSILIN(I,X(I),Y(I),xhat,yhat,PSI,PSI_Y,.FALSE.,.FALSE.)

         NX(I+1) = -PSI_X / SQRT(PSI_X**2 + PSI_Y**2)
         NY(I+1) = -PSI_Y / SQRT(PSI_X**2 + PSI_Y**2)
C
C------- set angle of wake panel normal
         APANEL(I) = ATAN2( PSI_Y , PSI_X )
C
   10 CONTINUE
C
C---- set wake presence flag and corresponding alpha
      LWAKE = .TRUE.
      AWAKE =  ALFA
C
C---- old source influence matrix is invalid for the new wake geometry
      LWDIJ = .FALSE.
C
      RETURN
      END



      SUBROUTINE STFIND
C-----------------------------------------
C     Locates stagnation point arc length
C     location SST and panel index IST.
C-----------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
      DO 10 I=1, N-1
        IF(GAM(I).GE.0.0 .AND. GAM(I+1).LT.0.0) GO TO 11
   10 CONTINUE
C
C      WRITE(*,*) 'STFIND: Stagnation point not found. Continuing ...'
      I = N/2
C
   11 CONTINUE
C
      IST = I
      DGAM = GAM(I+1) - GAM(I)
      DS = S(I+1) - S(I)
C
C---- evaluate so as to minimize roundoff for very small GAM(I) or GAM(I+1)
      IF(GAM(I) .LT. -GAM(I+1)) THEN
       SST = S(I)   - DS*(GAM(I)  /DGAM)
      ELSE
       SST = S(I+1) - DS*(GAM(I+1)/DGAM)
      ENDIF
C
C---- tweak stagnation point if it falls right on a node (very unlikely)
      IF(SST .LE. S(I)  ) SST = S(I)   + 1.0E-7
      IF(SST .GE. S(I+1)) SST = S(I+1) - 1.0E-7
C
      SST_GO = (SST  - S(I+1))/DGAM
      SST_GP = (S(I) - SST   )/DGAM
C
      RETURN
      END


      SUBROUTINE IBLPAN
C-------------------------------------------------------------
C     Sets  BL location -> panel location  pointer array IPAN
C-------------------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
C---- top surface first
      IS = 1
C
      IBL = 1
      DO 10 I=IST, 1, -1
        IBL = IBL+1
        IPAN(IBL,IS) = I
        VTI(IBL,IS) = 1.0
   10 CONTINUE
C
      IBLTE(IS) = IBL
      NBL(IS) = IBL
C
C---- bottom surface next
      IS = 2
C
      IBL = 1
      DO 20 I=IST+1, N
        IBL = IBL+1
        IPAN(IBL,IS) = I
        VTI(IBL,IS) = -1.0
   20 CONTINUE
C
C---- wake
      IBLTE(IS) = IBL
C
      DO 25 IW=1, NW
        I = N+IW
        IBL = IBLTE(IS)+IW
        IPAN(IBL,IS) = I
         VTI(IBL,IS) = -1.0
   25 CONTINUE
C
      NBL(IS) = IBLTE(IS) + NW
C
C---- upper wake pointers (for plotting only)
      DO 35 IW=1, NW
        IPAN(IBLTE(1)+IW,1) = IPAN(IBLTE(2)+IW,2)
         VTI(IBLTE(1)+IW,1) = 1.0
   35 CONTINUE
C
C
      IBLMAX = MAX(IBLTE(1),IBLTE(2)) + NW
      IF(IBLMAX.GT.IVX) THEN
C        WRITE(*,*) ' ***  BL array overflow.'
C        WRITE(*,*) ' ***  Increase IVX to at least', IBLMAX
        STOP
      ENDIF
C
      LIPAN = .TRUE.
      RETURN
      END


      SUBROUTINE XICALC
C-------------------------------------------------------------
C     Sets BL arc length array on each airfoil side and wake
C-------------------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
      IS = 1
C
      XSSI(1,IS) = 0.
C
      DO 10 IBL=2, IBLTE(IS)
        I = IPAN(IBL,IS)
        XSSI(IBL,IS) = SST - S(I)
   10 CONTINUE
C
C
      IS = 2
C
      XSSI(1,IS) = 0.
C
      DO 20 IBL=2, IBLTE(IS)
        I = IPAN(IBL,IS)
        XSSI(IBL,IS) = S(I) - SST
   20 CONTINUE
C
      IBL = IBLTE(IS) + 1
      XSSI(IBL,IS) = XSSI(IBL-1,IS)
C
      DO 25 IBL=IBLTE(IS)+2, NBL(IS)
        I = IPAN(IBL,IS)
        XSSI(IBL,IS) = XSSI(IBL-1,IS)
c     &               + SQRT((X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2)
     &               + SQRT( (X(I)-X(I-1))*(X(I)-X(I-1)) +
     &                       (Y(I)-Y(I-1))*(Y(I)-Y(I-1)))


   25 CONTINUE
C


C---- trailing edge flap length to TE gap ratio
      TELRAT = 2.50
C
C---- set up parameters for TE flap cubics
C
ccc   DWDXTE = YP(1)/XP(1) + YP(N)/XP(N)    !!! BUG  2/2/95
C
      CROSP = (XP(1)*YP(N) - YP(1)*XP(N))
     &      / SQRT(  (XP(1)**2 + YP(1)**2)
     &              *(XP(N)**2 + YP(N)**2) )
      DWDXTE = CROSP / SQRT(1.0 - CROSP**2)
C
C---- limit cubic to avoid absurd TE gap widths
      DWDXTE = MAX(DWDXTE,-3.0/TELRAT)
      DWDXTE = MIN(DWDXTE, 3.0/TELRAT)
C
      AA =  3.0 + TELRAT*DWDXTE
      BB = -2.0 - TELRAT*DWDXTE
C
      IF(SHARP) THEN
       DO 30 IW=1, NW
         WGAP(IW) = 0.
   30  CONTINUE
      ELSE
C----- set TE flap (wake gap) array
       IS = 2
       DO 35 IW=1, NW
         IBL = IBLTE(IS) + IW
         ZN = 1.0 - (XSSI(IBL,IS)-XSSI(IBLTE(IS),IS)) / (TELRAT*ANTE)
         WGAP(IW) = 0.
         IF(ZN.GE.0.0) WGAP(IW) = ANTE * (AA + BB*ZN)*ZN**2
   35  CONTINUE
      ENDIF
C
      RETURN
      END


      SUBROUTINE UICALC
C--------------------------------------------------------------
C     Sets inviscid Ue from panel inviscid tangential velocity
C--------------------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
      DO 10 IS=1, 2
        UINV  (1,IS) = 0.
        UINV_A(1,IS) = 0.
        DO 110 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
          UINV  (IBL,IS) = VTI(IBL,IS)*QINV  (I)
          UINV_A(IBL,IS) = VTI(IBL,IS)*QINV_A(I)
  110   CONTINUE
   10 CONTINUE
C
      RETURN
      END


      SUBROUTINE UECALC
C--------------------------------------------------------------
C     Sets viscous Ue from panel viscous tangential velocity
C--------------------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
      DO 10 IS=1, 2
        UEDG(1,IS) = 0.
        DO 110 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
          UEDG(IBL,IS) = VTI(IBL,IS)*QVIS(I)
  110   CONTINUE
   10 CONTINUE
C
      RETURN
      END


      SUBROUTINE QVFUE
C--------------------------------------------------------------
C     Sets panel viscous tangential velocity from viscous Ue
C--------------------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
          QVIS(I) = VTI(IBL,IS)*UEDG(IBL,IS)
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END


      SUBROUTINE QISET
C-------------------------------------------------------
C     Sets inviscid panel tangential velocity for
C     current alpha.
C-------------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
      DO 5 I=1, N+NW
        QINV  (I) =  COSA*QINVU(I,1) + SINA*QINVU(I,2)
        QINV_A(I) = -SINA*QINVU(I,1) + COSA*QINVU(I,2)
    5 CONTINUE
C
      RETURN
      END


      SUBROUTINE GAMQV
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
      DO 10 I=1, N
        GAM(I)   = QVIS(I)
        GAM_A(I) = QINV_A(I)
   10 CONTINUE
C
      RETURN
      END


      SUBROUTINE STMOVE
C---------------------------------------------------
C     Moves stagnation point location to new panel.
C---------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
C---- locate new stagnation point arc length SST from GAM distribution
      ISTOLD = IST
      CALL STFIND
C
      IF(ISTOLD.ceq.IST) THEN
C
C----- recalculate new arc length array
       CALL XICALC
C
      ELSE
C
CCC       WRITE(*,*) 'STMOVE: Resetting stagnation point'
C
C----- set new BL position -> panel position  pointers
       CALL IBLPAN
C
C----- set new inviscid BL edge velocity UINV from QINV
       CALL UICALC
C
C----- recalculate new arc length array
       CALL XICALC
C
C----- set  BL position -> system line  pointers
       CALL IBLSYS
C
       IF(IST.GT.ISTOLD) THEN
C------ increase in number of points on top side (IS=1)
        IDIF = IST-ISTOLD
C
        ITRAN(1) = ITRAN(1) + IDIF
        ITRAN(2) = ITRAN(2) - IDIF
C
C------ move top side BL variables downstream
        DO 110 IBL=NBL(1), IDIF+2, -1
          CTAU(IBL,1) = CTAU(IBL-IDIF,1)
          THET(IBL,1) = THET(IBL-IDIF,1)
          DSTR(IBL,1) = DSTR(IBL-IDIF,1)
          UEDG(IBL,1) = UEDG(IBL-IDIF,1)
  110   CONTINUE
C
C------ set BL variables between old and new stagnation point
        DUDX = UEDG(IDIF+2,1)/XSSI(IDIF+2,1)
        DO 115 IBL=IDIF+1, 2, -1
          CTAU(IBL,1) = CTAU(IDIF+2,1)
          THET(IBL,1) = THET(IDIF+2,1)
          DSTR(IBL,1) = DSTR(IDIF+2,1)
          UEDG(IBL,1) = DUDX * XSSI(IBL,1)
  115   CONTINUE
C
C------ move bottom side BL variables upstream
        DO 120 IBL=2, NBL(2)
          CTAU(IBL,2) = CTAU(IBL+IDIF,2)
          THET(IBL,2) = THET(IBL+IDIF,2)
          DSTR(IBL,2) = DSTR(IBL+IDIF,2)
          UEDG(IBL,2) = UEDG(IBL+IDIF,2)
  120   CONTINUE
C
       ELSE
C------ increase in number of points on bottom side (IS=2)
        IDIF = ISTOLD-IST
C
        ITRAN(1) = ITRAN(1) - IDIF
        ITRAN(2) = ITRAN(2) + IDIF
C
C------ move bottom side BL variables downstream
        DO 210 IBL=NBL(2), IDIF+2, -1
          CTAU(IBL,2) = CTAU(IBL-IDIF,2)
          THET(IBL,2) = THET(IBL-IDIF,2)
          DSTR(IBL,2) = DSTR(IBL-IDIF,2)
          UEDG(IBL,2) = UEDG(IBL-IDIF,2)
  210   CONTINUE
C
C------ set BL variables between old and new stagnation point
        DUDX = UEDG(IDIF+2,2)/XSSI(IDIF+2,2)
        DO 215 IBL=IDIF+1, 2, -1
          CTAU(IBL,2) = CTAU(IDIF+2,2)
          THET(IBL,2) = THET(IDIF+2,2)
          DSTR(IBL,2) = DSTR(IDIF+2,2)
          UEDG(IBL,2) = DUDX * XSSI(IBL,2)
  215   CONTINUE
C
C------ move top side BL variables upstream
        DO 220 IBL=2, NBL(1)
          CTAU(IBL,1) = CTAU(IBL+IDIF,1)
          THET(IBL,1) = THET(IBL+IDIF,1)
          DSTR(IBL,1) = DSTR(IBL+IDIF,1)
          UEDG(IBL,1) = UEDG(IBL+IDIF,1)
  220   CONTINUE
       ENDIF
C
      ENDIF
C
C---- set new mass array since Ue has been tweaked
      DO 50 IS=1, 2
        DO 510 IBL=2, NBL(IS)
          MASS(IBL,IS) = DSTR(IBL,IS)*UEDG(IBL,IS)
  510   CONTINUE
   50 CONTINUE
C
      RETURN
      END


      SUBROUTINE UESET
C---------------------------------------------------------
C     Sets Ue from inviscid Ue plus all source influence
C---------------------------------------------------------
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
C
          DUI = 0.
          DO 100 JS=1, 2
            DO 1000 JBL=2, NBL(JS)
              J  = IPAN(JBL,JS)
              UE_M = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
              DUI = DUI + UE_M*MASS(JBL,JS)
 1000       CONTINUE
  100     CONTINUE
C
          UEDG(IBL,IS) = UINV(IBL,IS) + DUI
C
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END


      SUBROUTINE DSSET
	use complexify
	implicit complex(a-h, o-z)
      include 'c_XFOIL.INC'
C
      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          DSTR(IBL,IS) = MASS(IBL,IS) / UEDG(IBL,IS)
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END
