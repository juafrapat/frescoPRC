MODULE modulo
  !***********************************************************************************************
  !* This module contains all functions and subroutines for frescoPRC program:                   *
  !*                                                                                             *
  !* Subroutine CLEBSCH: Compute C-G coefficients                                                *
  !*                                                                                             *
  !* Subroutine wigner: Compute matrix element of Wigner's functions accoding to B.15 from       *
  !*                    PRC 94 6 (2016), 064605 (even-even nuclei).                                                 *
  !*                                                                                             *
  !* Dispersive functions pack: Calculation of Analytical dispersive integrals.                  *
  !*                                                                                             *
  !* Subroutine FORMFACT: Numerical calculation of form factors using                            *
  !*                      Gauss-Legendre quadrature. Those form factors are part of              *
  !*                      FRESCO's input [fort.4].                                               *
  !*                                                                                             *
  !* Subroutine steps: Build the couplings between G.S band and excited bands for                *
  !*                   even-even nuclei according to equation B.12 from PRC 94 6 (2016), 064605. *
  !*                                                                                             *
  !***********************************************************************************************
CONTAINS

  SUBROUTINE CLEBSCH(AJ,BJ,CJ,AM,BM,CM,CG)

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION Q(100,100)
    DOUBLE PRECISION CG
    DOUBLE PRECISION AJ,BJ,CJ,AM,BM,CM
    INTEGER ZZ
    ZZ=MAX(2*AJ+1,2*BJ+1,2*CJ+1,AJ+BJ+CJ,AJ+AM,BJ+BM,CJ+CM)+2
    DO 2 I=1,ZZ
      Q(I,1)=1.d0
      Q(I,I)=1.d0
    2	CONTINUE
    DO 3 I=2,ZZ-1
    DO 3 K=2,I
      Q(I+1,K)=Q(I,K-1)+Q(I,K)
    3	CONTINUE
    CG=0.d0
    JA=AJ+AM+1.01d0
    MA=AJ-AM+1.01d0
    JB=BJ+BM+1.01d0
    MB=BJ-BM+1.01d0
    JC=CJ+CM+1.01d0
    MC=CJ-CM+1.01d0
    LA=BJ+CJ-AJ+1.01d0
    LB=CJ+AJ-BJ+1.01d0
    LC=AJ+BJ-CJ+1.01d0
    LT=AJ+BJ+CJ+1.01d0
    D=DABS(AM+BM-CM)-0.01d0
    IF (D) 10,10,20
      10	LD=MIN0(JA,JB,JC,MA,MB,MC,LA,LB,LC)
    IF (LD) 20,20,30
      30	JA2=AJ+AJ+AM+AM
      JB2=BJ+BJ+BM+BM
      JC2=CJ+CJ-CM-CM
      I2=JA2+JB2+JC2-JA2/2*2-JB2/2*2-JC2/2*2
    IF (I2) 20,40,20
      40	FN=Q(JA+MA-1,LC)/Q(LT,JC+MC-1)
    FN=FN*Q(JB+MB-1,LC)/Q(LT+1,2)
    FN=FN/Q(JA+MA-1,JA)
    FN=FN/Q(JB+MB-1,JB)
    FN=FN/Q(JC+MC-1,JC)
    K0=MAX(0,LC-JA,LC-MB)+1
    K1=MIN(LC,MA,JB)
    X=0.d0
    DO 50 K=K0,K1
      X=-X-Q(LC,K)*Q(LB,MA-K+1)*Q(LA,JB-K+1)
    50	CONTINUE
    IP=K1+LB+JC
    P=1-2*(IP-IP/2*2)
    CG=P*X*DSQRT(FN)
    CG=CG*DSQRT(2*CJ+1)*(-1)**IDNINT(AJ-BJ-CM)
    20	CONTINUE
    RETURN
  END SUBROUTINE CLEBSCH

  SUBROUTINE wigner(J1_val,J2_val,K1_val,K2_val,lamd_val,l1,l2, &
  reduced_matriz_element)
    DOUBLE PRECISION CG(4)
    DOUBLE PRECISION reduced_matriz_element
    REAL J1_val,J2_val,lamd_val
    INTEGER K1_val,K2_val,l1,l2
    INTEGER delta, factor,mu
    DOUBLE PRECISION f1
    INTEGER f2
    IF (K2_val==0) THEN
      delta=1
    ELSE
      delta=0
    ENDIF
    IF (K1_val==0) THEN
      delta1=1
    ELSE
      delta1=0
    ENDIF
    IF (K2_val==2 .OR. K1_val==2) THEN
      mu=2; factor=1
    ELSE
      mu=0; factor=2
    ENDIF
    CALL CLEBSCH(DBLE(J2_val),DBLE(lamd_val),DBLE(J1_val),DBLE(K2_val), &
    DBLE(mu),DBLE(K1_val),CG(1))
    CALL CLEBSCH(DBLE(J2_val),DBLE(lamd_val),DBLE(J1_val),-1d0*DBLE(K2_val), &
    DBLE(mu),DBLE(K1_val),CG(2))
    CALL CLEBSCH(DBLE(J2_val),DBLE(lamd_val),DBLE(J1_val),DBLE(K2_val), &
    DBLE(mu),-1d0*DBLE(K1_val),CG(3))
    CALL CLEBSCH(DBLE(J2_val),DBLE(lamd_val),DBLE(J1_val),-1d0*DBLE(K2_val), &
    DBLE(mu),-1d0*DBLE(K1_val),CG(4))
    f1= DSQRT(2*DBLE(J2_val)+1d0)/(DSQRT((1d0+DBLE(delta1))*(1d0+DBLE(delta))))
    f2= (1+(-1)**(lamd_val+K2_val+l1+l2))/2
    reduced_matriz_element=f1*f2*(CG(1) + ((-1)**(J2_val+l2))*CG(2) + &
    ((-1)**(J1_val+l1))*CG(3) + ((-1)**(J1_val+J2_val+l1+l2))*CG(4))
    reduced_matriz_element=reduced_matriz_element/factor
  END SUBROUTINE wigner

  !     *******************************************************
  !     START of dispersive PACK
  !     *******************************************************
  !==========================================================================
  !     AUTHOR: Dr. Roberto Capote Noy
  !
  !     e-mail: r.capotenoy@iaea.org ; rcapotenoy@yahoo.com;
  !
  !     DISPERSIVE OPTICAL MODEL POTENTIAL PACKAGE
  !
  !     Analytical dispersive integrals are included
  !     see Quesada JM, Capote R et al,
  !             Computer Physics Communications 153(2003) 97
  !             Phys. Rev. C67(2003) 067601
  !
  !     Dispersive integral's derivatives calculated by Dr.J.M.Quesada
  !
  REAL FUNCTION DOM_INT_Wv (Ef,Ep,Av,Bv,n,Einc,DerivIntWv)

    IMPLICIT NONE
    REAL Ef,Ep,Av,Bv,E,pi,Einc
    REAL E0,Ex,Eplus,Emin,Rs,ResEmin,ResEplus
    REAL DerEmin, DerEplus, Rds, DerivIntWv
    DOUBLE COMPLEX Pj,I,Zj,Ztmp
    DOUBLE COMPLEX Fs,Ds
    INTEGER N,j,IS

    DATA I/(0.d0,1.d0)/

    pi=4.d0*atan(1.d0)

    IS = 1
    E = Einc
    IF(Einc.LE.Ef) THEN
      E=2.d0*Ef-Einc
      IS = -1
    ENDIF

    E0 = Ep - Ef
    Ex = E  - Ef
    Eplus = Ex + E0
    Emin  = Ex - E0
    DOM_INT_Wv = 0.d0
    DerivIntWv = 0.d0

    ResEmin  =  Emin**n / (Emin**n + Bv**n)

    DerEmin  =  Emin**(n-1)* &
                ( Emin**n + Bv**n*(1.d0 + n*log(ABS(Emin)) ) ) &
                / (Emin**n + Bv**n)**2

    ResEplus = -Eplus**n / (Eplus**n + Bv**n)

    DerEplus = -Eplus**(n-1) * &
                   ( Eplus**n + Bv**n*(1.d0+n*log(Eplus)) ) &
                   / (Eplus**n + Bv**n)**2

    Fs = (0.d0,0.d0)
    Ds = (0.d0,0.d0)
    DO j=1,n
      Ztmp = I*(2*j-1)/dble(n)*pi
      Pj = Bv*exp(Ztmp)
      Zj = Pj * (2*Pj +Eplus -Emin) * Ex
      Zj = Zj / ( (Pj+E0) * (Pj+Eplus) * (Pj-Emin) )
      Fs = Fs + Zj*log(-Pj)
      Ds = Ds + 2*Pj*(Ex*Ex + (Pj+E0)**2)*log(-Pj) &
                  /( (Pj+Eplus)**2 * (Pj-Emin)**2 )
    ENDDO

    IF(ABS(IMAG(Fs)).gt.1.d-4) STOP 'Too big imag part in Wv'
    IF(ABS(IMAG(Ds)).gt.1.d-4) STOP 'Too big imag deriv in Wv'
    Rs  = REAL(Fs)
    Rds = REAL(Ds)
    DOM_INT_Wv = -Av/pi*IS* &
        (Rs/n  + ResEplus*log(Eplus) + ResEmin*log(ABS(Emin)))
    DerivIntWv =  Av/pi*IS*( Rds/n + DerEplus + DerEmin)

    RETURN
  END FUNCTION DOM_INT_Wv

  REAL FUNCTION DOM_INT_Ws (Ef,Ep,As,Bs,Cs,m,Einc,DerivIntWs)

    IMPLICIT NONE
    REAL Ef,Ep,As,Bs,Cs,E,Einc
    DOUBLE COMPLEX I,Pj,Zj,Ztmp
    REAL E0,Ex,pi
    REAL Rs,ResEmin,ResEplus
    REAL DerivIntWs,DerEmin,DerEplus,Rds
    INTEGER m,j,IS
    DOUBLE COMPLEX Fs,Ds
    REAL*8 Emin,Eplus

    DATA I/(0.d0,1.d0)/

    pi=4.d0*atan(1.d0)

    IS = 1
    E = Einc
    IF(Einc.LE.Ef) THEN
      E=2.d0*Ef-Einc
      IS = -1
    ENDIF

    E0 = Ep - Ef
    Ex = E  - Ef
    Eplus = Ex + E0
    Emin  = Ex - E0
    DOM_INT_Ws = 0.d0
    DerivIntWs = 0.d0
    ResEmin  =  Emin**m / (Emin**m + Bs**m)

    DerEmin  = -Emin**(m-1) * &
                 ( Emin**m + Bs**m + ( -Cs*Emin**(m+1) + &
                Bs**m *(-Cs*Emin+m) ) * exp(-Cs*Emin)*EIn(Cs*Emin) ) &
                / (Emin**m + Bs**m)**2

    ResEplus = -Eplus**m / (Eplus**m + Bs**m)

    DerEplus =  Eplus**(m-1) * &
                  ( Eplus**m + Bs**m + ( Cs*Eplus**(m+1) + &
                   Bs**m *(Cs*Eplus+m) ) * exp(Cs*Eplus)*EIn(-Cs*Eplus) ) &
                   / (Eplus**m + Bs**m)**2

    Fs = (0.d0,0.d0)
    Ds = (0.d0,0.d0)
    DO j=1,m
         Ztmp = I*(2*j-1)/dble(m)*pi
         Pj = Bs*exp(Ztmp)
         Zj = Pj * (2*Pj +Eplus -Emin) * Ex
         Zj = Zj / (Pj+E0) / (Pj+Eplus) / (Pj-Emin)
         Fs = Fs + Zj* zfi(-Pj*Cs)
         Ds = Ds + 2*Pj*(Ex*Ex + (Pj+E0)**2)*zfi(-Pj*Cs) &
                  /( (Pj+Eplus)**2 * (Pj-Emin)**2 )
    ENDDO

    IF(ABS(IMAG(Fs)).GT.1.d-4) STOP 'Too big imag part in Ws'
    IF(ABS(IMAG(Ds)).GT.1.d-4) STOP 'Too big imag deriv in Ws'
    Rs = REAL(Fs)
    Rds = REAL(Ds)

    DOM_INT_Ws = As/pi*IS*(Rs/m &
                         - ResEplus*exp(Cs*Eplus)*EIn(-Cs*Eplus) &
                         - ResEmin*exp(-Cs*Emin)*EIn(Cs*Emin) )
    RETURN
  END FUNCTION DOM_INT_Ws

  REAL FUNCTION WV(A,B,Ep,Ef,E,n)

    IMPLICIT NONE
    REAL  A,B,Ep,Ef,E,ee
    INTEGER n
    WV=0.d0
    IF(E.LE.Ef) E=2.d0*Ef-E
    IF(E.LT.Ep) RETURN
    ee=(E-Ep)**n
    WV=A*ee/(ee+B**n)
    RETURN
  END FUNCTION WV

  REAL FUNCTION WDD(A,B,C,Ep,Ef,E,m)
    IMPLICIT NONE
    REAL A,B,C,Ep,Ef,E,ee,arg
    INTEGER m
    WDD=0.d0
    IF(E.LE.Ef) E=2.d0*Ef-E
    IF(E.LT.Ep) RETURN
    arg=C*(E-Ep)
    IF(arg.GT.15) RETURN
    ee=(E-Ep)**m
    WDD=A*ee/(ee+B**m)*EXP(-arg)
    RETURN
  END FUNCTION WDD

  REAL FUNCTION DOM_int_T1(Ef,Ea,E)

    IMPLICIT NONE
    REAL E,Ea,Ef,Ex,Ea2,Eax,Pi,T11,T12,T13
    Pi=4.d0*ATAN(1.d0)
    Ex=E-Ef
    Ea2=Ea**2
    Eax=Ex+Ea
    T11 = 0.5d0*log(Ea)/Ex
    T12 =  ( (2*Ea+Ex)*log(Ea)+0.5d0*pi*Ex ) &
             /(2.*(Eax**2 + Ea2))
    T13 = -Eax**2*log(Eax)/(Ex*(Eax**2+Ea2))
    DOM_int_T1 = Ex/Pi*(T11+T12+T13)
    RETURN
  END FUNCTION DOM_int_T1

  REAL FUNCTION DOM_int_T2(Ef,Ea,E)

    IMPLICIT NONE
    REAL E,Ea,Ef,EL,Pi
    Pi=4.d0*ATAN(1.d0)
    EL=Ef+Ea
    DOM_int_T2= 1.d0 / Pi * ( &
             sqrt(abs(Ef)) * atan( (2*sqrt(EL*abs(Ef)))/(EL-abs(Ef)) ) &
        +    EL**1.5d0/(2*Ef)*log(Ea/EL) )
    IF(E.GT.EL) THEN
      DOM_int_T2 = DOM_int_T2 + 1.d0/Pi* ( &
         sqrt(E) * log( (sqrt(E)+sqrt(EL)) / (sqrt(E)-sqrt(EL)) ) + &
         1.5d0*sqrt(EL)*log((E-EL)/Ea) + EL**1.5d0/(2*E)*log(EL/(E-EL)) )
    ELSEIF(E.EQ.EL) THEN
      DOM_int_T2 = DOM_int_T2 + 1.d0/Pi*1.5d0*sqrt(EL) &
        *log((2**(4.d0/3.d0)*EL)/Ea)
    ELSEIF(E.GT.0.d0 .AND. E.LE.EL) THEN
      DOM_int_T2 = DOM_int_T2 + 1.d0/Pi * ( &
        sqrt(e) * log( (sqrt(E)+sqrt(EL)) / (sqrt(EL)-sqrt(E)) ) + &
        1.5d0*sqrt(EL)*log((EL-E)/Ea)+EL**1.5d0/(2.d0*E)*log(EL/(EL-E)) )
    ELSEIF(abs(E)<1e-10) then
      DOM_int_T2 = DOM_int_T2 + 1.d0/Pi*( 1.5*sqrt(EL) &
        * log(EL/Ea) + 0.5d0*sqrt(EL) )
    ELSE
      DOM_int_T2 = DOM_int_T2 + 1.d0/Pi * ( &
        -sqrt(abs(E))*atan( 2*(sqrt(EL*abs(E))) / (EL-abs(E)) ) + &
        1.5d0*sqrt(EL)*log((EL-E)/Ea)+EL**1.5d0/(2.d0*E)*log(EL/(EL-E)) )
    ENDIF
  	!write(101,*) E,DOM_int_T2,Ef,Ea
    RETURN
  END FUNCTION DOM_int_T2

  DOUBLE COMPLEX FUNCTION zfi(za)

    IMPLICIT NONE
    REAL aj
    DOUBLE COMPLEX za,y
    INTEGER m,i
    zfi=0.d0
    IF (za.EQ.0.d0) RETURN
    IF (abs(real(za)+18.5d0).GE.25.d0) GO TO 3
    IF (SQRT(625.d0-(REAL(za)+18.5d0)**2)/1.665d0.LT.ABS(imag(za))) GO TO 3
    zfi=-.57721566490153d0-log(za)
    y=1.d0
    DO 1 m=1,2000
      aj=m
      y=-y*za/aj
      IF (ABS(y).lt.1.d-15*ABS(zfi)) GO TO 2
      1 zfi=zfi-y/aj
      2 zfi=EXP(za)*zfi
        RETURN
      3 DO 4 i=1,20
          aj=21-i
          zfi=aj/(za+zfi)
          4 zfi=aj/(1.d0+zfi)
          zfi=1.d0/(zfi+za)
    RETURN
  END FUNCTION zfi

  REAL*8 FUNCTION EIn(X)

    IMPLICIT NONE
    REAL*8 FAC, H, X
    INTEGER N
    EIn = 0.57721566490153d0+LOG(ABS(X))
    FAC = 1.0
    DO N = 1,100
      H = FLOAT(N)
      FAC = FAC*H
      EIn = EIn + X**N/(H*FAC)
    ENDDO
    RETURN
  END FUNCTION EIn

  SUBROUTINE dispers2(A,Z,k,eopt, &
  v,rvv,avv, dv,drv,dav, dvs,drs,das, w,rw,aw, wd,rwd,awd, &
  vso,rvso,avso,dvso, wso,rwso,awso, &
  Vlin,Vdep,lambdaHF,Cviso,Vso0,lambdaso,Ccoul, &
  AAv,BBv,W0,BBs,CCs,Cwiso,Wso0,BBso, &
  Ea,alpha,eferm,Ades, &
  rHFl,rHFdep,aHFl,aHFdep,rv,avl,avdep, &
  rsl,rsdep,as, &
  rso,aso,rc,ac)

    REAL eopt,asym,eferm,f,Cviso,viso,Ccoul,Cwiso
    REAL lambdaHF,lambdaso,Ades
    pi = 4.0*atan(1.)
    !
    ! *** Parameters of Soukhovitskii, Capote, Quesada, Chiba and Martyanov (Nov 25, 2015) ***
    ! *** with Asymmetrical W energy-dependence
    ! *** dispers2:  T1 integral with correct coefficient
    ! *** done by Ian Thompson
    ! k         : designator for particle
    ! Z         : charge number of residual nucleus
    ! A         : mass number of residual nucleus
    ! eopt      : incident energy
    ! asym      : asymmetry parameter
    ! eferm     : Fermi energy
    ! f         : eopt-eferm

    Au = A-Ades
    asym=(A-2.*Z)/A
    V0 = Vlin + Vdep*Au
    rHF = rHFL + rHFdep * Au
    aHF = aHFl + aHFdep * Au
    av = avl + avdep * Au
    rs = rsl + rsdep * Au
    AAHF = V0 * (1 + (-1)**k * Cviso*asym/V0)
    AAs  = W0 * (1 + (-1)**k * Cwiso*asym/W0)
    eoffset = 0.
    IF (k==2) eoffset =  Ccoul * Z/A**(1./3.)
    Eeff = eopt - eoffset
    f = Eeff  - eferm
    v = AAHF * EXP(-lambdaHF*f)
    vso=Vso0*EXP(-lambdaso*f)

    ! sources of dispersive terms
    w = AAv * f*f/(f*f + BBv**2)
    IF(f < -Ea) THEN
      fe  = f + Ea
      w = w * (1 - fe*fe/(fe*fe + Ea*Ea))
    ELSE IF(f>Ea) THEN
      w = w + alpha * (SQRT(Eeff) + (eferm+Ea)**1.5d0/(2*Eeff) &
           	- 1.5d0*SQRT(eferm+Ea))
    endif
    wd = AAs * f*f/(f*f + BBs**2) * EXP( -CCs * ABS(f))
    wso=Wso0* f*f/(f*f + BBso**2)
    ! dispersive terms to add for real volume and real surface forms
    drv= rv;  dav = av
    drs = rs; das = as
    dvs = DOM_INT_Ws (eferm,eferm,AAs,BBs,CCs,2,Eeff,DerivIntWs)
    DWv = DOM_INT_Wv (eferm,eferm,AAv,BBv,2,Eeff,DerivIntWv)
    T1 = DOM_int_T1(eferm,Ea,Eeff) * AAv * Ea*Ea/(Ea*Ea + BBv**2)
    T2 = DOM_int_T2(eferm,Ea,Eeff) * alpha
    dv = DWv + T1 + T2
    dvso = DOM_INT_Wv (eferm,eferm,Wso0,BBso,2,Eeff,DerivIntWv)


    ! name translations
    rvv = rHF ; avv = aHF
    rvso = rso; avso = aso
    rwso = rso; awso = aso
    rw  = rv  ; aw   = av
    rwd = rs  ; awd  = as
    RETURN
  END SUBROUTINE dispers2
  ! ***********************************************************
  ! *                    END of dispersive                    *
  ! ***********************************************************

!----------------------------------------------------------------------------

  ! **************************************************************************
  ! *                                                                        *
  ! *  Numerical calculation of form factors using Gauss-Legendre quadrature *
  ! *                                                                        *
  ! **************************************************************************

  ! - Form factors calculated for volume and surface terms with axial deformations up to beta_{60}.
  ! - Additional 1/sqrt(4*pi) factor added to match FRESCO's convention.
  ! - Output (extension .form) prepared to be used as FRESCO's external forms input [fort.4].

  DOUBLE PRECISION FUNCTION arm(x,L)

    DOUBLE PRECISION x,P(0:L),pi
    INTEGER L,i
    pi=ACOS(-1d0)
    IF(L==0) THEN
      arm=DSQRT(1.d0/(4d0*pi))
    ENDIF
    IF(L==1) THEN
      arm=DSQRT(3.d0/(4d0*pi))*x
    ENDIF
    IF(L.GE.2) THEN
      P(0)=1d0
      P(1)=x
      DO i=2,L,1
        P(i)=((2d0*DBLE(i)-1d0)*x*P(i-1) - DBLE(i-1)*P(i-2))/(DBLE(i))
      END DO
      arm=dsqrt((2d0*DBLE(L)+1d0)/(4d0*pi))*P(L)
    ENDIF
    RETURN
  END FUNCTION arm

  DOUBLE PRECISION FUNCTION rdef(R,x,BETA2,BETA4,BETA6)

    REAL R,BETA2,BETA4,BETA6
    DOUBLE PRECISION x
    rdef=DBLE(R)*(1d0+DBLE(BETA2)*arm(x,2)+DBLE(BETA4)*arm(x,4)+DBLE(BETA6)*arm(x,6))
    RETURN
  END FUNCTION rdef

  DOUBLE PRECISION FUNCTION vws(r,v0,RR,a,x,L,BETA2,BETA4,BETA6)

    REAL r,v0,RR,a,BETA2,BETA4,BETA6
    DOUBLE PRECISION x
    INTEGER L
    vws=(DBLE(v0)/(1d0+EXP((DBLE(r)-rdef(RR,x,BETA2,BETA4,BETA6))/DBLE(a))))*arm(x,L)
    RETURN
  END FUNCTION vws

  DOUBLE PRECISION FUNCTION vs(r,v0,RR,a,x,L,BETA2,BETA4,BETA6)

    REAL r,v0,RR,a,BETA2,BETA4,BETA6
    DOUBLE PRECISION x,f
    INTEGER L
    f=EXP((DBLE(r)-rdef(RR,x,BETA2,BETA4,BETA6))/DBLE(a))
    vs=-4d0*DBLE(v0)*arm(x,L)*f/(1+f)**2
    RETURN
  END FUNCTION vs

  DOUBLE PRECISION FUNCTION gaussv(r,v0,RR,a,L,BETA2,BETA4,BETA6)

    REAL r,v0,RR,a,x(19),w(19),BETA2,BETA4,BETA6
    DOUBLE PRECISION pp,pi,p
    INTEGER j,L
    pi=ACOS(-1d0)
    x(1:10)=(/-.992407,-.960208,-.903156,-.822715,-.720966, &
        -.600545,-.464571,-.316564,-.160359,0.0/)
    w(1:10)=(/.0194618,.0448182,.0690445,.09149,.111567,.128754, &
                .142607,.152766,.158969,.161054/)
    DO i=1,9,1
      x(10+i)=-x(10-i)
      w(10+i)=w(10-i)
    END DO
    p=0d0
    DO j=1,19,1
      p=p+DBLE(w(j))*vws(r,v0,RR,a,DBLE(x(j)),L,BETA2,BETA4,BETA6)
    END DO
    gaussv=DSQRT(pi)*p !Factor required by FRESCO's convention.
    RETURN
  END FUNCTION gaussv

  DOUBLE PRECISION FUNCTION gausss(r,v0,RR,a,L,BETA2,BETA4,BETA6)

    REAL r,v0,RR,a,x(19),w(19),BETA2,BETA4,BETA6
    DOUBLE PRECISION pp,pi,p
    INTEGER j,L
    pi=ACOS(-1d0)
    x(1:10)=(/-.992407,-.960208,-.903156,-.822715,-.720966, &
      -.600545,-.464571,-.316564,-.160359,0.0/)
    w(1:10)=(/.0194618,.0448182,.0690445,.09149,.111567,.128754, &
              .142607,.152766,.158969,.161054/)
    DO i=1,9,1
      x(10+i)=-x(10-i)
      w(10+i)=w(10-i)
    END DO
    p=0d0
    DO j=1,19,1
      p=p+DBLE(w(j))*vs(r,v0,RR,a,DBLE(x(j)),L,BETA2,BETA4,BETA6)
    END DO
    gausss=-DSQRT(pi)*p !Factors required by FRESCO's convention
    RETURN
  END FUNCTION gausss

  SUBROUTINE FORMFACT(VR,RR,AR,dv,drv,dav,W,RW,AW, &
  VD,RVD,AVD,WD,RD,AD,N,rmax,BETA2,BETA4,BETA6,A,nexo,sump)

    DOUBLE PRECISION deltar,r(0:N)
    REAL rmax,BETA2,BETA4,BETA6,A
    REAL VR,RR,AR !Real volume
    REAL dv,drv,dav,W,RW,AW !Real/imaginary dispersive volume
    REAL VD,RVD,AVD,WD,RD,AD!Real/imaginary dispersive surface
    INTEGER N,i,nexo,sump
    deltar=DBLE(rmax)/(DBLE(N)-1d0)
    r(0)=0.0d0
    DO i=1,N-1,1
      r(i)=r(i-1)+deltar
    END DO
    IF(sump .NE. 0) THEN
      WRITE(94,*) '!Real volume'
      WRITE(94,20) N,deltar,R(0)
      DO i=0,N-1,1
        WRITE(94,10) gaussv(REAL(r(i)),-VR,RR*(A**(1./3.)),AR,2,BETA2,BETA4,BETA6)
      END DO
    ENDIF
    IF(nexo .NE. 0) THEN
      WRITE(94,*) '!Real volume'
      WRITE(94,20) N,deltar,R(0)
      DO i=0,N-1,1
        WRITE(94,10) gaussv(REAL(r(i)),-VR,RR*(A**(1./3.)),AR,2,BETA2,BETA4,BETA6)
      END DO
    ENDIF
    IF(sump .NE. 0) THEN
      WRITE(94,*) '!Dispersive real and imaginary volume'
      WRITE(94,20) N,deltar,R(0)
      DO i=0,N-1,1
        WRITE(94,10) gaussv(REAL(r(i)),-dv,drv*(A**(1./3.)),dav,2,BETA2,BETA4,BETA6)
        WRITE(94,10) gaussv(REAL(r(i)),-W,RW*(A**(1./3.)),AW,2,BETA2,BETA4,BETA6)
      END DO
    ENDIF
    IF(nexo .NE. 0) THEN
      WRITE(94,*) '!Dispersive real and imaginary volume'
      WRITE(94,20) N,deltar,R(0)
      DO i=0,N-1,1
        WRITE(94,10) gaussv(REAL(r(i)),-dv,drv*(A**(1./3.)),dav,2,BETA2,BETA4,BETA6)
        WRITE(94,10) gaussv(REAL(r(i)),-W,RW*(A**(1./3.)),AW,2,BETA2,BETA4,BETA6)
      END DO
    ENDIF
    IF(sump .NE. 0) THEN
      WRITE(94,*) '!Dispersive real and imaginary surface'
      WRITE(94,20) N,deltar,R(0)
      DO i=0,N-1,1
        WRITE(94,10) gausss(REAL(r(i)),-VD,RVD*(A**(1./3.)),AVD,2,BETA2,BETA4,BETA6)
        WRITE(94,10) gausss(REAL(r(i)),-WD,RD*(A**(1./3.)),AD,2,BETA2,BETA4,BETA6)
      END DO
    ENDIF
    IF(nexo .NE. 0) THEN
      WRITE(94,*) '!Dispersive real and imaginary surface'
      WRITE(94,20) N,deltar,R(0)
      DO i=0,N-1,1
        WRITE(94,10) gausss(REAL(r(i)),-VD,RVD*(A**(1./3.)),AVD,2,BETA2,BETA4,BETA6)
        WRITE(94,10) gausss(REAL(r(i)),-WD,RD*(A**(1./3.)),AD,2,BETA2,BETA4,BETA6)
      END DO
    ENDIF
    20 FORMAT(' ',i2,' ',f9.6,' ',f9.6)
    10 FORMAT(' ',f9.6)
    CLOSE(94)
    RETURN
  END SUBROUTINE FORMFACT

  ! **************************************************************************
  ! *                                                                        *
  ! *                        End of form factors pack                        *
  ! *                                                                        *
  ! **************************************************************************

!-------------------------------------------------------------------------------

  ! ****************************************************************************
  ! *                                                                          *
  ! *      Interband coupling between ground state band and excited bands for  *
  ! *                               even-even nuclei                           *
  ! *                                                                          *
  ! ****************************************************************************

  ! - Subroutine used to build  &STEPS part of FRESCO's input.

  SUBROUTINE steps(J1,J2,p1,p2,I1,I2,K1,K2,BETA_VAL)
    REAL J1,J2,lamd,BETA_VAL
    INTEGER I1,I2,K1,K2,lamd_ph_1,lamd_ph_2
    INTEGER p1,p2
    DOUBLE PRECISION wigner_val
    REAL jmin,jmax,jaux(40)
    lamd_ph_1=0 ! No phonons in the G.S band.
    IF (p2 .LT. 0) THEN
      lamd=3.0
      lamd_ph_2=3
    ELSE
      lamd=2.0
      lamd_ph_2=2
    ENDIF
     2 format('&STEP ia=',i2,' ib=',i2,' k=',i1,' str=',1f9.4,'/')
    IF (J1.GE.lamd) THEN
      jmin=J1-lamd
      jmax=J1+lamd
      DO i=0,INT(2*lamd),1
        jaux(i+1)=jmin+i
      ENDDO
      siz=2*lamd+1
    ELSE
      jmin=lamd-J1
      jmax=lamd+J1
      DO i=0,2*INT(J1),1
        jaux(i+1)=jmin+i
      ENDDO
      siz=2*J1+1
    ENDIF
    DO i=1,INT(siz)
      IF(jaux(i) .EQ. J2  .AND. J2 .LE. jmax ) THEN
        CALL wigner(J1,J2,K1,K2,lamd,lamd_ph_1,lamd_ph_2,wigner_val)
        WRITE(1,2) I1,I2,INT(lamd),((-1)**((J2-J1+ABS(J2-J1))/2))*wigner_val*BETA_VAL ! Phase accoding to FRESCO's convention.
        CALL wigner(J2,J1,K2,K1,lamd,lamd_ph_2,lamd_ph_1,wigner_val)
        WRITE(1,2) I2,I1,INT(lamd),((-1)**((J1-J2+ABS(J1-J2))/2))*wigner_val*BETA_VAL
      ENDIF
    ENDDO
    RETURN
  END SUBROUTINE steps

  ! - Stop in case of allocation error.
  SUBROUTINE error (errval,errtype)
    INTEGER errval,errtype
    IF(errval.NE.0.AND.errtype.EQ.1) THEN
      STOP "Allocation request denied"
    ENDIF
    IF(errval.NE.0.AND.errtype.EQ.2) THEN
      STOP "Deallocation request denied"
    ENDIF
    RETURN
  END SUBROUTINE error

  ! - Generate xmgrace's input in case of drawing cross section graphs at the end of the run.
  SUBROUTINE GRACE
    OPEN(94,form='formatted',file='graphs.gr')
    WRITE(94,*) '#Generate .eps archives with cross sections graphs.'
    WRITE(94,*) 'READ BLOCK "xsec_completo.out"'
    WRITE(94,*) 'FOCUS G0'
    WRITE(94,*) 'BLOCK XY "1:2"'
    WRITE(94,*) 's0 symbol 1'
    WRITE(94,*) 's0 symbol color 2'
    WRITE(94,*) 's0 legend "Elastic cross section"'
    WRITE(94,*) 'world ymax MAX(s0.y)+0.05*MAX(s0.y)'
    WRITE(94,*) 'world ymin MIN(s0.y)-0.05*MIN(s0.y)'
    WRITE(94,*) 'legend 0.7, 0.65'
    WRITE(94,*) 'title "PRC2016 Elastic C.S"'
    WRITE(94,*) 'xaxis label "Energy (MeV)"'
    WRITE(94,*) 'yaxis label "Cross section (mb)"'
    WRITE(94,*) 'VIEW 0.15, 0.15, 1.15, 0.85'
    WRITE(94,*) 'PRINT TO "elastic.eps"'
    WRITE(94,*) 'PRINT'
    WRITE(94,*) '##################################################'
    WRITE(94,*) 'KILL G0.s0'
    WRITE(94,*) 'FOCUS G1'
    WRITE(94,*) 'BLOCK XY "1:3"'
    WRITE(94,*) 's0 symbol 4'
    WRITE(94,*) 's0 symbol color 3'
    WRITE(94,*) 's0 legend "Absorption cross section"'
    WRITE(94,*) 'world ymax MAX(s0.y)+0.05*MAX(s0.y)'
    WRITE(94,*) 'world ymin MIN(s0.y)-0.05*MIN(s0.y)'
    WRITE(94,*) 'title "PRC2016 Absorption C.S"'
    WRITE(94,*) 'xaxis label "Energy (MeV)"'
    WRITE(94,*) 'yaxis label "Cross section (mb)"'
    WRITE(94,*) 'PRINT TO "absorption.eps"'
    WRITE(94,*) 'PRINT'
    WRITE(94,*) '###################################################'
    WRITE(94,*) 'KILL G1.s0'
    WRITE(94,*) 'KILL G0.s0'
    WRITE(94,*) 'FOCUS G2'
    WRITE(94,*) 'BLOCK XY "1:5"'
    WRITE(94,*) 's0 symbol 1'
    WRITE(94,*) 's0 symbol color 4'
    WRITE(94,*) 's0 legend "Total cross section"'
    WRITE(94,*) 'world ymax MAX(s0.y)+0.05*MAX(s0.y)'
    WRITE(94,*) 'world ymin MIN(s0.y)-0.05*MIN(s0.y)'
    WRITE(94,*) 'title "PRC2016 Total C.S"'
    WRITE(94,*) 'xaxis label "Energy (MeV)"'
    WRITE(94,*) 'yaxis label "Cross section (mb)"'
    WRITE(94,*) 'PRINT TO "total.eps"'
    WRITE(94,*) 'PRINT'
    CLOSE(94)
  END SUBROUTINE GRACE
END MODULE modulo
