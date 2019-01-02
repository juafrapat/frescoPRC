MODULE modulo
  !***********************************************************************************************
  !* This module contains all functions and subroutines for frescoPRC program:                   *
  !*                                                                                             *
  !* Subroutine CLEBSCH: Compute C-G coefficients with a couple extra factors                    *
  !*                     DSQRT(2*AJ+1) and ((-1)**((AJ-CJ+dabs(AJ-CJ))/2))                       *
  !*                     required by FRESCO's convention.                                        *
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
    CG=CG*((-1)**((AJ-CJ+dabs(AJ-CJ))/2))
    CG=CG*DSQRT(2*AJ+1)
    20	CONTINUE
    RETURN
  END SUBROUTINE CLEBSCH

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

  SUBROUTINE steps (mu,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
  nval,neg,neo,neb,nega,neax,jl)

    INTEGER neg,np,neo,neb,nega,neax,n,nval
    DOUBLE PRECISION CG
    REAL jl(nval)
    INTEGER J(neg),JO(neo),JB(neb),JG(nega),JAX(neax)
    REAL BETA3,BETA2,BETA2EFF,GAMMA2EFF,GAMMANAX
    CHARACTER*10 banda
    INTEGER i,k,s,f,l,p,o,siz,mu
    INTEGER IE(6),IEO(5),IEb(2),IEg(3),IEax(2),LAMDA(4),JMIN,JMAX,JAUX(40)
    2   format('&STEP ia=',i2,' ib=',i2,' k=',i1,' str=',1f9.4,'/')
    DO i=1,nval
      IF (i.le.neg) THEN
        J(i)=INT(jl(i))
      ELSE IF (i.gt.neg .and. i.le.neg+neo) THEN
        JO(i-neg)=INT(jl(i))
      ELSE IF (i.gt.neg+neo .and. i.le.neg+neo+neb) THEN
        JB(i-neg-neo)=INT(jl(i))
      ELSE IF (i.gt.neg+neo+neb .and. i.le.neg+neo+neb+nega) THEN
        JG(i-neg-neo-neb)=INT(jl(i))
      ELSE IF (i.gt.neg+neo+neb+nega .and. i.le.nval) THEN
        JAX(i-neg-neo-neb-nega)=INT(jl(i))
      ENDIF
    ENDDO
    LAMDA(1:4)=(/2,4,6,3/)
    np=size(LAMDA)
    ! G.S -> G.S (test case)-------------------------------------------------------
    IF (mu==1) THEN
      DO i=1,neg
        IE(i)=i
      END DO
      DO i=1,neg,1
        DO k=1,3,1
          IF (J(i).GE.LAMDA(k)) THEN
            JMIN=J(i)-LAMDA(k)
            JMAX=J(i)+ LAMDA(k)
            DO s=0,2*LAMDA(k)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*LAMDA(k))+1
          ELSE
            JMIN=LAMDA(k)-J(i)
            JMAX=J(i)+LAMDA(k)
            DO s=0,2*J(i)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*J(i))+1
          ENDIF
          DO l=1,siz,1
            DO p=1,neg,1
              IF (JAUX(l).EQ.J(p)  .AND. IE(i).LT.IE(p) .AND. J(p).LE.JMAX) THEN
                 CALL CLEBSCH(DBLE(J(i)),DBLE(LAMDA(k)),DBLE(J(p)),0d0,0d0,0d0,CG)
                 WRITE(1,2) IE(i),IE(p),LAMDA(k),CG
                 WRITE(1,2) IE(p),IE(i),LAMDA(k),CG
              ENDIF
            END DO
          END DO
        END DO
      END DO
    ENDIF
    ! G.S --> Octupole --------------------------------------------------------
    IF (mu==2) THEN
      k=4
      DO i=1,neo
        IEO(i)=i+neg
      END DO
      DO i=1,neg
        IE(i)=i
      END DO
      DO i=1,neg,1
        IF (J(i).GE.LAMDA(k)) THEN
          JMIN=J(i)-LAMDA(k)
          JMAX=J(i)+ LAMDA(k)
          DO s=0,2*LAMDA(k)
            JAUX(s+1)=JMIN+s
          END DO
          siz=(2*LAMDA(k))+1
        ELSE
          JMIN=LAMDA(k)-J(i)
          JMAX=J(i)+LAMDA(k)
          DO s=0,2*JO(i)
            JAUX(s+1)=JMIN+s
          END DO
          siz=(2*JO(i))+1
        ENDIF
        DO l=1,siz,1
          DO p=1,neo,1
            IF (JAUX(l).EQ.JO(p)  .AND. IE(i).LT.IEO(p) .AND. JO(p).LE.JMAX) THEN
                 CALL CLEBSCH(DBLE(J(i)),DBLE(LAMDA(k)),DBLE(JO(p)),0d0,0d0,0d0,CG)
                 WRITE(1,2) IE(i),IEO(p),LAMDA(k),CG*BETA3
                 WRITE(1,2) IEO(p),IE(i),LAMDA(k),CG*BETA3
            ENDIF
          END DO
        END DO
      END DO
    ENDIF
    ! G.S ---> Beta -----------------------------------------------------------
    IF (mu==3) THEN
      k=1
      DO i=1,neb
        IEb(i)=i+neo+neg
      END DO
      DO i=1,neg
        IE(i)=i
      END DO
      DO i=1,neg,1
        IF (J(i).GE.LAMDA(k)) THEN
          JMIN=J(i)-LAMDA(k)
          JMAX=J(i)+ LAMDA(k)
          DO s=0,2*LAMDA(k)
            JAUX(s+1)=JMIN+s
          END DO
          siz=(2*LAMDA(k))+1
        ELSE
          JMIN=LAMDA(k)-J(i)
          JMAX=J(i)+LAMDA(k)
          DO s=0,2*JB(i)
            JAUX(s+1)=JMIN+s
          END DO
          siz=(2*JB(i))+1
        ENDIF
        DO l=1,siz,1
          DO p=1,neb,1
            IF (JAUX(l).EQ.JB(p)  .AND. IE(i).LT.IEb(p) .AND. JB(p).LE.JMAX) THEN
              CALL CLEBSCH(DBLE(J(i)),DBLE(LAMDA(k)),DBLE(JB(p)),0d0,0d0,0d0,CG)
              WRITE(1,2) IE(i),IEb(p),LAMDA(k),CG*BETA2EFF
              WRITE(1,2) IEb(p),IE(i),LAMDA(k),CG*BETA2EFF
            ENDIF
          END DO
        END DO
      END DO
    ENDIF
    ! G.S ---> GAMMA ----------------------------------------------------------
    IF (mu==4) THEN
      k=1
      DO i=1,nega
        IEg(i)=i+neo+neg+neb
      END DO
      DO i=1,neg
        IE(i)=i
      END DO
      DO i=1,neg,1
        IF (J(i).GE.LAMDA(k)) THEN
          JMIN=J(i)-LAMDA(k)
          JMAX=J(i)+ LAMDA(k)
          DO s=0,2*LAMDA(k)
            JAUX(s+1)=JMIN+s
          END DO
          siz=(2*LAMDA(k))+1
        ELSE
          JMIN=LAMDA(k)-J(i)
          JMAX=J(i)+LAMDA(k)
          DO s=0,2*JG(i)
            JAUX(s+1)=JMIN+s
          END DO
          siz=(2*JG(i))+1
        ENDIF
        DO l=1,siz,1
          DO p=1,nega,1
            IF (JAUX(l).EQ.JG(p)  .AND. IE(i).LT.IEg(p) .AND. JG(p).LE.JMAX) THEN
              CALL CLEBSCH(DBLE(J(i)),DBLE(LAMDA(k)),DBLE(JG(p)),0d0,0d0,0d0,CG)
              WRITE(1,2) IE(i),IEg(p),LAMDA(k),CG*GAMMA2EFF
              WRITE(1,2) IEg(p),IE(i),LAMDA(k),CG*GAMMA2EFF
            ENDIF
          END DO
        END DO
      END DO
    ENDIF
    ! G.S ---> GAMMA NO AXIAL -------------------------------------------------
    IF (mu==5) THEN
      k=1
      DO i=1,neax
        IEax(i)=i+neo+neg+neb+nega
      END DO
      DO i=1,neg
        IE(i)=i
      END DO
      DO i=1,neg,1
        IF (J(i).GE.LAMDA(k)) THEN
          JMIN=J(i)-LAMDA(k)
          JMAX=J(i)+ LAMDA(k)
          DO s=0,2*LAMDA(k)
            JAUX(s+1)=JMIN+s
          END DO
          siz=(2*LAMDA(k))+1
        ELSE
          JMIN=LAMDA(k)-J(i)
          JMAX=J(i)+LAMDA(k)
          DO s=0,2*JAX(i)
            JAUX(s+1)=JMIN+s
          END DO
          siz=(2*JAX(i))+1
        ENDIF
        DO l=1,siz,1
          DO p=1,neax,1
            IF (JAUX(l).EQ.JAX(p)  .AND. IE(i).LT.IEax(p) .AND. JAX(p).LE.JMAX) THEN
              CALL CLEBSCH(DBLE(JAX(p)),DBLE(LAMDA(k)),DBLE(J(i)),-2d0,2d0,0d0,CG)
              WRITE(1,2) IE(i),IEax(p),LAMDA(k),sqrt(2.0)*((-1)**(DBLE(JAX(p))))*CG*GAMMANAX
              WRITE(1,2) IEax(p),IE(i),LAMDA(k),sqrt(2.0)*((-1)**(DBLE(JAX(p))))*CG*GAMMANAX
            ENDIF
          END DO
        END DO
      END DO
    ENDIF
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
END MODULE modulo
