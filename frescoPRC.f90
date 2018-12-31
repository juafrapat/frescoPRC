!   frescoPRC => Input generator to perform calculations based in the model presented in
!		PRC 94 6 (2016), 064605 for even-even actinides [effective couplings between
!		bands with dispersive corrections to the optical potential] with FRESCO.
  USE modulo
  PARAMETER (mxsym=100)
  CHARACTER*5 NAME
  CHARACTER*8 POTL
  DOUBLE PRECISION CG
  CHARACTER*100 fname,pname,h,potname
  CHARACTER*2 rela,SYMBOL(mxsym),PT(4)
  REAL elevels(20), jlevels(20),enodes(10),val(4,10,8), enalevels(20),jnalevels(20)
  REAL eolevels(20), jolevels(20)
  REAL eulevels(20), julevels(20)
  REAL eilevels(20), jilevels(20),lambdaHF,lambdaso
  INTEGER za,shape,kpnodes(10),fnodes,mu,i, sump
  CHARACTER*312 ELEMENTS,out
  CHARACTER tab
  CHARACTER*12 pottype(8)
  tab = char(9)

  data pottype / 'REAL_VOLUME', 'REAL_VOLUME', 'IMAG_VOLUME','REAL_SURFACE', &
     &  		'IMAG_SURFACE', 'REAL_SPINORB', 'IMAG_SPINORB', 'REAL_SPINORB' /
        ELEMENTS =						&
     & ' H HE LI BE  B  C  N  O  F NE NA MG AL SI  P  S CL AR  K CA ' // &
     & 'SC TI  V CR MN FE CO NI CU ZN GA GE AS SE BR KR RB SR  Y ZR ' // &
     & 'NB MO TC RU RH PD AG CD IN SN SB TE  I XE CS BA LA CE PR ND ' // &
     & 'PM SM EU GD TB DY HO ER TM YB LU HF TA  W RE OS IR PT AU HG ' // &
     & 'TL PB BI PO AT RN FR RA AC TH PA  U NP PU AM CM BK CF ES FM'
  READ (ELEMENTS,1021)(SYMBOL(i), i=1, mxsym)
  1021 FORMAT (300(A2,1X))
  OPEN(40,STATUS='old',FILE='input.INP')
  READ(40,*) Z,A,eferm
  READ(40,*) nexe,nexo,nexu,nexi,nna
  NAME = symbol(nint(Z))//'000'
  WRITE(6,*) 'Z,A,name =',nint(Z),nint(A),symbol(nint(Z))
  WRITE(NAME(3:5),'(i3.3)') nint(A)
  READ(40,*) jtmax,hcm0,rmatch,Ngrid
  absend=0.001
  elevels(1:10) =(/ 0., .04490, .14840, 0.3074, 0.5174, 0.7759, 1.0767, 1.4155, 1.7884, 2.1911/)  ! G.S band
  jlevels(1:10) =(/ 0., 2., 4., 6., 8., 10., 12., 14., 16., 18. /)
  eolevels(1:5) =(/ 0.6798, .7313, .8267, 0.9664, 1.1507  /) ! Octopule band - parity (k=0-)
  jolevels(1:5) =(/ 1., 3., 5., 7., 9. /)
  eulevels(1:2)=(/0.993, 1.0373/) !beta band
  julevels(1:2)=(/0., 2./)
  eilevels(1:3)=(/0.9270, .9663, 1.0564/) !gamma band
  jilevels(1:3)=(/0.,2.,4./)
  enalevels(1:2)=(/1.0603,1.1057/) !NAX band
  jnalevels(1:2)=(/2.0,3.0/)
  POTL= 'dispers2'
  escale=1.0
  nex = nexe + nexo + nexu + nexi + nna
  sump = nexu + nexi + nna
  sump2=sump+nexo
  READ(40,*) EMIN,EMAX,NE
  nodes = 7
  fnodes = 6
  enodes(1:7) = (/ 0.01, 1., 5., 10., 20., 50. , 200./)
  kpnodes(1:7) = (/ 10, 11, 15, 20, 30, 32, 34 /)
  kpp=1
  IF(NAME(1:1) ==' ') NAME(1:5)=NAME(2:5)//' '
  pname = POTL//'-'//trim(NAME)//'-parameters.txt'
  OPEN(10,form='formatted',file=trim(pname))
  OPEN(69,form='formatted',file='lista.txt' )
  WRITE(10,1) POTL,NAME
  1	FORMAT('####### OPTICAL PARAMETERS for ',A8,' ########'/'#'/'#     neutron on ',a5,/'#'/ &
        & '#Energy    V     rv    av     dV    drv   dav      W     rw    aw     ', &
        &           ' Vd   rvd   avd      Wd   rwd   awd     ', &
        &           'Vso   rvso  avso    dvso    Wso   rwso  awso  rc    ac')
  DE = 1.
  IF(NE>1) DE = (EMAX-EMIN)/(NE-1)
  NA = nint(A); NZ = nint(Z)
  ACroot = real(NA)**(1./3.)
  Ccoul=1.36
  rc=1.2894
  ac=0.547
                      !Dispersive parameters
  READ(40,*) Vlin,Vdep,lambdaHF,Cviso,Vso0,lambdaso,Ccoul
  READ(40,*) AAv,BBv,W0,BBs,CCs,Cwiso,Wso0,BBso
  READ(40,*) Ea,alpha,Ades
  READ(40,*) rHFl,rHFdep,aHFl,aHFdep,rv,avl,avdep
  READ(40,*) rsl,rsdep,as
  READ(40,*) rso,aso,rc,ac
            !Deformation parameters => Multiplied by BETA_{20} like in OPTMAN input.
  READ(40,*) BETA2,BETA4,BETA6
  READ(40,*) BETA2EFF,GAMMA2EFF,GAMMANAX,BETA3
  BETA3 = BETA3/(BETA2)!GS > OCTUPOLE BAND
  BETA2EFF=BETA2EFF/(BETA2)!GS > BETA BAND
  GAMMA2EFF=GAMMA2EFF/(BETA2)!GS > GAMMA BAND(mu=0)
  GAMMANAX=GAMMANAX/(BETA2) ! GS>NON-AXIAL GAMMA BAND (mu=2)
  DO IE=1,ABS(NE)
     IEN=IE-ABS(NE)
	   E=EMIN+(IE-1)*DE
	   hcm = hcm0
	   IF(E>50.0) hcm = hcm0/SQRT(E/50.0)
	   dv=0; drv=0; dav=0; dvso=0.
     NTYPE = 1 ! neutron!!
     CALL dispers2(A,Z,NTYPE,E,VR,RR,AR, dv,drv,dav, VD,RVD,AVD, &
                     W,RW,AW, WD,RD,AD, VSO,RSO,ASO, dvso, WSO,WRSO,WASO, &
                     Vlin,Vdep,lambdaHF,Cviso,Vso0,lambdaso,Ccoul, &
                     AAv,BBv,W0,BBs,CCs,Cwiso,Wso0,BBso, &
                     Ea,alpha,eferm,Ades, &
                     rHFl,rHFdep,aHFl,aHFdep,rv,avl,avdep, &
                     rsl,rsdep,as, &
                     rso,aso,rc,ac)
	   RVOL = ACroot * RR
	   RVOL2 = ACroot * RW
	   RSURF = ACroot * RD
	   WRITE(10,10) E,VR,RR,AR, dv,drv,dav, W,RW,AW, VD,RVD,AVD, WD,RD,AD,VSO,RSO,ASO, dvso,WSO,WRSO,WASO, RC,AC
	   10	FORMAT(f7.3, 6(f8.3,2f6.3),2f8.3,2f6.3,2f6.3)
     IF(IEN<=0) THEN
  	    fname='fresco-00-'//POTL//'-s'//CHAR(ICHAR('0')+nexe)//',o'//CHAR(ICHAR('0')+nexo)//'-E0000000.in'
        WRITE(fname(8:8),'(i1)') MOD(NINT(z),10)
  	    WRITE(fname(9:9),'(i1)') MOD(NINT(a),10)
        WRITE(fname(27:33),'(f7.3)') e
        WRITE(fname(27:29),'(i3.3)') INT(e)
        OPEN(1,FORM='formatted',FILE=TRIM(fname))
		    WRITE(0,*) 'Create file <'//TRIM(fname)//'>'
        IF (sump2 .NE. 0) THEN
    	     potname='fresco-00-'//POTL//'-s'//CHAR(ICHAR('0')+nexe)//',o'//CHAR(ICHAR('0')+nexo)//'-E0000000.form'
           WRITE(potname(8:8),'(i1)') MOD(NINT(z),10)
           WRITE(potname(9:9),'(i1)') MOD(NINT(a),10)
     	     WRITE(potname(27:33),'(f7.3)') e
           WRITE(potname(27:29),'(i3.3)') INT(e)
     	     WRITE(0,*) 'Create file <'//TRIM(potname)//'>'
     	     OPEN(94,FORM='formatted',FILE=TRIM(potname))
     	     WRITE(69,156) TRIM(fname),TRIM(potname)
        ELSE
	   	     WRITE(69,156) TRIM(fname)
 		    ENDIF
		    out  = 'n+'//TRIM(NAME)//' with '//TRIM(POTL)//', s='//CHAR(ICHAR('0')+nexe)//',o='
        out = out//CHAR(ICHAR('0')+nexo)//' at E ='//fname(27:33)//' rela='//rela//' Ex*'
		    156  FORMAT(' ',a38,' ',a38)
		    WRITE(1,'(a,f5.2)') out,escale
		    WRITE(1,'(a)') 'NAMELIST'
        IF(E.LE.1.20) THEN
    	     WRITE(1,75) hcm,rmatch
        ELSE
    	     WRITE(1,76) hcm,rmatch
  	    ENDIF
		    WRITE(1,'(a,i4,a,f8.6)') '    jtmin=   0.0 jtmax=',jtmax,' absend= ',absend
		    WRITE(1,14) nex
		    14	FORMAT('    thmin=0.0 thinc=2 thmax=000. iblock=',i3)
		    75  FORMAT(' &Fresco  hcm= ',f6.3,' rmatch= ',f6.3,' rela=''''')
		    76  FORMAT(' &Fresco  hcm= ',f6.3,' rmatch= ',f6.3,' rela=''3d''')
		    WRITE(1,'(a)') '    chans= 1 smats= 2 xstabl= 1'
		    WRITE(1,15) E
		    15	FORMAT('    elab=',f10.3,'  hort=1 /')
		    WRITE(1,*)
 		    WRITE(1,16) nex
		    16	FORMAT('&Partition namep=''n       '' massp=  1.008665 zp=  0 nex=',i3)
		    WRITE(1,17) NAME,A,Z
		    17	FORMAT('            namet=''',a8,''' masst=',f10.6,' zt=',f5.1,' qval=  0.000/')
 		    WRITE(1,'(a)') '&States jp= 0.5 ptyp= 1 ep=  0.000000  cpot=  1 jt= 0.0 ptyt= 1 et= 0.000000 KKt=0/'
 		    DO il = 2,nexe
   		     WRITE(1,21) kpp,jlevels(il),+1,elevels(il)*escale
		    ENDDO
 		    DO il = 1,nexo
   	       WRITE(1,21) kpp,jolevels(il),-2,eolevels(il)*escale
		    ENDDO
        DO il = 1,nexu
   		     WRITE(1,21) kpp,julevels(il),+3,eulevels(il)*escale
		    ENDDO
  	    DO il = 1,nexi
   		     WRITE(1,21) kpp,jilevels(il),+4,eilevels(il)*escale
		    ENDDO
  	    DO il = 1,nna
   		     WRITE(1,222) kpp,jnalevels(il),+5,enalevels(il)*escale
		    ENDDO
	      21	FORMAT('&States copyp= 1                       cpot=',i3,' jt=',f4.1,' ptyt=',i2,' et=',f8.4,' KKt=0/')
		    222	FORMAT('&States copyp= 1                       cpot=',i3,' jt=',f4.1,' ptyt=',i2,' et=',f8.4,' KKt=2/')
 		    WRITE(1,'(a)') '&Partition /'
		    WRITE(1,*)
  	    IF(sump2 .EQ. 0) GO TO 780
  	    CALL FORMFACT(VR,RR,AR,dv,drv,dav,W,RW,AW,VD,RVD,AVD,WD,RD,AD,Ngrid,rmatch,BETA2,BETA4,BETA6, &
                        real(NA),nexo,sump)
		    780  CONTINUE
		    WRITE(1,29) kpp,0,0,REAL(NA),0.0,RC,AC
		    29	FORMAT('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:4)=',4f9.4,'/')
		    30	FORMAT('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:3)=',3f9.4,'/')
		    31	FORMAT('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:6)=',6f9.4,'/')
		    32	FORMAT('&POT /'/)
		    34  FORMAT('&STEP /')
		    WRITE(1,30) kpp,1,0,VR,RR,AR
		    WRITE(1,31) kpp,11,13, 0.0, RVOL*BETA2,0.0, RVOL*BETA4, 0.0, RVOL*BETA6
        IF(sump2 .EQ. 0) GO TO 777
        IF(nexo .EQ. 0) THEN
	         WRITE(1,31) kpp,13,7, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
  	    ELSE IF(nexu .EQ. 0 .AND. nexi .EQ. 0 .AND. nna .EQ. 0 ) THEN
    	     WRITE(1,31) kpp,13,7, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
  	    ELSE
    	     WRITE(1,31) kpp,13,7, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0
  	    ENDIF
  	    IF(nexo .NE. 0) THEN
    	     CALL steps(2,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                        nexe,nexo,nexu,nexi,nna)
        ENDIF
        IF(nexu .NE. 0) THEN
    	     CALL steps(3,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                        nexe,nexo,nexu,nexi,nna)
        ENDIF
  	    IF(nexi .NE. 0) THEN
    	      CALL steps(4,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                        nexe,nexo,nexu,nexi,nna)
        ENDIF
        IF(nna .NE. 0) THEN
    	      CALL steps(5,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       	nexe,nexo,nexu,nexi,nna)
	      ENDIF
  	    WRITE(1,34)
		    777  continue
		    IF(ABS(W)+ABS(dv)>1e-10) THEN
			       WRITE(1,31) kpp,1,0,dv,drv,dav,W,RW,AW
			       WRITE(1,31) kpp,11,13, 0.,BETA2*RVOL2,0.0,BETA4*RVOL2, 0.,BETA6*RVOL2
  		       IF(sump2 .EQ. 0) GO TO 778
  		       IF(nexo .EQ. 0) THEN
	              WRITE(1,31) kpp,13,9, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
  		       ELSE IF(nexu .EQ. 0 .AND. nexi .EQ. 0 .AND. nna .EQ. 0 ) THEN
     		        WRITE(1,31) kpp,13,9, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
  		       ELSE
                WRITE(1,31) kpp,13,9, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0
             ENDIF
  		       IF(nexo .NE. 0) THEN
    	          CALL steps(2,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
               		             nexe,nexo,nexu,nexi,nna)
    	       ENDIF
  		       IF(nexu .NE. 0) THEN
    		        CALL steps(3,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
               	              nexe,nexo,nexu,nexi,nna)
  	         ENDIF
  		       IF(nexi .NE. 0) THEN
    	          CALL steps(4,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                	            nexe,nexo,nexu,nexi,nna)
  	         ENDIF
  		       IF(nna .NE. 0) THEN
  			        CALL steps(5,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                	            nexe,nexo,nexu,nexi,nna)
             ENDIF
    	       WRITE(1,34)
	       ENDIF
		     778  CONTINUE
		     WRITE(1,31) kpp,2,0,VD,RVD,AVD,WD,RD,AD
		     WRITE(1,31) kpp,11,13, 0.0,BETA2*RSURF, 0.0,BETA4*RSURF, 0.0,BETA6*RSURF
         IF(sump2 .EQ. 0) GO TO 779
  	     IF(nexo .EQ. 0) THEN
	 		     WRITE(1,31) kpp,13,9, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
  	     ELSE IF(nexu .EQ. 0 .AND. nexi .EQ. 0 .AND. nna .EQ. 0 ) THEN
           WRITE(1,31) kpp,13,9, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
  	     ELSE
    	     WRITE(1,31) kpp,13,9, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0
  	     ENDIF
  	     IF(nexo .NE. 0) THEN
           CALL steps(2,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
             		       nexe,nexo,nexu,nexi,nna)
  	     ENDIF
  	     IF(nexu .NE. 0) THEN
    	     CALL steps(3,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
             		       nexe,nexo,nexu,nexi,nna)
         ENDIF
  	     IF(nexi .NE. 0) THEN
    	     CALL steps(4,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
               	       nexe,nexo,nexu,nexi,nna)
         ENDIF
  	     IF(nna .NE. 0) THEN
    	      CALL steps(5,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
               	        nexe,nexo,nexu,nexi,nna)
  	     ENDIF
  	     WRITE(1,34)
		     779  CONTINUE
		     WRITE(1,31) kpp,3,0,VSO,RSO,ASO, WSO,WRSO,WASO
		     IF(abs(dvso)>1e-9) WRITE(1,30) kpp,3,0,dvso,WRSO,WASO
  	     WRITE(1,32)
		     WRITE(1,*) '&Overlap /'
		     WRITE(1,*) '&Coupling /'
		     CLOSE(1)
     ENDIF  
   ENDDO
   CLOSE(69)
   CLOSE(40)
   CLOSE(10)
   END
