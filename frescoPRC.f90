!   frescoPRC => Input generator to perform calculations based in the model presented in
!		PRC 94 6 (2016), 064605 for even and odd actinides [effective couplings between
!		bands with dispersive corrections to the optical potential] with FRESCO.
  USE modulo
  PARAMETER (mxsym=100)
  CHARACTER*5 NAME
  CHARACTER*8 POTL
  CHARACTER*100 fname,pname,h,potname
  CHARACTER*20 input_file
  CHARACTER*2 rela,SYMBOL(mxsym)
  !&Target and &fresco
  INTEGER Nenergy,nstat,NBAND,Ngrid,Jmax
  REAL Z,A,eferm
  REAL BETA20,BETA40,BETA60
  REAL hcm,rmatch
  REAL Et,Jval,Kval,COEFF
  ! Potential
  REAL v0a,v0b,lambdhf,cviso,vspo,lambdso,ccoul
  REAL av,bv,w0,bs,wspo,bso,ea,alphav,cs,cwiso,adv
  REAL rhfa,rhfb,ahfa,ahfb,rv,ava,avb,rsa,rsb,as,rso,aso,rc,ac 
  !!!!! 
  INTEGER za,shape,mu,i,sum_neg,sum_pos,gv,ii
  CHARACTER*312 ELEMENTS,out
  CHARACTER*12 pottype(8)
  REAL,ALLOCATABLE:: Ener_levels(:), J_val(:),BETA_EFF(:),BETA_PAR(:),KBAND(:)
  REAL,POINTER:: po
  INTEGER,ALLOCATABLE:: BAND(:),or_val(:), counts(:), indexx(:)
  INTEGER nexe,num_bands
  INTEGER err
  INTEGER pythonFlag/0/ !False by default 
  REAL,TARGET,ALLOCATABLE:: jgsval(:), jexc(:)
  INTEGER,ALLOCATABLE:: index_gs(:), index_exc(:)
  INTEGER,ALLOCATABLE:: gs_Bparity(:), exc_Bparity(:)
  REAL,ALLOCATABLE:: gs_KBAND(:), exc_KBAND(:), elab(:)
  !!!!!
  NAMELIST /target/ Z,A,eferm,BETA20,BETA40,BETA60,nstat
  NAMELIST /fresco/ Jmax,hcm,rmatch,Ngrid,pythonFlag, &
                    Nenergy  
  NAMELIST /energies/ elab
  NAMELIST /state/ Et,Jval,Kval,NBAND,COEFF
  NAMELIST /potential/ v0a,v0b,lambdhf,cviso,vspo,lambdso,ccoul, &
                      av,bv,w0,bs,wspo,bso,ea,alphav,cs,cwiso,adv, &
                      rhfa,rhfb,ahfa,ahfb,rv,ava,avb,rsa,rsb,as, &
                      rso,aso,rc,ac
  !------------------------------------------------------
  data pottype / 'REAL_VOLUME', 'REAL_VOLUME', 'IMAG_VOLUME','REAL_SURFACE', &
     &  		'IMAG_SURFACE', 'REAL_SPINORB', 'IMAG_SPINORB', 'REAL_SPINORB' /
        ELEMENTS =						&
     & ' H HE LI BE  B  C  N  O  F NE NA MG AL SI  P  S CL AR  K CA ' // &
     & 'SC TI  V CR MN FE CO NI CU ZN GA GE AS SE BR KR RB SR  Y ZR ' // &
     & 'NB MO TC RU RH PD AG CD IN SN SB TE  I XE CS BA LA CE PR ND ' // &
     & 'PM SM EU GD TB DY HO ER TM YB LU HF TA  W RE OS IR PT AU HG ' // &
     & 'TL PB BI PO AT RN FR RA AC TH PA  U NP PU AM CM BK CF ES FM'
  READ(ELEMENTS,1021)(SYMBOL(i), i=1, mxsym)
  1021 FORMAT (300(A2,1X))
  450 CONTINUE
  WRITE(6,*) 'Write the name of the input file:'
  READ(*,*) input_file
  OPEN(40,STATUS='old',FILE=input_file,IOSTAT=err)
  IF(err .NE. 0) THEN
    WRITE(6,*) 'Error reading the name of the input file.'
    WRITE(6,*) 'Please, check the name of the file [it must include the extension].'
    GO TO 450
  ENDIF

  !Reading input

  !READ(40,'(F5.1,F10.5,F7.4)') Z,A,eferm
  !READ(40,'(I2)') nstat
  READ (40, NML=target, END=761, IOSTAT=ios, ERR=761 )
761   IF (ios .ne. 0) then
        WRITE(*,*) 'Input read error while reading Target: ', ios
        STOP
      ENDIF
  ALLOCATE(Ener_levels(nstat),J_val(nstat),BAND(nstat),BETA_EFF(nstat), &
  KBAND(nstat),or_val(nstat),indexx(nstat), STAT=err)
  CALL error(err,1)
  DO i=1,nstat
    !READ(40,'(E12.5,F5.2,F4.2,I3,F7.5)') Ener_levels(i),J_val(i),KBAND(i),BAND(i),BETA_EFF(i)
    READ (40, NML=state, END=762, IOSTAT=iosss, ERR=762 )
762     IF (iosss .ne. 0) THEN
            WRITE(*,*) 'Input read error while reading States: ', iosss
            STOP
        ENDIF
    !WRITE(*,state)
    Ener_levels(i)=Et; J_val(i)=Jval; KBAND(i)=Kval; 
    BAND(i)=NBAND; BETA_EFF(i)=COEFF
    indexx(i) = i
  ENDDO
  Ener_levels = Ener_levels/1000 !Reading in KeV but FRESCO reads it in MeV.
  !READ(40,'(I3,F6.2,F7.2,I4,I3)') Jmax,hcm,rmatch,Ngrid,Nenergy
  READ (40, NML=fresco, END=763, IOSTAT=ioss, ERR=763 )
763   IF (ioss .ne. 0) then
        WRITE(*,*) 'Input read error while reading Fresco: ', ioss
        STOP 
      ENDIF
  
  ALLOCATE(elab(Nenergy))
  !READ(40,'(5E12.5)') (elab(i), i=1, Nenergy)
  READ (40, NML=energies, END=764, IOSTAT=iosss, ERR=764 )
764  IF (iosss .ne. 0) then
        WRITE(*,*) 'Input read error while reading Energies: ', iosss
        STOP
     ENDIF
  ccoul = 1.36 ! Default
  rc = 1.2894 !Default
  ac = 0.547 !Default 
                        !Dispersive parameters
  READ (40, NML=potential, END=765, IOSTAT=iosss, ERR=765 )
765   IF (iosss .ne. 0) then
        WRITE(*,*) 'Input read error while reading potential: ', iosss
        STOP
      ENDIF
  !WRITE(*,potential)
  !READ(40,'(7E12.4)') v0a,v0b,lambdhf,cviso,vspo,lambdso,ccoul
  !READ(40,'(6E12.4)') av,bv,w0,bs,Wso0,BBso
  !READ(40,'(5E12.4)') Ea,alpha,CCs,Cwiso,Ades
  !READ(40,'(7E12.4)') rHFl,rHFdep,aHFl,aHFdep,rv,avl,avdep
  !READ(40,'(3E12.4)') rsl,rsdep,as
  !!READ(40,'(4E12.4)') rso,aso,rc,ac
  !READ(40,'(3E12.4)') BETA20,BETA40,BETA60
  !READ(40,'(I1)') pythonFlag
  BETA_EFF = BETA_EFF/BETA20 ! input like OPTMAN, with x\beta_{20} factor.

  !Reading done

  NAME = symbol(NINT(Z))//'000'
  WRITE(6,*) 'Z,A,name =',NINT(Z),NINT(A),symbol(NINT(Z))
  WRITE(NAME(3:5),'(i3.3)') NINT(A)
  absend = 0.001
  POTL = 'DOMEIC16'
  kpp = 1
  IF(NAME(1:1) ==' ') NAME(1:5)=NAME(2:5)//' '
  pname = POTL//'-'//TRIM(NAME)//'-parameters.txt'
  OPEN(10,FORM='formatted',FILE=TRIM(pname))
  OPEN(69,FORM='formatted',FILE='lista.txt' ) !Auxiliar .txt to run FRESCO's inputs for different energies.
  WRITE(10,1) POTL,NAME
  1	FORMAT('####### OPTICAL PARAMETERS for ',A8,' ########'/'#'/'#     neutron on ',a5,/'#'/ &
        & '#Energy    V     rv    av     dV    drv   dav      W     rw    aw     ', &
        &           ' Vd   rvd   avd      Wd   rwd   awd     ', &
        &           'Vso   rvso  avso    dvso    Wso   rwso  awso  rc    ac')
  NA = NINT(A); NZ = NINT(Z)
  ACroot = real(NA)**(1./3.)
 
                      

  ! Separation of parameters for G.S band  and excited bands
  !///////////////////////////////////////////////////////////////////////////////////
  num_bands = 1
  or_val(1) = BAND(1) !Identify wich states are in a same band
  outer: DO i=2,nstat
    DO j=1,num_bands
      IF(or_val(j)==BAND(i)) THEN
        CYCLE outer
      ENDIF
    ENDDO
    num_bands = num_bands+1
    or_val(num_bands) = BAND(i) ! different NBAND values.
  ENDDO outer
  ALLOCATE(counts(num_bands),STAT=err)
  CALL error(err,1)
  counts = 0
  DO i=1,num_bands
    DO j=1,nstat
      IF(or_val(i)==BAND(j)) THEN
        counts(i) = counts(i)+1 !Number of states per band -> counts(1)= G.S
      ENDIF
    ENDDO
  ENDDO
  nexe = counts(1); n_exc = nstat-nexe
  ALLOCATE(jgsval(nexe),index_gs(nexe),gs_Bparity(nexe), &
  gs_KBAND(nexe),STAT=err)
  CALL error(1,err)
  ALLOCATE(jexc(n_exc),index_exc(n_exc),exc_Bparity(n_exc), &
  exc_KBAND(n_exc),BETA_PAR(n_exc),STAT=err)
  CALL error(1,err)
  gv = 1; ii = 1 !indices
  DO i=1,nstat
    IF(ABS(BAND(i))==1) THEN ! G.S band MUST be |NBAND|=1 in the input.
      po=>jgsval(gv)
      po = J_val(i)
      index_gs(gv) = indexx(i)
      gs_Bparity(gv) = BAND(i)
      gs_KBAND(gv) = KBAND(i)
      gv = gv+1
    ELSE !Saving all the infomation about excited band's states.
      po=>jexc(ii)
      po = J_val(i)
      index_exc(ii) = indexx(i)
      exc_Bparity(ii) = BAND(i)
      exc_KBAND(ii) = KBAND(i)
      BETA_PAR(ii) = BETA_EFF(i)
      ii = ii+1
    ENDIF
  ENDDO
  sum_pos = 0; sum_neg = 0 !Number of excited states with positive/negative parity.
  DO j=1,n_exc
    IF(exc_Bparity(j) .LT. 0) THEN
      sum_neg = sum_neg+1
    ELSE IF(exc_Bparity(j) .GT. 0) THEN
      sum_pos = sum_pos+1
    ENDIF
  ENDDO
  !///////////////////////////////////////////////////////////////////////////////////
  DO k=1,Nenergy
     E = elab(k)
     hcmv = hcm
     IF(E>50.0) hcmv = hcm/SQRT(E/50.0)
     dv = 0; drv = 0; dav = 0; dvso = 0.
     NTYPE = 1 ! neutron only!!
     CALL dispers2(A,Z,NTYPE,E,VR,RR,AR, dv,drv,dav, VD,RVD,AVD, &
                     W,RW,AW, WD,RD,AD, VSO,RSO,ASO, dvso, WSO,WRSO,WASO, &
                     v0a,v0b,lambdhf,cviso,vspo,lambdso,ccoul, &
                     av,bv,w0,bs,cs,cwiso,wspo,bso, &
                     ea,alphav,eferm,adv, &
                     rhfa,rhfb,ahfa,ahfb,rv,ava,avb, &
                     rsa,rsb,as, &
                     rso,aso,rc,ac)
     RVOL = ACroot * RR
     RVOL2 = ACroot * RW
     RSURF = ACroot * RD
     WRITE(10,10) E,VR,RR,AR, dv,drv,dav, W,RW,AW, VD,RVD,AVD, WD,RD,AD,VSO,RSO,ASO, dvso,WSO,WRSO,WASO, RC,AC
     10	FORMAT(f7.3, 6(f8.3,2f6.3),2f8.3,2f6.3,2f6.3)
     fname = 'fresco-00-'//POTL//'-s'//CHAR(ICHAR('0')+nexe)//',o'//CHAR(ICHAR('0')+sum_neg)//'-E0000000.in'
     WRITE(fname(8:9),'(i2)') NINT(Z) ! Z=>10
     WRITE(fname(27:33),'(f7.3)') e
     WRITE(fname(27:29),'(i3.3)') INT(e)
     OPEN(1,FORM='formatted',FILE=TRIM(fname))
     WRITE(0,*) 'Create file <'//TRIM(fname)//'>'
     IF(n_exc .NE. 0) THEN !If the number of states in excited bands is 0 then .form files are not necessary.
       potname='fresco-00-'//POTL//'-s'//CHAR(ICHAR('0')+nexe)//',o'//CHAR(ICHAR('0')+sum_neg)//'-E0000000.form'
       WRITE(potname(8:9),'(i2)') NINT(Z) !Z=>10
       WRITE(potname(27:33),'(f7.3)') e
       WRITE(potname(27:29),'(i3.3)') INT(e)
       WRITE(0,*) 'Create file <'//TRIM(potname)//'>'
       OPEN(94,FORM='formatted',FILE=TRIM(potname))
       WRITE(69,156) TRIM(fname),TRIM(potname)
     ELSE
       WRITE(69,156) TRIM(fname)
     ENDIF
     out  = 'n+'//TRIM(NAME)//' with '//TRIM(POTL)//''
     156 FORMAT(' ',a38,' ',a38)
     WRITE(1,'(a38)') out
     WRITE(1,'(a)') 'NAMELIST'
!    frxy6j version of FRESCO has some problems using RELA=3d for low energies. After some tests, 1.20 MeV limit seems to work.
!    Feb/2019 -> Antonio fixed the problem with low energies when rela=3d for frxy6j and included it in FRESCO v3.
!    New version of FRESCO "v3.2" can be used with rela=c for any energy but "hort" still not available for this version.
     IF(E.LE.1.20) THEN
       WRITE(1,75) hcmv,rmatch ! 76 --> 75 if using not fixed frxy6j.
     ELSE
       WRITE(1,76) hcmv,rmatch
     ENDIF
!!!!!!!!!
     WRITE(1,'(a,i4,a,f8.6)') '    jtmin=   0.0 jtmax=',Jmax,' absend= ',absend
     WRITE(1,14) nstat
     14	FORMAT('    thmin=0.0 thinc=2 thmax=000. iblock=',i3)
     75 FORMAT(' &Fresco  hcm= ',f6.3,' rmatch= ',f6.3,' rela=''''')
     76 FORMAT(' &Fresco  hcm= ',f6.3,' rmatch= ',f6.3,' rela=''3d''') !RELA=3d(frxy6j) or RELA=c (v32) is neccesary to agree with OPTMAN.
     !WRITE(1,'(a)') '    chans= 1 smats= 2 xstabl= 1' !extra output information.
     WRITE(1,15) E
     15	FORMAT('    elab=',f10.3,' /') ! hort might be neccesary for large number of coupled states -> only in frxy6j for now.
     WRITE(1,*)
     WRITE(1,16) nstat
     16	FORMAT('&Partition namep=''n       '' massp=  1.008665 zp=  0 nex=',i3)
     WRITE(1,17) NAME,A,Z
     17	FORMAT('            namet=''',a8,''' masst=',f10.5,' zt=',f5.1,' qval=  0.000/')
     451 FORMAT('&States jp= 0.5 ptyp= 1 ep=  0.000000  cpot=  1 jt=',f4.1,' ptyt=',i2,' et=',f8.4,' KKt=',f3.1,'/')
     WRITE(1,451) J_val(1),BAND(1),Ener_levels(1),KBAND(1)
     DO i=2,nstat ! First state (i=1) must be target's ground state.
       WRITE(1,21) 1,J_val(i),BAND(i),Ener_levels(i),KBAND(i)
     ENDDO
     21	FORMAT('&States copyp= 1                       cpot=',i3,' jt=',f4.1,' ptyt=',i2,' et=',f8.4,' KKt=',f3.1,'/')
     WRITE(1,'(a)') '&Partition /'
     WRITE(1,*)
     IF(n_exc .EQ. 0) GO TO 780 !if the number of states in excited bands is 0 then .form files are not generated.
     CALL FORMFACT(VR,RR,AR,dv,drv,dav,W,RW,AW,VD,RVD,AVD,WD,RD,AD,Ngrid,rmatch,BETA20,BETA40,BETA60, &
                  REAL(NA),sum_neg,sum_pos)
     780 CONTINUE
     WRITE(1,29) kpp,0,0,REAL(NA),0.0,RC,AC
     29	FORMAT('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:4)=',4f9.4,'/')
     30 FORMAT('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:3)=',3f9.4,'/')
     31 FORMAT('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:6)=',6f9.4,'/')
     32 FORMAT('&POT /'/)
     34 FORMAT('&STEP /')
     WRITE(1,30) kpp,1,0,VR,RR,AR
     WRITE(1,31) kpp,11,13, 0.0, RVOL*BETA20,0.0, RVOL*BETA40, 0.0, RVOL*BETA60
     IF(n_exc .EQ. 0) GO TO 777
     IF(sum_neg .EQ. 0) THEN   !&pot type=13 is different in function of the number of positive and negative parity states.
       WRITE(1,31) kpp,13,7, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
     ELSE IF(sum_pos .EQ. 0) THEN
       WRITE(1,31) kpp,13,7, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
     ELSE
       WRITE(1,31) kpp,13,7, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0
     ENDIF
     IF(MOD(INT(A),2) .EQ. 0) THEN  ! Checking if we are in the case of even or odd target -> different expressions.
       DO i=1,nexe
         DO j=1,n_exc
           CALL steps_even(jgsval(i),jexc(j),gs_Bparity(i), exc_Bparity(j), &
                           index_gs(i),index_exc(j),gs_KBAND(i),exc_KBAND(j),BETA_PAR(j))
         ENDDO
       ENDDO
     ELSE
       DO i=1,nexe
         DO j=1,n_exc
           CALL steps_odd(jgsval(i),jexc(j),gs_Bparity(i), exc_Bparity(j), &
                          index_gs(i),index_exc(j),gs_KBAND(i),exc_KBAND(j),BETA_PAR(j))
         ENDDO
       ENDDO
     ENDIF
     WRITE(1,34)
     777 CONTINUE
     IF(ABS(W)+ABS(dv)>1e-10) THEN
       WRITE(1,31) kpp,1,0,dv,drv,dav,W,RW,AW
       WRITE(1,31) kpp,11,13, 0.,BETA20*RVOL2,0.0,BETA40*RVOL2, 0.,BETA60*RVOL2
       IF(n_exc .EQ. 0) GO TO 778
       IF(sum_neg .EQ. 0) THEN !&pot type=13 is different in function of the number of positive and negative parity states.
         WRITE(1,31) kpp,13,9, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
       ELSE IF(sum_pos .EQ. 0) THEN
         WRITE(1,31) kpp,13,9, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
       ELSE
         WRITE(1,31) kpp,13,9, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0
       ENDIF
       IF(MOD(INT(A),2) .EQ. 0) THEN
         DO i=1,nexe
           DO j=1,n_exc
             CALL steps_even(jgsval(i),jexc(j),gs_Bparity(i), exc_Bparity(j), &
                             index_gs(i),index_exc(j),gs_KBAND(i),exc_KBAND(j),BETA_PAR(j))
           ENDDO
         ENDDO
       ELSE
         DO i=1,nexe
           DO j=1,n_exc
             CALL steps_odd(jgsval(i),jexc(j),gs_Bparity(i), exc_Bparity(j), &
                            index_gs(i),index_exc(j),gs_KBAND(i),exc_KBAND(j),BETA_PAR(j))
           ENDDO
         ENDDO
       ENDIF
       WRITE(1,34)
     ENDIF
     778 CONTINUE
     WRITE(1,31) kpp,2,0,VD,RVD,AVD,WD,RD,AD
     WRITE(1,31) kpp,11,13, 0.0,BETA20*RSURF, 0.0,BETA40*RSURF, 0.0,BETA60*RSURF
     IF(n_exc .EQ. 0) GO TO 779
     IF(sum_neg .EQ. 0) THEN !&pot type = 13 is different in function of the number of positive and negative parity states.
       WRITE(1,31) kpp,13,9, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
     ELSE IF(sum_pos .EQ. 0 ) THEN
       WRITE(1,31) kpp,13,9, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
     ELSE
       WRITE(1,31) kpp,13,9, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0
     ENDIF
     IF(MOD(INT(A),2) .EQ. 0) THEN
       DO i=1,nexe
         DO j=1,n_exc
           CALL steps_even(jgsval(i),jexc(j),gs_Bparity(i), exc_Bparity(j), &
                           index_gs(i),index_exc(j),gs_KBAND(i),exc_KBAND(j),BETA_PAR(j))
         ENDDO
       ENDDO
     ELSE
       DO i=1,nexe
         DO j=1,n_exc
           CALL steps_odd(jgsval(i),jexc(j),gs_Bparity(i), exc_Bparity(j), &
                          index_gs(i),index_exc(j),gs_KBAND(i),exc_KBAND(j),BETA_PAR(j))
         ENDDO
       ENDDO
     ENDIF
     WRITE(1,34)
     779 CONTINUE
     WRITE(1,31) kpp,3,0,VSO,RSO,ASO, WSO,WRSO,WASO
     IF(ABS(dvso)>1e-9) WRITE(1,30) kpp,3,0,dvso,WRSO,WASO
     WRITE(1,32)
     WRITE(1,*) '&Overlap /'
     WRITE(1,*) '&Coupling /'
     CLOSE(1)
   ENDDO
   IF(pythonFlag .EQ. 1) THEN !generating graphs.py in case of graphs generation at the end of the run using runall.sh
     CALL Graphs
     WRITE(6,*) 'graphs.py created: C.S graphs will be generated using Python (matplotlib)'
   ENDIF
   DEALLOCATE(Ener_levels,J_val,BAND,BETA_EFF,KBAND,or_val, &
   counts,indexx,STAT=err)
   CALL error(err,2)
   DEALLOCATE(jgsval,index_gs,jexc,index_exc,BETA_PAR, &
   gs_Bparity,exc_Bparity,gs_KBAND,exc_KBAND,elab,STAT=err)
   CALL error(err,2)
   CLOSE(69)
   CLOSE(40)
   CLOSE(10)
   END
