!          neutron  on actinides
	parameter (mxsym=100)
        character*5 NAME
	character*8 POTL
        DOUBLE PRECISION CG 
	character*100 fname,pname,h,potname
	character*2 rela,SYMBOL(mxsym),PT(4)
	real elevels(20), jlevels(20),enodes(10),val(4,10,8), enalevels(20),jnalevels(20)
	real eolevels(20), jolevels(20)
        real eulevels(20), julevels(20)
        real eilevels(20), jilevels(20),lambdaHF,lambdaso
	integer za,shape,kpnodes(10),fnodes,mu,i
	character*312 ELEMENTS,out
	character tab
	character*12 pottype(8)
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
1021	FORMAT (300(A2,1X))
        open(40,status='old',file='input.INP')
	read(40,*) Z,A,eferm
        read(40,*) nexe,nexo,nexu,nexi,nna
        NAME = symbol(nint(Z))//'000'       
	write(6,*) 'Z,A,name =',nint(Z),nint(A),symbol(nint(Z))
        write(NAME(3:5),'(i3.3)') nint(A)
        read(40,*) jtmax,hcm0,rmatch,Ngrid 
	absend=0.001
	elevels(1:10) =(/ 0., .04490, .14840, 0.3074, 0.5174, 0.7759, 1.0767, 1.4155, 1.7884, 2.1911/)  ! G.S band 
	jlevels(1:10) =(/ 0., 2., 4., 6., 8., 10., 12., 14., 16., 18. /)
	eolevels(1:5) =(/ 0.6798, .7313, .8267, 0.9664, 1.1507  /) ! Octopule band - parity (k=0-)
	jolevels(1:5) =(/ 1., 3., 5., 7., 9. /)
        eulevels(1:2)=(/0.993, 1.0373/) !beta band
        julevels(1:2)=(/0., 2./)
        eilevels(1:3)=(/0.9270, .9663, 1.0564/) !gamma band
        jilevels(1:3)=(/0.,2.,4./)
        enalevels(1:2)=(/1.0603,1.1057/)
        jnalevels(1:2)=(/2.0,3.0/)
        POTL= 'dispers2'
	!nexe =6
	!nexo =5
        !nexu=2
        !nexi=3
        !nna=2
        escale=1.0
	nex = nexe + nexo + nexu + nexi + nna
        sump = nexu + nexi + nna !auxiliar para FORM
        read(40,*) EMIN,EMAX,NE          
        nodes = 7
	fnodes = 6 
        enodes(1:7) = (/ 0.01, 1., 5., 10., 20., 50. , 200./)
        kpnodes(1:7) = (/ 10, 11, 15, 20, 30, 32, 34 /)
	kpp=1
	if(NAME(1:1) ==' ') NAME(1:5)=NAME(2:5)//' '
	
	pname = POTL//'-'//trim(NAME)//'-parameters.txt'
        open(10,form='formatted',file=trim(pname))
        open(69,form='formatted',file='lista.txt' )
	write(10,1) POTL,NAME
1	format('####### OPTICAL PARAMETERS for ',A8,' ########'/'#'/'#     neutron on ',a5,/'#'/ &
     & '#Energy    V     rv    av     dV    drv   dav      W     rw    aw     ', & 
     &           ' Vd   rvd   avd      Wd   rwd   awd     ', &
     &           'Vso   rvso  avso    dvso    Wso   rwso  awso  rc    ac')		
	DE = 1.
	IF(NE>1) DE = (EMAX-EMIN)/(NE-1)
	IF(NE<-1) DE = (log(EMAX)-log(EMIN))/(abs(NE)-1)
      	NA = nint(A); NZ = nint(Z)
	ACroot = real(NA)**(1./3.)
        Ccoul=1.36
        rc=1.2894
        ac=0.547
        !Lectura de parámetros dispersivos
        read(40,*) Vlin,Vdep,lambdaHF,Cviso,Vso0,lambdaso,Ccoul
        read(40,*) AAv,BBv,W0,BBs,CCs,Cwiso,Wso0,BBso
        read(40,*) Ea,alpha,Ades
        read(40,*) rHFl,rHFdep,aHFl,aHFdep,rv,avl,avdep
        read(40,*) rsl,rsdep,as
        read(40,*) rso,aso,rc,ac
        !Lectura de deformaciones => MULTIPLICADAS POR BETA2 COMO EL INPUT DE OPTMAN
        read(40,*) BETA2,BETA4,BETA6
        read(40,*) BETA2EFF,GAMMA2EFF,GAMMANAX,BETA3
	BETA3 = BETA3/(BETA2)!Paper _!!!!!!
	BETA2EFF=BETA2EFF/(BETA2)!GS>BETA-0+
        GAMMA2EFF=GAMMA2EFF/(BETA2)!GS>GAMMA BAND(mu=0)
        GAMMANAX=GAMMANAX/(BETA2) !!  
     	DO IE=1,abs(NE)!+nodes
	 IEN = IE-abs(NE)
	if(IE<=abs(NE)) then
	  if(NE>1) E = EMIN + (IE-1)*DE
	  if(NE<-1) E = exp(log(EMIN) + (IE-1)*DE)
	 else
	  E = enodes(IEN)
	  kp = kpnodes(IEN)
	 endif
	hcm = hcm0
	if(E>50.0) hcm = hcm0/sqrt(E/50.0)
	dv=0; drv=0; dav=0; dvso=0.
        NTYPE = 1 ! neutron!!
                         ! Parámetros del potencial !! (incluidas las correcciones dispersivas)
        call dispers2(A,Z,NTYPE,E,VR,RR,AR, dv,drv,dav, VD,RVD,AVD, &
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
	write(10,10) E,VR,RR,AR, dv,drv,dav, W,RW,AW, VD,RVD,AVD, WD,RD,AD,VSO,RSO,ASO, dvso,WSO,WRSO,WASO, RC,AC
10	format(f7.3, 6(f8.3,2f6.3),2f8.3,2f6.3,2f6.3)

      if(IEN<=0) then  

      fname='fresco-00-'//POTL//'-s'//CHAR(ICHAR('0')+nexe)//',o'//CHAR(ICHAR('0')+nexo)//'-E0000000.in'
      potname='fresco-00-'//POTL//'-s'//CHAR(ICHAR('0')+nexe)//',o'//CHAR(ICHAR('0')+nexo)//'-E0000000.form'
        write(fname(8:8),'(i1)') mod(nint(z),10)
        write(fname(9:9),'(i1)') mod(nint(a),10)
        write(fname(27:33),'(f7.3)') e
        write(fname(27:29),'(i3.3)') int(e)
        write(0,*) 'Create file <'//trim(fname)//'>'
        write(potname(8:8),'(i1)') mod(nint(z),10)
        write(potname(9:9),'(i1)') mod(nint(a),10)
        write(potname(27:33),'(f7.3)') e
        write(potname(27:29),'(i3.3)') int(e)
        write(0,*) 'Create file <'//trim(potname)//'>'
        open(1,form='formatted',file=trim(fname))
        open(94,form='formatted',file=trim(potname))
156     format(' ',a38,' ',a38) 
        write(69,156) trim(fname),trim(potname)
	out  = 'n+'//trim(NAME)//' with '//trim(POTL)//', s='//CHAR(ICHAR('0')+nexe)//',o='
        out = out//CHAR(ICHAR('0')+nexo)//' at E ='//fname(27:33)//' rela='//rela//' Ex*'
	write(1,'(a,f5.2)') out,escale
	write(1,'(a)') 'NAMELIST'
        if(E.le.1.20) then 
     	     !write(1,'(a,f6.3,a)') ' &Fresco  hcm= ',hcm,' rmatch=  20.000 rela='''''
             write(1,75) hcm,rmatch
        else
             !write(1,'(a,f6.3,a)') ' &Fresco  hcm= ',hcm,' rmatch=  20.000 rela=''3d'''
             write(1,76) hcm,rmatch
        endif	
	write(1,'(a,i4,a,f8.6)') '    jtmin=   0.0 jtmax=',jtmax,' absend= ',absend
	write(1,14) nex
14	format('    thmin=0.0 thinc=2 thmax=000. iblock=',i3)
75      format (' &Fresco  hcm= ',f6.3,' rmatch= ',f6.3,' rela=''''')
76      format (' &Fresco  hcm= ',f6.3,' rmatch= ',f6.3,' rela=''3d''')
	write(1,'(a)') '    chans= 1 smats= 2 xstabl= 1'
	write(1,15) E
15	format('    elab=',f10.5,'  pel=1 exl=1 lab=1 lin=1 lex=1 hort=1 /')
	write(1,*) 
 	write(1,16) nex
16	format('&Partition namep=''n       '' massp=  1.008665 zp=  0 nex=',i3)
	write(1,17) NAME,A,Z
17	format('            namet=''',a8,''' masst=',f10.6,' zt=',f5.1,' qval=  0.000/')
 	write(1,'(a)') '&States jp= 0.5 ptyp= 1 ep=  0.000000  cpot=  1 jt= 0.0 ptyt= 1 et= 0.000000 KKt=0/'
 	do il = 2,nexe
   	  write(1,21) kpp,jlevels(il),+1,elevels(il)*escale
	enddo
 	do il = 1,nexo
   	  write(1,21) kpp,jolevels(il),-2,eolevels(il)*escale
	enddo
        do il = 1,nexu
   	  write(1,21) kpp,julevels(il),+3,eulevels(il)*escale
	enddo
        do il = 1,nexi
   	  write(1,21) kpp,jilevels(il),+4,eilevels(il)*escale
	enddo
        do il = 1,nna
   	  write(1,222) kpp,jnalevels(il),+5,enalevels(il)*escale
	enddo
21	format('&States copyp= 1                       cpot=',i3,' jt=',f4.1,' ptyt=',i2,' et=',f8.4,' KKt=0/')
222	format('&States copyp= 1                       cpot=',i3,' jt=',f4.1,' ptyt=',i2,' et=',f8.4,' KKt=2/')
 	write(1,'(a)') '&Partition /'
	write(1,*) 
            call FORMFACT(VR,RR,AR,dv,drv,dav,W,RW,AW,VD,RVD,AVD,WD,RD,AD,Ngrid,rmatch,BETA2,BETA4,BETA6, &
                          real(NA),nexo,sump)
	write(1,29) kpp,0,0,real(NA),0.0,RC,AC !COULOMB
29	format('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:4)=',4f9.4,'/')
30	format('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:3)=',3f9.4,'/')
31	format('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:6)=',6f9.4,'/')
32	format('&POT /'/)
34      format('&STEP /')
	write(1,30) kpp,1,0,VR,RR,AR !NUCLEAR VOLUMEN (campo medio)
	write(1,31) kpp,11,13, 0.0, RVOL*BETA2,0.0, RVOL*BETA4, 0.0, RVOL*BETA6
        if(nexo .eq. 0) then 
	     write(1,31) kpp,13,7, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
        else if(nexu .eq. 0 .and. nexi .eq. 0 .and. nna .eq. 0 ) then
             write(1,31) kpp,13,7, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
        else
             write(1,31) kpp,13,7, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0
        endif
        if(nexo .ne. 0) then
            call steps(2,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        if(nexu .ne. 0) then
            call steps(3,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        if(nexi .ne. 0) then
            call steps(4,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        if(nna .ne. 0) then
            call steps(5,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        write(1,34)
	if(abs(W)+abs(dv)>1e-10) then
	write(1,31) kpp,1,0,dv,drv,dav,W,RW,AW !VOLUMEN real (DISPERSIVO) + imaginario
	write(1,31) kpp,11,13, 0.,BETA2*RVOL2,0.0,BETA4*RVOL2, 0.,BETA6*RVOL2 !DEFORMADO ROTOR VOLUMEN
        if(nexo .eq. 0) then 
	     write(1,31) kpp,13,9, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
        else if(nexu .eq. 0 .and. nexi .eq. 0 .and. nna .eq. 0 ) then
             write(1,31) kpp,13,9, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
        else
             write(1,31) kpp,13,9, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0
        endif
        if(nexo .ne. 0) then
            call steps(2,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        if(nexu .ne. 0) then
            call steps(3,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        if(nexi .ne. 0) then
            call steps(4,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        if(nna .ne. 0) then
            call steps(5,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        write(1,34)
	endif
	write(1,31) kpp,2,0,VD,RVD,AVD,WD,RD,AD !SUPERFICIE real (DISPERSIVO) + imaginario
	write(1,31) kpp,11,13, 0.0,BETA2*RSURF, 0.0,BETA4*RSURF, 0.0,BETA6*RSURF !DEFORMADO ROTOR VOLUMEN
        if(nexo .eq. 0) then 
	     write(1,31) kpp,13,9, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
        else if(nexu .eq. 0 .and. nexi .eq. 0 .and. nna .eq. 0 ) then
             write(1,31) kpp,13,9, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
        else
             write(1,31) kpp,13,9, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0
        endif     
        if(nexo .ne. 0) then
            call steps(2,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        if(nexu .ne. 0) then
            call steps(3,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        if(nexi .ne. 0) then
            call steps(4,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        if(nna .ne. 0) then
            call steps(5,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                       nexe,nexo,nexu,nexi,nna)
        endif
        write(1,34)
	write(1,31) kpp,3,0,VSO,RSO,ASO, WSO,WRSO,WASO !S.O real (DISPERSIVO) + imaginario
	if(abs(dvso)>1e-9) write(1,30) kpp,3,0,dvso,WRSO,WASO    
        write(1,32) 
	write(1,*) '&Overlap /'
	write(1,*) '&Coupling /'
	close(1)
	endif  !IEN>0
	enddo
        close(69)
        close(40)
	end
