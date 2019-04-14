PROGRAM test
INTEGER pythonFlag/0/ !False by default
!&Target and &fresco
INTEGER Nenergy,nstat,NBAND
REAL Z,A,eferm
REAL BETA20,BETA40,BETA60
REAL Jmax,hcm,rmatch,Ngrid
REAL Et,Jval,Kval,COEFF
!!!!
! Potential
REAL v0a,v0b,lambdhf,cviso,vspo,lambdso,ccoul
REAL av,bv,w0,bs,wspo,bso,ea,alphav,cs,cwiso,adv
REAL rhfa,rhfb,ahfa,ahfb,rv,ava,avb,rsa,rsb,as,rso,aso,rc,ac 
!!!!! 

! State
REAL,ALLOCATABLE:: elab(:)
REAL,ALLOCATABLE:: Ener_levels(:),J_val(:),KBAND(:),BETA_EFF(:)
INTEGER,ALLOCATABLE:: BAND(:)
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

OPEN(40,FILE='fort.out',STATUS='old') 

!------------------------------------------------------

READ (40, NML=target, END=761, IOSTAT=ios, ERR=761 )
761  IF (ios .ne. 0) then
        WRITE(*,*) 'Input read error while reading Target: ', ios
        STOP
    ENDIF
WRITE(*,target)
ALLOCATE(Ener_levels(nstat),J_val(nstat),KBAND(nstat), &
         BAND(nstat), BETA_EFF(nstat))
DO i=1,nstat
    READ (40, NML=state, END=762, IOSTAT=iosss, ERR=762 )
    762  IF (iosss .ne. 0) THEN
            WRITE(*,*) 'Input read error while reading States: ', iosss
            STOP
        ENDIF
    WRITE(*,state)
    Ener_levels(i)=Et; J_val(i)=Jval; KBAND(i)=Kval; 
    BAND(i)=NBAND; BETA_EFF(i)=COEFF
ENDDO
READ (40, NML=fresco, END=763, IOSTAT=ioss, ERR=763 )
763  IF (ioss .ne. 0) then
        WRITE(*,*) 'Input read error while reading Fresco: ', ioss
        STOP 
    ENDIF
WRITE(*,fresco)

ALLOCATE(elab(Nenergy))
READ (40, NML=energies, END=764, IOSTAT=iosss, ERR=764 )
764  IF (iosss .ne. 0) then
        WRITE(*,*) 'Input read error while reading Energies: ', iosss
        STOP
    ENDIF
WRITE(*,energies)

READ (40, NML=potential, END=765, IOSTAT=iosss, ERR=765 )
765  IF (iosss .ne. 0) then
        WRITE(*,*) 'Input read error while reading potential: ', iosss
        STOP
    ENDIF
WRITE(*,potential)
!-----------------------------------------------------
CLOSE(40)
WRITE(6,*) "Valores de Et: ", Ener_levels
WRITE(6,*) "Valores de J: ", J_val
WRITE(6,*) "Valores de K: ", KBAND
WRITE(6,*) "Valores de N: ", BAND
WRITE(6,*) "Valores de C: ", BETA_EFF 
DEALLOCATE(elab)
DEALLOCATE(Ener_levels,J_val,KBAND,BAND,BETA_EFF)
END PROGRAM test