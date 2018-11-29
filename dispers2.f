      SUBROUTINE dispers2(A,Z,k,eopt,
     #  v,rvv,avv, dv,drv,dav, dvs,drs,das, w,rw,aw, wd,rwd,awd,
     #  vso,rvso,avso,dvso, wso,rwso,awso, 
     #  Vlin,Vdep,lambdaHF,Cviso,Vso0,lambdaso,Ccoul,
     #  AAv,BBv,W0,BBs,CCs,Cwiso,Wso0,BBso,
     #  Ea,alpha,eferm,Ades,
     #  rHFl,rHFdep,aHFl,aHFdep,rv,avl,avdep,
     #  rsl,rsdep,as,
     #  rso,aso,rc,ac)
      	real    eopt,asym,eferm,f,Cviso,viso,Ccoul,Cwiso
	real lambdaHF,lambdaso,Ades
    	pi = 4.0*atan(1.)
c
c *** Parameters of Soukhovitskii, Capote, Quesada, Chiba and Martyanov (Nov 25, 2015) ***
c *** with Asymmetrical W energy-dependence
c *** dispers2:  T1 integral with correct coefficient
c *** by Ian Thompson 
c k         : designator for particle
c Z         : charge number of residual nucleus
c A         : mass number of residual nucleus
c eopt      : incident energy
c asym      : asymmetry parameter
c eferm     : Fermi energy
c f         : eopt-eferm
	
	Au = A-Ades
	asym=(A-2.*Z)/A

	V0 = Vlin + Vdep*Au 
	
	rHF = rHFL + rHFdep * Au 
	aHF = aHFl + aHFdep * Au  			
	av = avl + avdep * Au
	rs = rsl + rsdep * Au
	AAHF = V0 * (1 + (-1)**k * Cviso*asym/V0)
	AAs  = W0 * (1 + (-1)**k * Cwiso*asym/W0)

! Energy-dependent quantities;
	
	eoffset = 0.
	if (k==2) eoffset =  Ccoul * Z/A**(1./3.)
    	Eeff = eopt - eoffset
	f = Eeff  - eferm

! Non-dispersive 
	v = AAHF * exp(-lambdaHF*f)
	vso=Vso0*exp(-lambdaso*f)

! sources of dispersive terms
	w = AAv * f*f/(f*f + BBv**2)
    	if(f < -Ea) then         
        	fe  = f + Ea
        	w = w * (1 - fe*fe/(fe*fe + Ea*Ea))
    	else if(f>Ea) then       
        	w = w + alpha * (sqrt(Eeff) + (eferm+Ea)**1.5d0/(2*Eeff)
     x          	- 1.5d0*sqrt(eferm+Ea))
   	endif
	wd = AAs * f*f/(f*f + BBs**2) * exp( -CCs * abs(f))
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

        !write(7,100) eopt,DWv,T1,T2
        !write(8,100) eopt,w,wd,wso
        !write(9,100) eopt,dv,dvs,dvso
100     format(f12.3, 8f10.5)
	return
	end
