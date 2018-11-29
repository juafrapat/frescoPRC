  double precision function arm(x,L) !L=lambda, x=cos(tetha) !correcto 24/11
        !funciona¡ => Armónico esférico L,x con m=0
	double precision x,P(0:L),pi
        integer L,i
        pi=acos(-1d0)
        if(L==0) then
	   arm=dsqrt(1.d0/(4d0*pi))
        endif
        if(L==1) then
	   arm=dsqrt(3.d0/(4d0*pi))*x
        endif
        if(L.ge.2) then	 
	    P(0)=1d0
	    P(1)=x
      	    DO i=2,L,1
              P(i)=((2d0*DBLE(i)-1d0)*x*P(i-1) - DBLE(i-1)*P(i-2))/(DBLE(i))
 	    END DO
            arm=dsqrt((2d0*DBLE(L)+1d0)/(4d0*pi))*P(L)
	endif
       return
       end		
 double precision function rd(R,x,BETA2,BETA4,BETA6)
       real R,BETA2,BETA4,BETA6
       double precision x,arm
       rd=DBLE(R)*(1d0+DBLE(BETA2)*arm(x,2)+DBLE(BETA4)*arm(x,4)+DBLE(BETA6)*arm(x,6))
       return
       end
 double precision function vws(r,v0,RR,a,x,L,BETA2,BETA4,BETA6)
        real r,v0,RR,a,BETA2,BETA4,BETA6
        double precision x,rd,arm
        integer L
       vws=(DBLE(v0)/(1d0+exp((DBLE(r)-rd(RR,x,BETA2,BETA4,BETA6))/DBLE(a))))*arm(x,L)
	return
        end
 double precision function vs(r,v0,RR,a,x,L,BETA2,BETA4,BETA6)
        real r,v0,RR,a,BETA2,BETA4,BETA6
        double precision x,rd,arm,f
        integer L
        f=exp((DBLE(r)-rd(RR,x,BETA2,BETA4,BETA6))/DBLE(a))
	vs=-4d0*DBLE(v0)*arm(x,L)*f/(1+f)**2
        return
        end
 double precision function gaussv(r,v0,RR,a,L,BETA2,BETA4,BETA6)
        real r,v0,RR,a,x(19),w(19),BETA2,BETA4,BETA6
        double precision pp,vws,pi,p
        integer j,L
        pi=acos(-1d0)
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
        gaussv=dsqrt(pi)*p !Factor para coincidir con FRESCO
        return
        end
 double precision function gausss(r,v0,RR,a,L,BETA2,BETA4,BETA6)
        real r,v0,RR,a,x(19),w(19),BETA2,BETA4,BETA6
        double precision pp,vs,pi,p
        integer j,L
        pi=acos(-1d0)
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
        gausss=-dsqrt(pi)*p !signo y factor por convenio para coincidir con FRESCO
        return
        end
  SUBROUTINE FORMFACT(VR,RR,AR,dv,drv,dav,W,RW,AW,VD,RVD,AVD,WD,RD,AD,N,rmax,BETA2,BETA4,BETA6,A,nexo,sump)
   double precision vws,gaussv,gausss,deltar,r(0:N)  
   real rmax,BETA2,BETA4,BETA6,A 
   real VR,RR,AR !Volumen campo medio solo real
   real dv,drv,dav,W,RW,AW !Volumen dispersivo real/imaginario
   real VD,RVD,AVD,WD,RD,AD
   integer N,i,nexo,sump !Si nexo=0 o sump=0 y nexo =/= 0 entonces cambian los archivos .form
   deltar=DBLE(rmax)/(DBLE(N)-1d0)
   r(0)=0.0d0
   DO i=1,N-1,1
	r(i)=r(i-1)+deltar
   END DO
   if(sump .ne. 0) then
        write(94,*) '!Volumen real 2'
        write(94,20) N,deltar,R(0)
        DO i=0,N-1,1
	   write(94,10) gaussv(REAL(r(i)),-VR,RR*(A**(1./3.)),AR,2,BETA2,BETA4,BETA6)
        END DO
   end if
   if(nexo .ne. 0) then
        write(94,*) '!Volumen real 2'
        write(94,20) N,deltar,R(0)
        DO i=0,N-1,1
	   write(94,10) gaussv(REAL(r(i)),-VR,RR*(A**(1./3.)),AR,2,BETA2,BETA4,BETA6)
        END DO
   endif
   if(sump .ne. 0) then
        write(94,*) '!Volumen dispersivo real e imaginario 2'
        write(94,20) N,deltar,R(0)
        DO i=0,N-1,1
	   write(94,10) gaussv(REAL(r(i)),-dv,drv*(A**(1./3.)),dav,2,BETA2,BETA4,BETA6)
           write(94,10) gaussv(REAL(r(i)),-W,RW*(A**(1./3.)),AW,2,BETA2,BETA4,BETA6)     
        END DO
   endif
   if(nexo .ne. 0) then
        write(94,*) '!Volumen dispersivo real e imaginario 2'
        write(94,20) N,deltar,R(0)
        DO i=0,N-1,1
	   write(94,10) gaussv(REAL(r(i)),-dv,drv*(A**(1./3.)),dav,2,BETA2,BETA4,BETA6)
           write(94,10) gaussv(REAL(r(i)),-W,RW*(A**(1./3.)),AW,2,BETA2,BETA4,BETA6)     
        END DO
   endif
   if(sump .ne. 0) then	
        write(94,*) '!Superficie dispersivo real e imaginario 2'
        write(94,20) N,deltar,R(0)
        DO i=0,N-1,1
	   write(94,10) gausss(REAL(r(i)),-VD,RVD*(A**(1./3.)),AVD,2,BETA2,BETA4,BETA6)
           write(94,10) gausss(REAL(r(i)),-WD,RD*(A**(1./3.)),AD,2,BETA2,BETA4,BETA6) 
        END DO
   endif
   if(nexo .ne. 0) then
        write(94,*) '!Superficie dispersivo real e imaginario 2'
        write(94,20) N,deltar,R(0)
        DO i=0,N-1,1
	   write(94,10) gausss(REAL(r(i)),-VD,RVD*(A**(1./3.)),AVD,2,BETA2,BETA4,BETA6)
           write(94,10) gausss(REAL(r(i)),-WD,RD*(A**(1./3.)),AD,2,BETA2,BETA4,BETA6) 
        END DO	
   endif
20 format(' ',i2,' ',f9.6,' ',f9.6)
10 format(' ',f9.6)
   close(94)
   return
  end
