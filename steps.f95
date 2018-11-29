    subroutine steps (mu,BETA2,BETA3,BETA2EFF,GAMMA2EFF,GAMMANAX, &
                     neg,neo,neb,nega,neax)
!   Los genera todos bien¡¡ // Los acoplamientos de G.S con G.S están como test de comprobación => Aprobado el 08/11/2018
    INTEGER neg,np,neo,neb,nega,neax,n
    DOUBLE PRECISION CG
    REAL BETA3,BETA2,BETA2EFF,GAMMA2EFF,GAMMANAX
    CHARACTER*10 banda
    INTEGER i,k,s,f,l,p,o,siz,mu ! mu=identificador de banda
    INTEGER IE(6),IEO(5),IEb(2),IEg(3),IEax(2),J(6),LAMDA(4),JMIN,JMAX,JAUX(40),JO(5),JB(2),JG(3),JAX(2)  
!   ne,np --> número de estados J's y número de multipolaridades distintas
!   i,k,s,l índices de bucles
!   IE,J,LAMDA --> Índices para steps,momento angular del estado, multipolaridad
!   JMIN,JMAX,JAUX --> momentos angulares para determinar los posibles acoplamientos J+LAMDA
2   format('&STEP ia=',i2,' ib=',i2,' k=',i1,' str=',1f9.4,'/')
    J(1:6)=(/0,2,4,6,8,10/) !Banda G.S
    JO(1:5)=(/1,3,5,7,9/) !Banda octupolar
    JB(1:2)=(/0,2/)  !Banda beta
    JG(1:3)=(/0,2,4/) ! Banda gamma
    JAX(1:2)=(/2,3/) !Banda anómala
    LAMDA(1:4)=(/2,4,6,3/) ! Multipolaridades
    !neg=size(J)
    !neo=size(JO)
    !neb=size(JB)
    !nega=size(JG)
    !neax=size(JAX)
    np=size(LAMDA)
! G.S -> G.S -----------------------------------------------------------------------------------------------------------------------
     if (mu==1) then
     DO i=1,neg
        IE(i)=i !Asignación  índice <------> J
     END DO
     DO i=1,neg,1
       DO k=1,3,1
         if (J(i).ge.LAMDA(k)) then !Elijo y almaceno en JAUX los posibles valores del acoplamiento de cada J con cada K
            JMIN=J(i)-LAMDA(k)
            JMAX=J(i)+ LAMDA(k)
            DO s=0,2*LAMDA(k)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*LAMDA(k))+1
         else
            JMIN=LAMDA(k)-J(i)
            JMAX=J(i)+LAMDA(k)
            DO s=0,2*J(i)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*J(i))+1
        endif
           DO l=1,siz,1 !Recorro todos los valores de JAUX
             DO p=1,neg,1 ! Recorro todos los valores de J
                 if (JAUX(l).eq.J(p)  .and. IE(i).lt.IE(p) .and. J(p).le.JMAX) then  ! Todas las transiciones menos reorientación       
                   !write(1,2) J(i),J(p),LAMDA(k) !para imprimir en steps momentos angular en vez de índices
                   !write(1,2) J(p),J(i),LAMDA(k)
                   call CLEBSCH(DBLE(J(i)),DBLE(LAMDA(k)),DBLE(J(p)),0d0,0d0,0d0,CG) !cálculo de CG*sqrt(2*IA + 1) = steps falta el beta_eff
                   write(1,2) IE(i),IE(p),LAMDA(k),CG !imprime los índices en steps
                   write(1,2) IE(p),IE(i),LAMDA(k),CG
                 endif
!                if (JAUX(l).eq.J(p) .and. IE(i).eq.IE(p) .and. J(p).le.JMAX) then
!                   write(1,2) J(i),J(p),LAMDA(k) ! Incluyo los términos de reorientación
!                endif 
             END DO
           END DO   
       END DO
     END DO 
     endif
! G.S --> Octupolar ------------------------------------------------------------------------------------------------------------------
     if (mu==2) then
        k=4 !eligo lambda=3
     DO i=1,neo
        IEO(i)=i+neg !Asignación  índice <------> J Octupolar
     END DO
     DO i=1,neg
        IE(i)=i !Asignación  índice <------> J G.S
     END DO
     DO i=1,neg,1
         if (J(i).ge.LAMDA(k)) then !Elijo y almaceno en JAUX los posibles valores del acoplamiento de cada J con cada K
            JMIN=J(i)-LAMDA(k)
            JMAX=J(i)+ LAMDA(k)
            DO s=0,2*LAMDA(k)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*LAMDA(k))+1
         else
            JMIN=LAMDA(k)-J(i)
            JMAX=J(i)+LAMDA(k)
            DO s=0,2*JO(i)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*JO(i))+1
        endif
           DO l=1,siz,1 !Recorro todos los valores de JAUX
             DO p=1,neo,1 ! Recorro todos los valores de J
                 if (JAUX(l).eq.JO(p)  .and. IE(i).lt.IEO(p) .and. JO(p).le.JMAX) then  ! Todas las transiciones menos reorientación       
                   !write(1,2) J(i),JO(p),LAMDA(k),1.0
                   !write(1,2) JO(p),J(i),LAMDA(k),1.0
                   call CLEBSCH(DBLE(J(i)),DBLE(LAMDA(k)),DBLE(JO(p)),0d0,0d0,0d0,CG)
                   write(1,2) IE(i),IEO(p),LAMDA(k),CG*BETA3
                   write(1,2) IEO(p),IE(i),LAMDA(k),CG*BETA3
                 endif  
             END DO
           END DO                                 
     END DO 
     endif 
! G.S ---> BETA ----------------------------------------------------------------------------------------------------------------------------
     if (mu==3) then
        k=1
     DO i=1,neb
        IEb(i)=i+neo+neg !Asignación  índice <------> J BETA
     END DO
     DO i=1,neg
        IE(i)=i !Asignación  índice <------> J G.S
     END DO
     DO i=1,neg,1
         if (J(i).ge.LAMDA(k)) then !Elijo y almaceno en JAUX los posibles valores del acoplamiento de cada J con cada K
            JMIN=J(i)-LAMDA(k)
            JMAX=J(i)+ LAMDA(k)
            DO s=0,2*LAMDA(k)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*LAMDA(k))+1
         else
            JMIN=LAMDA(k)-J(i)
            JMAX=J(i)+LAMDA(k)
            DO s=0,2*JB(i)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*JB(i))+1
        endif
           DO l=1,siz,1 !Recorro todos los valores de JAUX
             DO p=1,neb,1 ! Recorro todos los valores de J
                 if (JAUX(l).eq.JB(p)  .and. IE(i).lt.IEb(p) .and. JB(p).le.JMAX) then  ! Todas las transiciones menos reorientación       
                   !write(1,2) J(i),JB(p),LAMDA(k),1.0
                   !write(1,2) JB(p),J(i),LAMDA(k),1.0
                   call CLEBSCH(DBLE(J(i)),DBLE(LAMDA(k)),DBLE(JB(p)),0d0,0d0,0d0,CG)
                   write(1,2) IE(i),IEb(p),LAMDA(k),CG*BETA2EFF
                   write(1,2) IEb(p),IE(i),LAMDA(k),CG*BETA2EFF
                 endif  
             END DO
           END DO                                 
     END DO 
     endif
! G.S ---> GAMMA ----------------------------------------------------------------------------------------------------------------------------
     if (mu==4) then
        k=1
     DO i=1,nega
        IEg(i)=i+neo+neg+neb !Asignación  índice <------> J GAMMA
     END DO
     DO i=1,neg
        IE(i)=i !Asignación  índice <------> J G.S
     END DO
     DO i=1,neg,1
         if (J(i).ge.LAMDA(k)) then !Elijo y almaceno en JAUX los posibles valores del acoplamiento de cada J con cada K
            JMIN=J(i)-LAMDA(k)
            JMAX=J(i)+ LAMDA(k)
            DO s=0,2*LAMDA(k)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*LAMDA(k))+1
         else
            JMIN=LAMDA(k)-J(i)
            JMAX=J(i)+LAMDA(k)
            DO s=0,2*JG(i)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*JG(i))+1
        endif
           DO l=1,siz,1 !Recorro todos los valores de JAUX
             DO p=1,nega,1 ! Recorro todos los valores de J
                 if (JAUX(l).eq.JG(p)  .and. IE(i).lt.IEg(p) .and. JG(p).le.JMAX) then  ! Todas las transiciones menos reorientación       
                   !write(1,2) J(i),JG(p),LAMDA(k),1.0
                   !write(1,2) JG(p),J(i),LAMDA(k),1.0
                   call CLEBSCH(DBLE(J(i)),DBLE(LAMDA(k)),DBLE(JG(p)),0d0,0d0,0d0,CG)
                   write(1,2) IE(i),IEg(p),LAMDA(k),CG*GAMMA2EFF
                   write(1,2) IEg(p),IE(i),LAMDA(k),CG*GAMMA2EFF
                 endif  
             END DO
           END DO                                 
     END DO 
     endif
! G.S ---> GAMMA NO AXIAL --------------------------------------------------------------------------------------------------------------------
     if (mu==5) then
        k=1
     DO i=1,neax
        IEax(i)=i+neo+neg+neb+nega !Asignación  índice <------> J GAMMA
     END DO
     DO i=1,neg
        IE(i)=i !Asignación  índice <------> J G.S
     END DO
     DO i=1,neg,1
         if (J(i).ge.LAMDA(k)) then !Elijo y almaceno en JAUX los posibles valores del acoplamiento de cada J con cada K
            JMIN=J(i)-LAMDA(k)
            JMAX=J(i)+ LAMDA(k)
            DO s=0,2*LAMDA(k)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*LAMDA(k))+1
         else
            JMIN=LAMDA(k)-J(i)
            JMAX=J(i)+LAMDA(k)
            DO s=0,2*JAX(i)
              JAUX(s+1)=JMIN+s
            END DO
            siz=(2*JAX(i))+1
        endif
           DO l=1,siz,1 !Recorro todos los valores de JAUX
             DO p=1,neax,1 ! Recorro todos los valores de J
                 if (JAUX(l).eq.JAX(p)  .and. IE(i).lt.IEax(p) .and. JAX(p).le.JMAX) then  ! Todas las transiciones menos reorientación       
                   !write(1,2) J(i),JG(p),LAMDA(k),1.0
                   !write(1,2) JG(p),J(i),LAMDA(k),1.0
                   call CLEBSCH(DBLE(JAX(p)),DBLE(LAMDA(k)),DBLE(J(i)),-2d0,2d0,0d0,CG)
                   write(1,2) IE(i),IEax(p),LAMDA(k),sqrt(2.0)*((-1)**(DBLE(JAX(p))))*CG*GAMMANAX !Expresión exclusiva de los elementos para
                   write(1,2) IEax(p),IE(i),LAMDA(k),sqrt(2.0)*((-1)**(DBLE(JAX(p))))*CG*GAMMANAX !la banda NAX.
                 endif  
             END DO
           END DO                                 
     END DO 
     endif
     return
     end
!!!
   
