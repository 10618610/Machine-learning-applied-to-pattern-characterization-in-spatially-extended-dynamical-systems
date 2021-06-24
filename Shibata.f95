program d1
!sen(0)=sen(pi)=sen(2pi)=0 .. en(+ -pi/2)=1
integer::i,N,N1,t,tmax,seed,j,k,P,m,passos,delta,cont,epsi,ai
parameter(N=100)
real*8:: q,na,eps,alfa,a,b,d,e,sig,pi,x,xi=0.0d0,xf=1d0
real*8:: z(N,N), u1(N+2*n),u2(N+2*n),soma(4*n),na1(N+2*n),s1,s2,u3(N+2*n),f !funcao dopla precisao tem que declarar em real*8
real*8:: umin=0.d0,umax=0.d0,psi(N+2*n),phi(N+2*n),h,y(N+2*n),R(N+2*n),preal(N+2*n),pcompl(N+2*n),plato
real*8:: epsmax,epsmin,epspasso,amax,DETER,DET,Q1,q2,lambda,lamb(100*N),lya(100*N),ly(100*N),FindDet
external f,ran0,FindDet
external DET,der
open(1,file='shibat.dat')


tmax =5000
seed=2
b=0.1d0
epsmin=0.d0
epsmax = 1.0d0
P=1
amin=2.0
amax = 3.0d0
passos=500
!condição inicial
!do i=1,N
 ! u1(i)=ran0(seed)     !aleatoria
!end do
!a =2.5d0
! do while(a.LE.amax)
!eps=0.d0
! do while(eps.LE.epsmax)

apasso=(amax-amin)/float(passos-1)
  epspasso=(epsmax-epsmin)/dfloat(passos-1)
eps=0.1d0
  do 1000 ai=0,passos-1
  a=amin+apasso*float(ai)
   !do 1100 epsi=0,passos-1
   !eps=epsmin+epspasso*float(epsi)
q1 = 0.d0
do i=1,N
u1(i)=ran0(seed)
end do
            
   do i=1,n
   u1(i+N)=u1(i) ! u(i+n)=u(i) ...condicao de contorno periódica
   end do
   do i=1,N
   u1(i-N)=u1(i)  !u(i-n)=u(i) ... condicao de contorno periódica
   end do 

              do i=1,N
		do j=1,n
		   Z(i,j) = 0.d0
              end do
            end do

sig=2*real(p)
!P=1
 !do while(P.LE.Pmax)

  
      t=1
     do while(t.LE.tmax)

       do i=1,N 
       s1=0.d0       
         do j=i-P,i+P
         s2=0.d0
         s2=s2+f(u1(j),a)-f(u1(i),a)             
         s1=s1+s2      
         end do 
      u2(i)=f(u1(i),a)+(eps/sig)*s1
      end do
             q1=0.d0
            if(t.gt.tmax-100) then !i.gt.a = a menor que i
		q2 = 0.d0

             

                             do i=1,N
                               do k=1,N
                                 if(k.eq.i) then
                                   Z(i,k)=(1-eps)*(1+a-2*b*u1(k))
                                  else 
                                    Z(i,k)=0.d0
                                  end if
                              if(k+1.eq.i .or.k-1.eq.i) then
                             Z(i,K)=(eps/2)*(1+a-2*b*u1(k))
             
                                end if
           
          
                               end do
                               end do
                         DETER = FindDet(Z,N)
                    
                      DETER2=sqrt(deter*deter) 
                       lamb(t)=log(deter2)/float(n)

                  q1 = q1+lamb(t)
                    
                  
                       
                      
			
                      
                  
               end if
         
               
   do i=1,N
   u1(i)=u2(i)
   end do
   do i=1,N
   u1(n+i)=u1(i) ! u(-1)=u(1) ...
   end do
   do i=1,N
   u1(i-N)=u1(i)  !condicao de contorno 
   end do   
                                    
    

                      
                  !condicao de contorno 
             


                t=t+1
 end do 
write(1,*) a,q1
!1100 continue
! write(1,*)
1000 continue  
!write(1,*)  q1/100

!write(1,*) a,eps,q1
!eps=eps+0.0007518797
!end do
!write(1,*)
 !a=a+0.0012531328
!end do 
 
! P=P+1
!end do     
   close(1)

end program d1

function f(u,a)
IMPLICIT NONE
real*8::f,u,a,b=0.1d0
f=u+a*u-b*u*u
return
end
function der(u,a)
real::der,u,b=0.10
der=1+a-2*b*u
return
end
function ran0(idum)
integer:: idum,ia,im,iq,ir,mask
real:: ran0,am
parameter(ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,mask=123459876)
integer:: k
idum=ieor(idum,mask)
k=idum/iq
idum=ia*(idum-k*iq)-ir*k
if (idum.lt.0) idum=idum+im
ran0=am*idum
idum=ieor(idum,mask)
return
end
REAL*8 FUNCTION DET(A,V)
 ! A general purpose function written in FORTRAN77 to calculate determinant of a square matrix Passed parameters:
! A = the matrix
! N = dimension of the square matrix
 
IMPLICIT REAL*8 (A-H,O-Z)
INTEGER V
REAL*8 ELEM(V,V),A(V,V)
REAL*8 M, TEMP
INTEGER I, J, K, L
LOGICAL DETEXISTS
DO I=1,V
 DO J=1,V
 ELEM(I,J)=A(I,J)
 END DO
END DO
 DETEXISTS = .TRUE.
 L = 1
  !CONVERT TO UPPER TRIANGULAR FORM
DO K = 1, N-1
 IF (DABS(ELEM(K,K)).LE.1.0D-20) THEN
 DETEXISTS = .FALSE.
  DO I = K+1, V
  IF (ELEM(I,K).NE.0.0) THEN
  DO J = 1, V
 TEMP = ELEM(I,J)
 ELEM(I,J)= ELEM(K,J)
 ELEM(K,J) = TEMP
 END DO
 DETEXISTS = .TRUE.
 L=-L
 EXIT
 END IF
 END DO
 IF (DETEXISTS .EQV. .FALSE.) THEN
 GETDET = 0
 RETURN
 END IF
END IF
DO J = K+1, V
 M = ELEM(J,K)/ELEM(K,K)
 DO I = K+1, V
 ELEM(J,I) = ELEM(J,I) - M*ELEM(K,I)
 END DO
 END DO
END DO
!CALCULATE DETERMINANT BY FINDING PRODUCT OF DIAGONAL ELEMENTS
DET = L
DO I = 1, V
 DET = DET * ELEM(I,I)
END DO
END

REAL*8 FUNCTION FindDet(matrix, n)
    IMPLICIT NONE
    REAL*8, DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL*8 :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
   
END FUNCTION FindDet

