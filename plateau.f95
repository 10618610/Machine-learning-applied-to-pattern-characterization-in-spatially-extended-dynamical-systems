program d1
!sen(0)=sen(pi)=sen(2pi)=0 .. sen(+ -pi/2)=1
integer::i,N,N1,t,tmax,seed,j,k,P,m,passos,delta,cont
parameter(N=100)
integer:: ai,epsi
real*8:: q,na,eps,alfa,a,b,d,e,sig,pi,x,xi=0.0d0,xf=1d0
real*8:: u1(N+2*n),u2(N+2*n),soma(4*n),na1(N+2*n),s1,s2,u3(N+2*n),f !funcao dopla precisao tem que declarar em real*8
real*8:: umin=0.d0,umax=0.d0,psi(N+2*n),phi(N+2*n),h,y(N+2*n),R(N+2*n),preal(N+2*n),pcompl(N+2*n),plato
real*8:: epsmax,epsmin,epspasso
external f,ran0
open(1,file='plateau.dat')

tmax=5000
seed=2
P=1
amin = 2.0d0
amax=3.0
epsmin=0.d0
epsmax=0.300
passos=500
!condição inicial
eps = 0.1d0



 apasso=(amax-amin)/float(passos-1)
  epspasso=(epsmax-epsmin)/dfloat(passos-1)

  do 1000 ai=0,passos-1
  a=amin+apasso*float(ai)
  ! do 1100 epsi=0,passos-1
  ! eps=epsmin+epspasso*float(epsi)
do i=1,N
  u1(i)=ran0(seed)     !aleatoria
end do          
   do i=1,n
   u1(i+N)=u1(i) ! u(i+n)=u(i) ...condicao de contorno periódica
   end do
   do i=1,N
   u1(i-N)=u1(i)  !u(i-n)=u(i) ... condicao de contorno periódica
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
            if(t.gt.tmax-1) then !i.gt.a = a menor que i

              umin=u1(1)
                do k=2,n
                  if(u1(1).gt.u1(k)) then
                   umin=u1(k)
                   u1(1)=u1(k)
                   end if
                 end do
                umin=u1(1)
                      umax=u1(1)
                 do k=2,n
                   if(u1(1).lt.u1(k)) then !.lt. <         .gt. >
                    umax=u1(k)
                    u1(1)=u1(k)
                    end if
                 end do

                       do i=2,n
                        y(i)=(2*u1(i)-umax-umin)/(umax-umin)
                        if(umin.eq.umax) y(i)=0
                       
                         psi(i)=asin(y(i))
                        end do
                          do i=1,N
                          psi(n+i)=psi(i) 
                          end do
                          do i=1,N
                          psi(i-N)=psi(i)  !condicao de contorno 
                          end do  
                           delta=1
                           e=2*delta
                           cont=0.d0
             
                             do i=1,N
                               do j=i-delta,i+delta
                               preal(j)=0.d0
                               pcompl(j)=0.d0
                               preal(j)=preal(i)+cos(psi(j))
                               pcompl(j)=pcompl(i)+sin(psi(j))  
                             R(i)=(1/e)*sqrt(preal(j)*preal(j)+pcompl(j)*pcompl(j))
                               end do
                                if(R(i)>0.999d0) then
                                 cont=cont+1.d0
                                end if
                             end do
                             
                             plato=cont
                             plato=plato/N  


                  

                  
                  
                  
                      write(1,*) a, plato  
                  
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
!1100 continue
 !write(1,*)
1000 continue             
 close(1)
stop
end program d1    
   



function f(u,a)
IMPLICIT NONE
real*8::f,u,a,b=0.1d0
f=u+a*u-b*u*u
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
