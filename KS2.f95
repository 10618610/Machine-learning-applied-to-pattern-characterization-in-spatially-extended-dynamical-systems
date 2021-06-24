program d1
!esse programa variamos os parametros de controle para cada tempo dado
integer::i,N,t,trans,seed,k,l,mmax,m
integer:: passos,ai,epsi
parameter(N=100)
parameter(Np=N)
real::u1(N),u2(N),soma(N),x(N),a,lambda(N),wi(np),wr(np),w(n),y(n),v(n),q1,q2,KS,h(N)
real:: S(N,N),D(N,N),Z(N,N)
real:: amax,amin,epsmax,epsmin,apasso,epspasso
CHARACTER(len=15) :: FN
external f,ran0,der,balanc,elmhes,hqr, GRAMSCHMIDT, tred2
open(1,file='Lya_max.dat')

mmax=10 !numero de jacobianas multiplicadas
trans=5000
!amin=2.45
amax=3.0
!epsmin=0.0
!epsmax=0.300
eps=0.1d0
seed=2
passos=400

!condição inicial

  apasso=(amax-amin)/float(passos-1)
  epspasso=(epsmax-epsmin)/dfloat(passos-1)

 ! do 1000 ai=0,passos-1
 ! a=amin+apasso*float(ai)
  ! do 1100 epsi=0,passos-1
  ! eps=epsmin+epspasso*float(epsi)
   !Calculo do transiente
a=2.0d0
DO while(a .LE. amax) 
do i=1,N
   u1(i)=ran0(seed)
 end do
q1=0.0
  do 1001 t=1,trans
    soma(1)=0.0d0
    soma(N)=0.0d0
    soma(1)=soma(1)+f(u1(2),a)+f(u1(N),a)
    u2(1)=(1-eps)*f(u1(1),a)+(eps/2)*soma(1)
    
           do i=2,N-1
           soma(i)=0.0
           soma(i)=soma(i)+f(u1(i-1),a)+f(u1(i+1),a)
           u2(i)=(1-eps)*f(u1(i),a)+(eps/2)*soma(i)
           end do
        soma(N)=soma(N)+f(u1(N-1),a)+f(u1(1),a)
        u2(N)=(1-eps)*f(u1(N),a)+(eps/2)*soma(N)
            
                
             do i=1,N
             u1(i)=u2(i)   !valor de U1(i) depois do transiente
             end do

 1001 continue

  !calculo do matriz jacobiana    
       do k=1,N
        do l=1,N
        S(k,l)=0.0
        end do
       end do

       do k=1,N
       S(k,k)=1.0      !definindo a matriz identidade S(i,j)=I
       end do

       m=0       
  do while(m.LE.mmax)
      do k=1,N
        do l=1,N
         D(k,l)=0.0    !zerando todos elementos da matriz D(i,j)
        end do
       end do
  soma(1)=0.0d0
  soma(N)=0.0d0
  soma(1)=soma(1)+f(u1(2),a)+f(u1(N),a)
  u2(1)=(1-eps)*f(u1(1),a)+(eps/2)*soma(1)
    
    do i=2,N-1
    soma(i)=0.0d0
    soma(i)=soma(i)+f(u1(i-1),a)+f(u1(i+1),a)
    u2(i)=(1-eps)*f(u1(i),a)+(eps/2)*soma(i)
    end do
     soma(N)=soma(N)+f(u1(N-1),a)+f(u1(1),a)
     u2(N)=(1-eps)*f(u1(N),a)+(eps/2)*soma(N)
            
                
       do i=1,N
       u1(i)=u2(i)
       end do     
              
        do i=1,N
         do k=1,N
          if(k.eq.i) then
           Z(i,k)=(1-eps)*der(u1(k),a)
          else 
            Z(i,k)=0
          end if
            if(k+1.eq.i .or.k-1.eq.i) then
            Z(i,K)=(eps/2)*der(u1(k),a)
            
          end if
           
          
         end do
        end do

        do k=1,N
         do l=1,N
          do i=1,N
         D(k,l)=D(k,l)+S(k,i)*Z(i,l)
          end do
         end do
        end do


           do k=1,N
            do l=1,N
            S(k,l)=D(k,l)
            end do
           end do
         !D(K,L) FORNECE O RESULTADO DO PRODUTO DAS JACOBIANAS 
  m=m+1
end do

  !CALCULO DOS AUTOVALORES   
   !call GRAMSCHMIDT(s,N,M,NMAX)    
   call balanc(S,N,np)
   call elmhes(S,N,np)
   call hqr(s,N,np,wr,wi)
   !call tred2(s,n,np,y,v)
! CALCULO DO EXPOENTE DE LYAPUNOV
     
       do i=1,N
       w(i)=sqrt(wr(i)*wr(i)+wi(i)*wi(i))
       end do
          do i=1,n
        lambda(i)=(log(w(i)))/(mmax)
         end do
  !ORDENAR lambda(i) EM ORDEM DECRESCENTE 

            do l=1,N
                do i=1,N
                 if(lambda(i+1)>=lambda(i)) then
                  temp=lambda(i)
                  lambda(i)=lambda(i+1)
                  lambda(i+1)=temp
                 end if
                end do
             end do

    !imprimindo as valores de lambda
            !do i=1,n
   !write(1,*)i, lambda(i)
         
        !end do
!calculando KS


  do i=1,N
   if(lambda(i) >= 0.0) then
    q1=q1+lambda(i)
    KS=q1/N
   end if
end do
         
write(1,*) a, lambda(1)
a = a+0.005
end do

!1100 continue
 !write(1,*)
!1000 continue             
 close(1)
stop
end program d1


function f(u,a)
real::f,u,b=0.10
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
!Dado uma matriz n por n, armazenada em uma matriz de dimensões físicas np por np, esta rotina o substitui por uma matriz !equilibrada com autovalores idênticos. Uma matriz simétrica já está equilibrada e não é descrita por este procedimento. O  !parâmetro RADIX deve ser a base de ponto flutuante da máquina

SUBROUTINE balanc(s,n,np)
INTEGER n,np
REAL S(np,np),RADIX,SQRDX
PARAMETER (RADIX=2.0,SQRDX=RADIX**2)
INTEGER i,j,last
REAL d,f,g,r,W
1   continue
last=1
do 14 i=1,n !Calcule as normas de linha e coluna
r=0
d=0
do 11 j=1,n
 if(j.ne.i)then !.ne. eh diferente a.ne.b = a#b
 d=d+abs(s(j,i))
 r=r+abs(s(i,j))
 end if
11 continue

if(d.ne. 0.0 .and.r.ne. 0.0) then   ! se ambos são diferente de zero
g=r/RADIX
f=1.0
W=d+r
2 if(d.lt.g) then   !a.lt.b = a<b
f=f*RADIX
d=d*SQRDX
go to 2
end if
g=r*RADIX
3 if(d.gt.g)then !a.gt.b = a>b
f=f/RADIX
d=d/SQRDX
go to 3
end if
 if((d+r)/f .lt. 0.95*W)then
 last=0
 g=1.0/f
 do 12 j=1,n     !Apply  similarity  transformation.
 s(i,j)=s(i,j)*g
 12 continue
do 13 j=1,n
s(j,i)=s(j,i)*f
 13 continue
end if
end if
14 continue
if(last.eq. 0)goto 1
return
END
!Redução para a forma de Hessenberg pelo método de eliminação. A matriz real, não-simétrica, n por n a, armazenada em uma matriz de dimensões físicas np por np, é substituído por uma matriz superior de Hessenberg com autovalores idênticos. Recomendado, mas não obrigatório, é que esta rotina seja precedida por balanc. Na saída, a matriz de Hessenberg está nos elementos a (i, j) com i<= j + 1. Os elementos com i> j + 1 devem ser pensados como zero, mas são retornados com valores aleatórios
SUBROUTINE elmhes(s,n,np)
INTEGER n,np
REAL s(np,np)
INTEGER i,j,m
REAL x,y
do 17 m=2,n-1                !m é chamado r+1 no texto
x=0.0
i=m
do 11 j=m,n                  !encontrar o pivo
 if(abs(s(j,m-1)).gt.abs(x))then
 x=s(j,m-1)
 i=j
 end if
11 continue

 if(i.ne.m)then               !Intercâmbio de linhas e colunas.
  do 12 j=m-1,n
  y=s(i,j)
  s(i,j)=s(m,j)
  s(m,j)=y
12 continue
   do 13 j=1,n
   y=s(j,i)
   s(j,i)=s(j,m)
   s(j,m)=y
13 continue
   end if
    if(x.ne. 0.0)then   !Execute a eliminação.

      do 16 i=m+1,n
      y=s(i,m-1)
       if(y.ne. 0.0)then
       y=y/x
       s(i,m-1)=y
         do 14 j=m,n
         s(i,j)=s(i,j)-y*s(m,j)
         14 continue
           do 15 j=1,n
           s(j,m)=s(j,m)+y*s(j,i)
           15 continue
         endif
     16 continue
endif
17 continue
return
END

!Encontra todos os autovalores de uma matriz de Hessenberg superior n por n e que é armazenada em uma matriz np por np. Na !entrada, um pode ser exatamente como saída de elmhes §11.5; na saída é destruída. As partes reais e imaginárias dos autovalores são retornadas em wr e wi, respectivamente.
SUBROUTINE hqr(s,n,np,wr,wi)
INTEGER n,np
REAL s(np,np),wi(np),wr(np)

INTEGER i,its,j,k,l,m,nn
REAL anorm,p,q,r,a,t,u,v,w,x,y,z         
anorm=0.
                  !Compute a norma da matriz para possível uso na localização de um único elemento subdiagonal pequeno.
  do 12 i=1,n
  do 11 j=max(i-1,1),n
  anorm=anorm+abs(s(i,j))
 11 continue
 12 continue
 nn=n
 t=0.0

1  if(nn.ge.1)then        !Comece a procurar o próximo autovalor.  
                                   
   its=0
2   do 13 l=nn,2,-1             !Comece a iteração: procure um único elemento 
                             !pequeno subdiagonal.                              
    a=abs(s(l-1,l-1))+abs(s(l,l))

     if(a.eq.0.)a=anorm
     if(abs(s(l,l-1))+a.eq.a)goto 3
    13 continue
    l=1
3   x=s(nn,nn)
     if(l.eq.nn)then
                                             !Uma raiz encontrada.
    wr(nn)=x+t
    wi(nn)=0.0
    nn=nn-1
    else
    y=s(nn-1,nn-1)
    w=s(nn,nn-1)*s(nn-1,nn)
    if(l.eq.nn-1)then                                    ! Duas raizes encontradas...
    p=0.5*(y-x)
    q=p**2+w
    z=sqrt(abs(q))
   x=x+t
  if(q.ge.0.)then                                       !...a real pair.
  z=p+sign(z,p)
  wr(nn)=x+z
  wr(nn-1)=wr(nn)
  if(z.ne.0.)wr(nn)=x-w/z
  wi(nn)=0.
  wi(nn-1)=0.
   else                                                    !...a complex pair.
  wr(nn)=x+p
  wr(nn-1)=wr(nn)
  wi(nn)=z
  wi(nn-1)=-z
  endif
  nn=nn-2
  else                                                       !No roots found. Continue iteration.
  if(its.eq.300)pause 'too many iterations in hqr'
  if(its.eq.10.or.its.eq.20)then                              !Frm exceptional shift.
  t=t+x
  do 14 i=1,nn
  s(i,i)=s(i,i)-x
 14 continue
 a=abs(s(nn,nn-1))+abs(s(nn-1,nn-2))
 x=0.75*a
 y=x
 w=-0.4375*a**2
 end if
 its=its+1
 do 15 m=nn-2,l,-1            !Form shift and then look for 2 consecu-
 z=s(m,m)                      !tive small subdiagonal elements.
 r=x-z
 a=y-z
 p=(r*a-w)/s(m+1,m)+s(m,m+1)  
 q=s(m+1,m+1)-z-r-a
 r=s(m+2,m+1)
 a=abs(p)+abs(q)+abs(r)         !Scale to prevent overflow or underflow.
 p=p/a
 q=q/a
 r=r/a
 if(m.eq.l)goto 4

  u=abs(s(m,m-1))*(abs(q)+abs(r))
  v=abs(p)*(abs(s(m-1,m-1))+abs(z)+abs(s(m+1,m+1)))
 if((u+v).eq.v)goto 4                                    
15 continue

4  do 16 i=m+2,nn
 s(i,i-2)=0.
 if (i.ne.m+2) s(i,i-3)=0.
16 continue
 do 19 k=m,nn-1                              !Double QR step on rows l to nn and
 if(k.ne.m)then                                !columns m to nn.
 p=s(k,k-1)                                     !Begin setup of Householder vector.
 q=s(k+1,k-1)
 r=0.0
 if(k.ne.nn-1)r=s(k+2,k-1)
 x=abs(p)+abs(q)+abs(r)
 if(x.ne.0.)then
 p=p/x                                             !Scale to prevent overflow or underflow.
 q=q/x
 r=r/x
 endif
 endif
 a=sign(sqrt(p**2+q**2+r**2),p)
 if(a.ne.0.)then
 if(k.eq.m)then
 if(l.ne.m)s(k,k-1)=-s(k,k-1)
 else
 s(k,k-1)=-a*x
 endif
 p=p+a           
 x=p/a
 y=q/a
 z=r/a
 q=q/p
 r=r/p
 do 17 j=k,nn            !Row modification.
 p=s(k,j)+q*s(k+1,j)
  if(k.ne.nn-1)then
  p=p+r*s(k+2,j)
  s(k+2,j)=s(k+2,j)-p*z
  endif
 s(k+1,j)=s(k+1,j)-p*y
 s(k,j)=s(k,j)-p*x
 17 continue
 do 18 i=l,min(nn,k+3)       !Column modification.
 p=x*s(i,k)+y*s(i,k+1)
 if(k.ne.nn-1)then
 p=p+z*s(i,k+2)
 s(i,k+2)=s(i,k+2)-p*r
 endif
 s(i,k+1)=s(i,k+1)-p*q
 s(i,k)=s(i,k)-p
 18 continue
 endif
 19 continue
 goto 2                      !...for next iteration on current eigenvalue.
 endif
 endif
 goto 1                      !...for next eigenvalue.
 endif
 return
 END

SUBROUTINE GRAMSCHMIDT(s,N,M,NMAX)
IMPLICIT NONE
INTEGER NMAX,N,M
DOUBLE PRECISION s(NMAX,*)
INTEGER I,J,L
DOUBLE PRECISION PI
DO J=1,M
DO L=1,J-1


PI = 0.0D0
DO I=1,N
PI = PI + s(I,L)*s(I,J)
END DO

DO I=1,N
s(I,J) = s(I,J)-PI*s(I,L)
END DO
END DO

PI = 0.0D0
DO I=1,N
PI = PI + s(I,J)**2
END DO
PI = SQRT(PI)

DO I=1,N
s(I,J) = s(I,J)/PI
END DO
END DO
RETURN
END

SUBROUTINE tred2(s,n,np,y,v)
INTEGER n,np
REAL s(np,np),y(np),v(np)
INTEGER i,j,k,l
REAL f,g,h,hh,scale
do 18 i=n,2,-1
l=i-1
h=0.
scale=0.
if(l.gt.1)then
do 11 k=1,l
scale=scale+abs(s(i,k))
 11 continue
if(scale.eq.0.)then
v(i)=s(i,l)
else
do 12 k=1,l
s(i,k)=s(i,k)/scale
h=h+s(i,k)**2
 12 continue
f=s(i,l)

g=-sign(sqrt(h),f)
v(i)=scale*g
h=h-f*g
s(i,l)=f-g

f=0.
do 15 j=1,l
s(j,i)=s(i,j)/h
g=0.
do 13 k=1,j
g=g+s(j,k)*s(i,k)
 13 continue
do 14 k=j+1,l
g=g+s(k,j)*s(i,k)
 14  continue
v(j)=g/h

f=f+v(j)*s(i,j)

 15 continue
hh=f/(h+h)

do 17 j=1,l

f=s(i,j)
g=v(j)-hh*f
v(j)=g
do 16 k=1,j

s(j,k)=s(j,k)-f*v(k)-g*s(i,k)
 16  continue
 17  continue
endif
else
v(i)=s(i,l)
endif
y(i)=h
 18 continue

y(1)=0.0
v(1)=0.0
do 24 i=1,n


l=i-1
if(y(i).ne.0.)then

do 22 j=1,l
g=0.0
do 19 k=1,l

g=g+s(i,k)*s(k,j)
 19 continue
do 21 k=1,l
s(k,j)=s(k,j)-g*s(k,i)
21 continue 
22 continue
endif

y(i)=s(i,i)

s(i,i)=1.

do 23 j=1,l

s(i,j)=0.
s(j,i)=0.
 23 continue
 24 continue
return
END

