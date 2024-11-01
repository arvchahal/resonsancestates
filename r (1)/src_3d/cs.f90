use linalg
implicit none
double precision,allocatable  :: A(:,:),AI(:,:),AA(:,:),te(:,:),ve(:,:), &
  he(:,:),oe(:,:),ci(:,:),cr(:,:),h(:,:),om0(:,:),sm(:,:),pm(:,:),hm0(:,:),ev(:,:),eva(:), &
  hmm(:,:),omm(:,:),hmmm(:,:),wm(:,:)
complex*16,allocatable        :: hm(:,:),om(:,:),c(:,:),v(:,:),e(:)
integer                       :: i,j,n,k,m,nm
double precision              :: theta,th
complex*16,parameter          :: zi=(0.d0,1.d0)
real*4                        :: t1,t2
complex*16                    :: p1,p2,p3,p4
complex*16,dimension(:,:),allocatable   :: VL,VR
complex*16,dimension(:),allocatable     :: alpha,beta,work
real*8,dimension(:),allocatable         :: RWORK
CHARACTER                               :: JOBVL, JOBVR
INTEGER                                 :: INFO,LWORK
integer,dimension(:),allocatable        :: indx
real*8,dimension(:),allocatable         :: ener
  
open(1,file='cs.dat')
open(2,file='gr.dat')
read(1,*)nm,n
allocate(A(N,N),AI(N,N),AA(N,N),hm(N,N),om(N,N),c(n,n),v(n,n),e(n))
allocate(ve(N,N),te(N,N),he(n,n),oe(n,n))
allocate(ci(N,N),cr(N,N),h(n,n),wm(n,n))
allocate(hm0(N,N),om0(N,N),sm(N,N),pm(n,n))
allocate(hmm(N,N),omm(N,N),hmmm(n,n))
allocate(ev(N,N),eva(n))
allocate(ener(n),indx(n))
  
JOBVL='V'
JOBVR='V'
LWORK=8*N
allocate(RWORK(8*N),ALPHA(N),BETA(N),VL(N,N),VR(N,N),WORK(LWORK))

do k=0,nm
  read(1,*)theta
  write(6,*)k
  do i=1,n
    do j=1,n
      read(1,*)hmmm(j,i),hm0(j,i),te(j,i),ve(j,i),sm(j,i),om0(j,i),pm(j,i),wm(i,j)
    end do
  end do
  if(k==0) then
    omm=om0
    call diag(hmmm,omm,n,eva,ev)
    write(6,*)eva(1)
  endif  

  om=om0
! rotations
  p1=exp(-zi*theta)
  p2=exp(-2.d0*zi*theta)
  p3=exp(zi*theta)
  p4=exp(2.d0*zi*theta)
! te : kinetic energy
! ve : Coulomb
! pm : polarization
! sm : self polarization
! wm : omega term
  do i=1,n
    do j=1,n
      hm(j,i)=te(j,i)*p2+ve(j,i)*p1+sm(j,i)*p4+pm(j,i)*p3+wm(i,j)
      write(99,*)hm(j,i),j,i
      write(99,*)hmmm(j,i)
    end do
  end do
  hmm=hm
  omm=om
  call diag(hmm,omm,n,eva,ev)
  write(6,*)eva(1)
    
  call ZGGEV(JOBVL,JOBVR,N,HM,N,OM,N,ALPHA,BETA,VL,N,VR,N,WORK,LWORK,RWORK,INFO)
  do i=1,n
    if(beta(i).ne.(0.d0,0.d0)) then
      ener(i)=real(alpha(i)/beta(i))
    else
      ener(i)=1000000.d0
    endif
    write(10,*)i,ener(i)
    write(10,*)i,alpha(i)/beta(i),beta(i)
  end do
  call indexx(n,ener,indx)  
  do i=1,100
    write(2,*)real(alpha(indx(i))/beta(indx(i))),imag(alpha(indx(i))/beta(indx(i)))
  end do
  write(2,*)
end do
close(2)

end



      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL*8 a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
