
!  e_best ---> z_best

MODULE CORRELATED_GAUSSIANS
  implicit none
  real*8,parameter                    :: h2m=0.5d0
  integer                             :: N_particle,N_basis_final,N_pairs,N_basis_start
  real*8,allocatable                  :: Alpha(:,:),S_vector(:,:,:)
  integer                             :: Number_of_atoms
  real*8,allocatable                  :: Position(:,:),Nuclear_Charge(:),Charge(:),Mass(:)
  real*8,parameter                    :: pi=3.1415926535897932384d0
  integer,allocatable                 :: ipr(:,:),perm(:),spins(:,:),center(:),photons(:),parity(:)
  real*8,allocatable                  :: Spin_coef(:),spin_ov(:),Lambda(:,:),pkind(:)
  integer                             :: Nperm,N_spin_conf
  integer,parameter                   :: N_particle_max=12
  real*8,allocatable                  :: hm(:,:),om(:,:),hr(:),or(:),or0(:),pr0(:),pm(:,:), &
&   hm0(:,:),om0(:,:),hr0(:),sr0(:),sm(:,:),tr0(:),vr0(:),tm0(:,:),vm0(:,:),wr0(:),wm(:,:)
  real*8,allocatable                  :: op(:),hp(:)
! norm of the basis function
  double precision,allocatable        :: norm(:)
! eigenvector and eigenvalue
  double precision,allocatable        :: eigvec(:,:)
  double precision,allocatable        :: eigval(:)
  double precision,allocatable        :: eigvec0(:,:)
  double precision,allocatable        :: eigval0(:)
  double precision, allocatable       :: qqv(:),eqv(:)
  integer                             :: n_current

  double precision                    :: E_ion_ion
!
  integer                             :: nt1,nt2
  double precision                    :: bm1,bm2,bms(3),bmin
  integer                             :: iseed,ico,n_refine
  real*8,allocatable                  :: ac(:,:),sc(:,:,:)
  real*8,allocatable                  :: hc(:,:)
  integer,allocatable                 :: pc(:),pcharge(:)
  integer                             :: n_opt
  real*8                              :: e_low,cs
  real*8,allocatable                  :: s_trial(:,:,:)
  integer                             :: N_s_trial
!  real*8,parameter                    :: lambda_pol(2)=(/2.98374508d0,2.98374508d0/)
  real*8                              :: lambda_pol(2)
  real*8                              :: lambda_norm(2)
  real*8                              :: omega,omega0
!  real*8                              :: omega=0.124991895d0
!  real*8,parameter                    :: lambda_pol(2)=(/1.,1./)
!  real*8,parameter                    :: omega=1.
  logical                             :: density,full_matrix,cm,exc
  integer                             :: nexc,excmode
! starting photon number
  integer                             :: photon_min
! dimension of basis for a given photon
  integer                             :: photon_dim
  real*8 :: td,tm
  real*8,allocatable                  :: umat(:,:),cm_mat(:,:)
! one body potential type
  integer                             :: onebody




CONTAINS

subroutine initialize
  implicit none
  integer    :: i,k,j,id
  real*8     :: rii,xm,dd,xk


open(1,file='cg.inp')
read(1,*)
read(1,*)N_particle
allocate(Charge(N_particle),mass(N_particle),pkind(n_particle))
read(1,*)
read(1,*)(charge(i),i=1,N_particle)
read(1,*)(pkind(i),i=1,N_particle)
read(1,*)
read(1,*)(mass(i),i=1,N_particle)
read(1,*)
read(1,*)N_spin_conf
allocate(Spin_coef(N_spin_conf),spins(N_particle,N_spin_conf))
read(1,*)
do k=1,N_spin_conf
  read(1,*)Spin_coef(k),(Spins(i,k),i=1,N_particle)
end do
read(1,*)
read(1,*)N_basis_final
read(1,*)
read(1,*)Number_of_atoms
allocate(Position(3,Number_of_atoms),Nuclear_Charge(Number_of_atoms))
read(1,*)
do i=1,Number_of_atoms
  read(1,*)(Position(k,i),k=1,3),Nuclear_Charge(i)
end do
write(6,*)'?1'
read(1,*)
read(1,*)nt1,nt2
read(1,*)
read(1,*)bm1,bm2,(bms(i),i=1,3),bmin
read(1,*)
read(1,*)iseed
read(1,*)
write(6,*)'?1.1'
read(1,*)ico,n_refine
read(1,*)
read(1,*)density
write(6,*)'?1.2'
read(1,*)
read(1,*)omega
write(6,*)'?1.3'
read(1,*)
read(1,*)dd,lambda_pol(1),lambda_pol(2)
lambda_pol=lambda_pol*dd


!lambda_pol(1)=2.98374508d0*dd
!lambda_pol(2)=2.98374508d0*dd
! normalized direction
dd=lambda_pol(1)**2+lambda_pol(2)**2
lambda_norm=0.d0
if(dd>0.d0) lambda_norm=lambda_pol/sqrt(dd)
read(1,*)
read(1,*)photon_min,photon_dim
read(1,*)
read(1,*)cm
read(1,*)
write(6,*)'?2'
read(1,*)onebody
read(1,*)
read(1,*)exc,nexc,excmode
read(1,*)
read(1,*)omega0


close(1)
!write(6,*)'start'
N_pairs=(N_particle*(N_particle-1))/2
allocate(pcharge(N_pairs))
j=0
  do i=1,N_particle
    do k=i+1,N_particle
      j=j+1
      pcharge(j)=charge(i)*charge(k)
    end do
  end do

  allocate(hm(N_basis_final,N_basis_final),om(N_basis_final,N_basis_final))
  allocate(wm(N_basis_final,N_basis_final))
  allocate(eigvec(N_basis_final,N_basis_final),eigval(N_basis_final))

allocate(Alpha(N_particle+N_pairs,N_basis_final),S_vector(3,N_particle,N_basis_final))
allocate(pm(N_basis_final,N_basis_final))
allocate(sm(N_basis_final,N_basis_final))
allocate(hm0(N_basis_final,N_basis_final))
allocate(om0(N_basis_final,N_basis_final))
allocate(tm0(N_basis_final,N_basis_final))
allocate(vm0(N_basis_final,N_basis_final))
allocate(hr(N_basis_final),norm(N_basis_final))
allocate(eigvec0(N_basis_final,N_basis_final),eigval0(N_basis_final))
allocate(ac(N_particle+N_pairs,N_basis_final),sc(3,N_particle,N_basis_final))
allocate(hc(N_basis_final,N_basis_final))
allocate(photons(N_basis_final))
allocate(parity(N_basis_final))
allocate(pc(N_basis_final))
allocate(hr0(N_basis_final))
allocate(wr0(N_basis_final))
allocate(or(N_basis_final))
allocate(tr0(N_basis_final))
allocate(vr0(N_basis_final))
allocate(sr0(N_basis_final))
allocate(or0(N_basis_final))
allocate(pr0(N_basis_final))
allocate(qqv(N_basis_final))
allocate(eqv(N_basis_final))

E_ion_ion=0.d0
do i=1,Number_of_atoms
  do k=i+1,Number_of_atoms
    rii=sqrt((Position(1,i)-Position(1,k))**2+(Position(2,i)-Position(2,k))**2+ &
&     (Position(3,i)-Position(3,k))**2)
    E_ion_ion=E_ion_ion+Nuclear_Charge(i)*Nuclear_Charge(k)/rii
  end do
end do

call antisymmetry

!  lambda matrix for kinetic energy
  allocate(lambda(N_particle,N_particle))
  allocate(umat(N_particle,N_particle))
  allocate(cm_mat(N_particle,N_particle))


  umat=0.d0
  xm=sum(mass)
  do i=1,N_particle
    do j=1,N_particle
      dd=0.d0
      if(i==j) dd=1.d0
      umat(i,j)=dd-mass(j)/xm
    end do
  end do

  cm_mat=0.d0
  do i=1,n_particle
    do j=1,n_particle
      cm_mat(i,j)=mass(i)*mass(j)/xm**2
    end do
  end do
  

  if(cm) then
!    xk=-1.d0
    lambda=0.d0
    do i=1,n_particle
      do j=1,n_particle
        do k=1,n_particle
          lambda(i,j)=lambda(i,j)+umat(i,k)*umat(j,k)/mass(k)
        end do
      end do
    end do  
  else
    xk=0.d0
    do i=1,n_particle
      do j=1,n_particle
        lambda(i,j)=xk/dfloat(n_particle)/mass(i)
      end do
      lambda(i,i)=(1.d0+xk/dfloat(n_particle))/mass(i)
    end do  
  endif
  if(N_particle==1) lambda=1.d0/mass(1)



end subroutine initialize

subroutine antisymmetry
  implicit none
  integer   :: i,fac,ip,ia(N_particle_max),k,i1,i2,ov,ok
  real*8    :: su
  fac=1
  do i=1,N_particle
    fac=fac*i
  end do
  allocate(ipr(N_particle_max,fac))
  allocate(perm(fac),spin_ov(fac))

  do i=1,fac
    call permut(i,N_particle,ia)
    call pape(ia,N_particle,ip)
    perm(i)=ip
    ipr(:,i)=ia(:)
  end do
  nperm=fac

  do i=1,nperm
    if(mod(i,1000)==0) write(6,*)i
    su=0.d0
    ok=1
    do k=1,N_particle
      if(pkind(k).ne.pkind(ipr(k,i))) ok=0
    end do
    if(ok==1) then
      do i1=1,N_spin_conf
        do i2=1,N_spin_conf
          ov=1
          do k=1,N_particle
            if(spins(k,i1).ne.spins(ipr(k,i),i2)) ov=0
            !          if(pkind(k).ne.pkind(ipr(k,i))) ov=0
          end do
          su=su+Spin_coef(i1)*Spin_coef(i2)*ov
        end do
      end do
      spin_ov(i)=perm(i)*su
    endif
  end do
  
end subroutine antisymmetry


subroutine partial_diag(m,n,ener)
  use linalg
  implicit none
  integer :: i,j,m,n,fail
  double precision,allocatable   :: h_m(:,:),o_m(:,:),eigval_m(:),eigvec_m(:,:)
  double precision               :: ener

  allocate(h_m(m,m),o_m(m,m),eigval_m(m),eigvec_m(m,m))
  do i=1,m
    do j=1,m
      h_m(j,i)=hm(j,i)
      o_m(j,i)=om(j,i)
    end do
  end do
  if(m<20) then
    call diag1(h_m,o_m,m,eigval_m,eigvec_m)
  else
    fail=0
    call diag(h_m,o_m,m,min(m,n),eigval_m,eigvec_m)
    if(fail>0) call diag(h_m,o_m,m,eigval_m,eigvec_m)
  endif
  do i=1,min(n,m)
    eigval(i)=eigval_m(i)
    do j=1,m
      eigvec(j,i)=eigvec_m(j,i)
    end do
  end do

  open(16,file='ener.dat')
  do i=1,n
    write(16,*)i,eigval(i)
  end do
  close(16)

  if(exc) then
    ener=0.d0
    do i=1,nexc
      ener=ener+(eigval(i)-e_low)**2
    end do
    if(excmode==1) ener=sqrt(ener)
    if(excmode==2) ener=eigval(nexc)
  else
    ener=eigval(1)
  endif
  
  deallocate(h_m,o_m,eigval_m,eigvec_m)

end subroutine partial_diag

subroutine lowest_eigval(ener,n)
! to calculate the lowest eigenvalue 
implicit none
  double precision :: ener,ener0,ener1
  integer          :: n  
  double precision :: psi(n),hpsi(n)
  double precision :: xnb,s,x1,x2,xacc,brac,fa,fb
  integer          :: i,j,err
  n_current=n
  psi=0.d0 
  hpsi=0.d0
! overlap of the new basis state with 
! the previously selected orthogonalized basis
  do i=1,n-1
    s=0.d0
    do j=1,n-1
      s=s+eigvec0(j,i)*om(j,n)
    end do
    psi(i)=s
  end do
! norm of the new basis state 
  s=om(n,n)
  do i=1,n-1     
    s=s-psi(i)*psi(i)
  end do
  if(s.le.0.d0) then
    ener=1.99999d+12
    return
  endif
  xnb=dsqrt(s)
! matrix element of the hamiltonian between the previously
! selected and the new basis state
  
  do i=1,n-1
    s=0.d0
    do j=1,n-1
      s=s+eigvec0(j,i)*hm(j,n)
    end do
    hpsi(i)=s
  end do
  do i=1,n-1
    qqv(i)=(hpsi(i)-eigval0(i)*psi(i))/(xnb)
  end do
  s=hm(n,n)
  do i=1,n-1     
    s=s+eigval0(i)*psi(i)*psi(i)-2.d0*hpsi(i)*psi(i)
  end do
  eqv(n)=s/(xnb**2)
  do i=1,n-1
    eqv(i)=eigval0(i)
  end do
  if(n.eq.1) then
    ener=hm(1,1)/om(1,1)
    return
  endif
  xacc=1.0d-16

  call locate_root(eigval0(1),x1,x2,err)
  if(err>0) then
    ener=eigval0(1)
    ener=3.1415926d+12
    return
  endif
  ener0=zbrent(x1,x2,xacc,xacc,err)
  if(err==0) then
    ener=ener0
  else
    ener=12345.d+12
  endif
  if(exc.and.n>nexc) then
    if(excmode==1) then
      do i=1,nexc
        fa=func(eigval(i)+1.d-12)
        fb=func(eigval(i+1)-1.d-12)
        if(fa*fb<0.d0) then
          ener1=zbrent(eigval(i)+1.d-12,eigval(i+1)-1.d-12,xacc,xacc,err)
        else
          ener1=eigval(i+1)
        endif
        ener=ener+(ener1-e_low)**2
      end do
      ener=sqrt(ener)
    endif   
    if(excmode==2) then
      i=nexc-1  
      fa=func(eigval(i)+1.d-12)
      fb=func(eigval(i+1)-1.d-12)
      if(fa*fb<0.d0) then
        ener1=zbrent(eigval(i)+1.d-12,eigval(i+1)-1.d-12,xacc,xacc,err)
      else
        ener1=eigval(i+1)
      endif
      ener=ener1
    endif
    
  endif

end subroutine lowest_eigval




subroutine locate_root(e,a,b,err)
! e       : lowest eigenvalue of n-1
! [x1,x2] : location of the root
implicit none
double precision :: a,b,e
integer          :: err
double precision, parameter :: h=0.1d0,f=1.d-08
double precision :: x
integer          :: i
x=f
a=e
do while (func(a)<0)
  a=a-h
end do
b=e-x
i=0
err=0
do while (func(b)>0)
  x=x/2.d0
  b=e-x
  i=i+1
  if(i>20) then
    err=1
    exit
  endif
end do

!b=0.5d0*(e+a)
!do while (func(b)>0)
!  b=0.5d0*(e+b)
!end do




end subroutine locate_root



function func(x)
implicit none
  double precision :: x,prod,w,ww,func,x1,x2
  integer          :: n,j
  n=n_current
  w=1.d0
  w=w*(eqv(n)-x)
  do j=1,n-1
    ww=qqv(j)*qqv(j)/(eqv(j)-x)
    w=w-ww
  end do
  func=w
end function func

      
FUNCTION zbrent(x1,x2,tol,epsiln,err)
  implicit none
  double precision           :: zbrent,tol,epsiln
  integer, PARAMETER         :: ITMAX=100
  INTEGER                    :: iter,err
  double precision           ::  a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm,x1,x2

  err=0
  a=x1
  b=x2
  fa=func(a)
  fb=func(b)


  if((fa.gt.0.d0.and.fb.gt.0.d0).or.(fa.lt.0.d0.and.fb.lt.0.d0)) then
    write(6,*)'root must be bracketed for zbrent'
    stop
  endif   
  c=b
  fc=fb
  do iter=1,ITMAX
    if((fb.gt.0.d0.and.fc.gt.0.d0).or.(fb.lt.0.d0.and.fc.lt.0.d0))then
      c=a
      fc=fa
      d=b-a
      e=d
    endif
    if(abs(fc).lt.abs(fb)) then
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    endif
    tol1=2.d0*epsiln*abs(b)+0.5d0*tol
    xm=.5d0*(c-b)
    if(abs(xm).le.tol1 .or. fb.eq.0.d0) then
      zbrent=b
      return
    endif
    if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
      s=fb/fa
      if(a.eq.c) then
        p=2.*xm*s
        q=1.d0-s
      else
        q=fa/fc
        r=fb/fc
        p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
        q=(q-1.d0)*(r-1.d0)*(s-1.d0)
      endif
      if(p.gt.0.d0) q=-q
      p=abs(p)
      if(2.d0*p .lt. min(3.d0*xm*q-abs(tol1*q),abs(e*q))) then
        e=d
        d=p/q
      else
        d=xm
        e=d
      endif
    else
      d=xm
      e=d
    endif
    a=b
    fa=fb
    if(abs(d) .gt. tol1) then
      b=b+d
    else
      b=b+sign(tol1,xm)
    endif
    fb=func(b)
  end do
  err=1
  zbrent=b
  return
END function zbrent


subroutine matrix_elements(z1,z2,s1,s2,ome,hme,kme,pme,polm,sp)
  implicit none
  real*8     :: z1(N_pairs+N_particle),z2(N_pairs+N_particle),s1(3,N_particle),s2(3,N_particle)
  real*8     :: ome,hme,pme,kme
!
  real*8     :: a1(N_particle,N_particle),a2(N_particle,N_particle)
  real*8     :: b(N_particle,N_particle),bi(N_particle,N_particle),v(3,N_particle),y(3,N_particle)
  real*8     :: a1b(N_particle,N_particle),a2b(N_particle,N_particle)
  real*8     :: t1(3,N_particle),t2(3,N_particle)
  real*8     :: t(3,N_particle)
  real*8     :: aba(N_particle,N_particle)
  real*8     :: abal(N_particle,N_particle)
  real*8     :: bv(3,N_particle)
  real*8     :: d,vv,tr,yy,ww,s,c,ss,tme,vme,ume,ovme,polm,self_polm,ss1,ss2,sp,ddo
  integer    :: i,k,ip,j,ipos,i1,i2

  
  ome=0.d0 
  hme=0.d0
  kme=0.d0 
  pme=0.d0
  polm=0.d0
  sp=0.d0
  call trcorr(z1,a1,s1,ss1,t)
  t1=t

do ip=1,Nperm
  if(spin_ov(ip).ne.0) then
  call trcorr(z2,a2,s2,ss2,t)
  do i=1,N_particle
    j=ipr(i,ip)
    t2(:,i)=t(:,j)
  end do
  ss=ss1+ss2


  b=a2
  do i=1,N_particle
    do j=1,N_particle
      a2(j,i)=b(ipr(j,ip),ipr(i,ip))
    end do
  end do
  
  b=a1+a2

  call inv(b,N_particle,bi)
  call det(b,N_particle,d)

  v=(t1+t2)*cs
  
  a1b=matmul(a1,bi)
  a2b=matmul(a2,bi)
  aba=matmul(a1b,a2)
  abal=matmul(aba,Lambda)

  y=0.d0
  bv=0.d0
  do k=1,N_particle
    do i=1,N_particle
      y(:,k)=y(:,k)+a2b(k,i)*t1(:,i)-a1b(k,i)*t2(:,i)
      bv(:,k)=bv(:,k)+bi(k,i)*v(:,i)
    end do
  end do


  vv=0.d0
  yy=0.d0
  tr=0.d0
  do i=1,N_particle
    vv=vv+v(1,i)*bv(1,i)+v(2,i)*bv(2,i)+v(3,i)*bv(3,i)
    do j=1,N_particle
      yy=yy+lambda(i,j)*(y(1,i)*y(1,j)+y(2,i)*y(2,j)+y(3,i)*y(3,j))
    end do
    tr=tr+abal(i,i)
  end do

  
! overlap
  ovme=exp(0.5d0*(vv-ss))/d**1.5d0
! kinetic
  tme=h2m*(3.d0*tr-yy)*ovme
! one-body potential
  if(density) then
    do i=1,N_particle         
      c=bi(i,i)
      write(50,'(4d16.8,2i5)')1.d0/c,bv(1,i),bv(2,i),bv(3,i),i,ip
    end do
    do i=1,N_particle
      do k=1,N_particle
        write(52,*)bi(i,k)
      end do
      write(52,*)bv(1,i),bv(2,i),bv(3,i)
    end do
  endif


  ww=0.d0
  do i=1,N_particle
    do k=1,Number_of_atoms
      s=sqrt((bv(1,i)-Position(1,k))**2+(bv(2,i)-Position(2,k))**2+ &
&       (bv(3,i)-Position(3,k))**2)
      c=bi(i,i)
      if(onebody==1) ww=ww+potc(c,s)*Nuclear_Charge(k)*charge(i)
      if(onebody==2) ww=ww+potring(c,s)
      if(onebody==3) ww=ww+pot_ho(c,s)*omega0**2
    end do
  end do

  ume=ww*ovme
! two-body Coulomb potential
  ww=0.d0
  do i=1,N_particle
    do k=i+1,N_particle
      s=sqrt((bv(1,i)-bv(1,k))**2+(bv(2,i)-bv(2,k))**2+(bv(3,i)-bv(3,k))**2)
      c=bi(i,i)+bi(k,k)-2.d0*bi(i,k)
      ww=ww+potc(c,s)*charge(ipr(i,ip))*charge(k)
! adding ho
      if(onebody==4) ww=ww+pot_ho(c,s)*omega0**2/N_particle
    end do
  end do
  vme=ww*ovme
! polariton
  ww=0.d0
  do i=1,N_particle
    ww=ww+lambda_pol(1)*charge(i)*bv(1,i)+lambda_pol(2)*charge(i)*bv(2,i)
  end do
  polm=polm+ww*ovme*spin_ov(ip)
  ww=0.d0
  do i=1,N_particle
    do j=1,N_particle
      ww=ww+charge(ipr(i,ip))*charge(j)*( &
&       lambda_pol(1)**2*(bi(i,j)+bv(1,i)*bv(1,j))+ &
&       lambda_pol(2)**2*(bi(i,j)+bv(2,i)*bv(2,j))+ &
&       2.d0*lambda_pol(1)*lambda_pol(2)*bv(1,i)*bv(2,j))
    end do
  end do
  self_polm=0.5d0*ww*ovme
  sp=sp+self_polm*spin_ov(ip)
  ome=ome+spin_ov(ip)*ovme
  hme=hme+spin_ov(ip)*(tme+ume+vme+ovme*E_ion_ion+self_polm)
  pme=pme+spin_ov(ip)*(ume+vme)
  kme=kme+spin_ov(ip)*tme

  
!  write(11,*)polm
  
!  write(66,*)ip,spin_ov(ip)
!  write(66,*)vv,ss,d
!  write(66,*)'o',ovme
!  write(66,*)'t',tme/ovme,tme
!  write(66,*)tr,yy
! write(66,*)'u',ume/ovme,ume
!  write(66,*)'v',vme/ovme,vme
!  write(66,*)'o',ome
!  write(66,*)'v',hme

!  write(66,*)'p',polm/ovme
!  write(66,*)'s',self_polm/ovme
!  write(66,*)'ome',ome
!  write(66,*)'hme',hme
!  write(66,*)'spi',spin_ov(ip)

  if(density) write(50,*)spin_ov(ip)*ovme,ip
  if(density) write(52,*)spin_ov(ip)*ovme
endif
end do

!  write(9999,*)hme
!  write(9999,*)kme+pme+sp


end subroutine matrix_elements

subroutine matrix_elements_check(z1,z2,s1,s2,ww)
  implicit none
  real*8     :: z1(N_pairs+N_particle),z2(N_pairs+N_particle),s1(3,N_particle),s2(3,N_particle)
  real*8     :: ome,hme
!
  real*8     :: a1(N_particle,N_particle),a2(N_particle,N_particle)
  real*8     :: b(N_particle,N_particle),bi(N_particle,N_particle),v(3,N_particle)
  real*8     :: t1(3,N_particle),t2(3,N_particle)
  real*8     :: bv(3,N_particle)
  real*8     :: d,vv,tr,yy,ww,s,c,ss,ss1,ss2
  integer    :: i,k,ip,j

  call trcorr(z1,a1,s1,ss1,t1)
  call trcorr(z2,a2,s2,ss2,t2)
  
  b=a1+a2
  call inv(b,N_particle,bi)
  call det(b,N_particle,d)
  v=t1+t2
  
  do k=1,N_particle
    do i=1,N_particle
      bv(:,k)=bv(:,k)+bi(k,i)*v(:,i)
    end do
  end do

  ss=ss1+ss2
  vv=0.d0
  yy=0.d0
  tr=0.d0
  do i=1,N_particle
    vv=vv+v(1,i)*bv(1,i)+v(2,i)*bv(2,i)+v(3,i)*bv(3,i)
  end do
  ww=vv-ss

end subroutine matrix_elements_check


function potc(c,s)
  implicit none
  real*8 :: c,s,potc

  if(s.eq.0.d0) then
    potc=sqrt(0.5d0/c)*(2.d0/sqrt(pi))
  else
    potc=erf(sqrt(0.5d0/c)*s)/s
  endif

end function potc




function potring(cc,s)
  implicit none
  real*8 :: cc,s,z,c,ho,a,ga,potring

  c=1.d0/cc

  a=1.d0/0.997d0**2
  ho=3.d0/c+s**2
  ga=c**1.5d0/(c+2.d0*a)**1.5d0*exp(-1.d0*a*c/(c+2.d0*a)*s**2)
  
  potring=0.5d0*0.7827d0*ho+17.70d0*ga
  
end function potring


function pot_ho(cc,s)
  implicit none
  real*8 :: cc,s,z,c,ho,a,ga,pot_ho

  c=1.d0/cc

  ho=3.d0/c+s**2
  pot_ho=0.5d0*ho
  
end function pot_ho

subroutine trcorr(z,a,s,ss,t)
! transform the nonlinear parameters
  implicit none
  real*8   :: a(N_particle,N_particle),z(N_pairs+N_particle)
  real*8   :: b(N_particle,N_particle),s(3,N_particle),t(3,N_particle),ss
  integer  :: i,j,k

  
  k=0
  a=0.d0
  do i=1,N_particle
    do j=i+1,N_particle
      k=k+1
      a(i,j)=-z(k)
      a(j,i)=-z(k)
    end do
  end do
  do i=1,N_particle
    do j=1,N_particle
      if(i.ne.j) a(i,i)=a(i,i)-a(i,j)
    end do
  end do

  if(cm) then
    b=0.d0
    do i=1,N_particle
      do j=1,N_particle
        do k=1,N_particle
          b(k,j)=b(k,j)+umat(i,k)*umat(i,j)*z(n_pairs+i)
        end do
      end do
    end do

    ss=0.d0
    t=0.d0
    do i=1,n_particle
      do j=1,n_particle
        ss=ss+(s(1,i)*s(1,j)+s(2,i)*s(2,j)+s(3,i)*s(3,j))*b(i,j)
        t(:,j)=t(:,j)+b(j,i)*s(:,i)
      end do
    end do
  else
    b=0.d0
    do i=1,n_particle
      b(i,i)=z(n_pairs+i)
    end do
    ss=0.d0
    t=0.d0
    do i=1,n_particle
      ss=ss+(s(1,i)*s(1,i)+s(2,i)*s(2,i)+s(3,i)*s(3,i))*b(i,i)
      t(:,i)=t(:,i)+b(i,i)*s(:,i)
    end do
  endif
  a=a+b

  if(cm)  a=a+0.5d0*cm_mat



  
end subroutine trcorr




subroutine det(am,n,d)
  implicit none
  integer   :: n
  real*8    :: am(n,n),d
!
  integer   :: i,j,k,ii,n1
  real*8    :: a(n,n)

  a=am
  if(n.eq.1) then
    d=a(1,1)
    return
  endif

  n1=n-1
  do i=1,n1
    ii=i+1
    do  k=ii,n
      d=a(i,k)/a(i,i)
      do j=ii,n
        a(j,k)=a(j,k)-a(j,i)*d
      end do
    end do
  end do
  d=1.d0
  do i=1,n
    d=d*a(i,i)
  end do
end subroutine det

subroutine inv(e,n,a)
!     inverse (real symmetric matrix)
  implicit none
  integer  :: n
  real*8   :: a(n,n),e(n,n)
!
  real*8   :: b(n,2*n),x
  integer  :: i,j,m,j2,k

  if(n.eq.1) then
    a(1,1)=1.d0/e(1,1)
    return
  endif

  b=0.d0
  j2=n*2
  do i=1,n
    b(i,i+n)=1.d0
  end do
  do i=1,n
    do j=1,n
      b(j,i)=e(j,i)
    end do
  end do
  do i=1,n
    x=b(i,i)
    do k=i,j2
      b(i,k)=b(i,k)/x
    end do
    do m=1,n
      x=b(m,i)
      if (m.ne.i) then
        do k=i,j2
          b(m,k)=b(m,k)-b(i,k)*x
        end do
      endif
    end do
  end do
  do j=1,n
    do i=1,n
      a(i,j)=b(i,j+n)
    end do
  end do

end  subroutine inv

      SUBROUTINE PERMUT(NRP,N,IA)
!     permutations 
      implicit none
      integer   :: nrp,n,i,in,io,m
      integer   :: ia(N_particle_max),ifct(0:N_particle_max),iv(N_particle_max+1)
      DATA (ifct(i),i=0,N_particle_max) /1,1,2,6,24,120,720,5040,40320, &
&     362880,3628800,39916800,479001600/

       do i=1,n
         iv(i)=i
       end do

       IO=NRP-1
       DO M=N-1,1,-1
         IN=IO/IFCT(M)+1
         IO=MOD(IO,IFCT(M))
         IA(N-M)=IV(IN)
         DO I=IN,M
           IV(I)=IV(I+1)
         end do
       end do
       IA(N)=IV(1)
      return
      end subroutine permut
                       
subroutine pape(ia,n,ip)
!     parity of a given permutation
      implicit none
      integer :: i,n,ip,ii,mm,j,ibf
      integer :: ia(N_particle_max),ib(N_particle_max)
      ibf=1
      do  i=1,n
        ib(i)=ia(i)
      end do
      ii=0
      do i=1,n
        do j=i,n
          if(i.eq.ib(j).and.i.ne.j) then
            mm=ib(i)
            ib(i)=ib(j)
            ib(j)=mm
            ii=ii+1
          endif    
        end do
      end do      
      ip=1
      if(ibf.eq.1) ip=(-1)**ii  
      return
end subroutine pape

subroutine input_basis
implicit none
  integer    :: i,j,k
  logical                :: lex
  
  inquire(file='cg.out',exist=lex)
    if(lex) then
      open(1,file='cg.out')
      read(1,*)N_basis_start
      do j=1,N_basis_start
        write(6,*)j
        read(1,*)i,photons(j),parity(j)
        read(1,*)(Alpha(i,j),i=1,N_pairs+N_particle)
        do i=1,N_particle
          read(1,*)(S_vector(k,i,j),k=1,3)
        end do
      end do
      close(1)
    endif
end  subroutine input_basis

subroutine output_basis(n)
implicit none
  integer    :: i,j,n,k
  open(1,file='cg.out')
  write(1,*)N
  do j=1,N
    write(1,*)j,photons(j),parity(j)
    write(1,*)(Alpha(i,j),i=1,N_pairs+N_particle)
    do i=1,N_particle
      write(1,*)(S_vector(k,i,j),k=1,3)
    end do
  end do
  close(1)  
end  subroutine output_basis


subroutine symm(alpha1,alpha2,s1,s2,lambda_norm,ome,hme,tme,pme,polm,self_pol,i1,i2)
! symmetryzing around lambda (lambda_norm is normalized)
implicit none
  double precision :: s1(3,n_particle),s2(3,n_particle),lambda_norm(2)
  double precision :: alpha1(n_particle),alpha2(n_particle)
  double precision :: ome,hme,polm,self_pol,pme,tme
  double precision :: s(3),sp(3),s1p(3,n_particle),s2p(3,n_particle)
  double precision :: o,h,p,sep,t,v
  integer          :: i,j,k,nn,i1,i2,ip
  integer          :: b(32)

      hme=0.d0
      ome=0.d0
      pme=0.d0
      tme=0.d0
      polm=0.d0
      self_pol=0.d0

      s1p=s1
      s2p=s2
      call matrix_elements(alpha1,alpha2,s1p,s2p,o,h,t,v,p,sep)    
      ome=ome+o
      hme=hme+h
      pme=pme+v
      tme=tme+t
      polm=polm+p
      self_pol=self_pol+sep
    
!      write(999,*)hme
!      write(999,*)pme+tme+self_pol
      
end subroutine symm


subroutine integer2binary(i,b)
implicit none
  integer,intent(in) :: i
  integer :: b(32)
  integer k,j
  b=0
  j=i
  do k=1,size(b)
    b(k)=mod(j,2)
    j=j/2
  enddo
end subroutine integer2binary

subroutine reflect(x,d,xp)
!reflecting a vector around a normalized vector 
implicit none
  double precision :: x(2),d(2),xp(2),xd

  xd=x(1)*d(1)+x(2)*d(2)
  xp(1)=-x(1)+2.d0*xd*d(1)
  xp(2)=-x(2)+2.d0*xd*d(2)

end subroutine reflect

subroutine calculate_row(n)
implicit none 
  real*8     :: ome,hme,polm,self_pol,tme,pme
  real*8     :: s1(3,n_particle),s2(3,n_particle),d(2)
  integer    :: i,j,n,n1,n2,i1,i2

  d=lambda_norm

  
  
  i1=parity(n)
  
  n1=photons(n)
  do i=1,n

!    write(66,*)'--------------------------------',i,n
    write(50,*)i,n
    n2=photons(i)
    i2=parity(i)

    hr(i)=0.d0
    or(i)=0.d0
    pr0(i)=0.d0
    tr0(i)=0.d0
    sr0(i)=0.d0
    vr0(i)=0.d0
    wr0(i)=0.d0
!    if(full_matrix) then
!      s1(1,:)=s_vector(1,:,i)
!      s1(2,:)=s_vector(2,:,i)
!      s1(3,:)=s_vector(3,:,i)
!      s2(1,:)=s_vector(1,:,n)
!      s2(2,:)=s_vector(2,:,n)
!      s2(3,:)=s_vector(3,:,n)
!      call symm(Alpha(:,i),Alpha(:,n),S1,S2,d,ome,hme,tme,pme,polm,self_pol,i1,i2)
!      hr0(i)=hme
!      or0(i)=ome
!      sr0(i)=self_pol
!      tr0(i)=tme
!      vr0(i)=pme
!    endif
    if(abs(n1-n2)<2) then
      s1(1,:)=s_vector(1,:,i)
      s1(2,:)=s_vector(2,:,i)
      s1(3,:)=s_vector(3,:,i)
      s2(1,:)=s_vector(1,:,n)
      s2(2,:)=s_vector(2,:,n)
      s2(3,:)=s_vector(3,:,n)
      if(density.ne..true.) call symm(Alpha(:,i),Alpha(:,n),S1,S2,d,ome,hme,tme,pme,polm,self_pol,i1,i2)
      if(n1==n2) then
        hr(i)=hme+omega*n1*ome
        wr0(i)=omega*n1*ome
        hr0(i)=hme
        or(i)=ome
        sr0(i)=self_pol
        tr0(i)=tme
        vr0(i)=pme
      endif
      if(n1==n2+1) then
        hr(i)=sqrt(omega/2.d0)*polm*sqrt(1.d0*(n1))
        pr0(i)=sqrt(omega/2.d0)*polm*sqrt(1.d0*(n1))
        or(i)=0.d0
        hr0(i)=0.d0
        wr0(i)=0.d0
        tr0(i)=0.d0
        sr0(i)=0.d0
        vr0(i)=0.d0
      endif
      if(n1==n2-1) then
        hr(i)=sqrt(omega/2.d0)*polm*sqrt(1.d0*(n1+1))
        pr0(i)=sqrt(omega/2.d0)*polm*sqrt(1.d0*(n1+1))
        or(i)=0.d0
        hr0(i)=0.d0
        wr0(i)=0.d0
        tr0(i)=0.d0
        sr0(i)=0.d0
        vr0(i)=0.d0
      endif
    end if
  end do
  
end  subroutine calculate_row

subroutine calculate_matrix_elements(n)
! calculate matrix elements up to n
  implicit none
  integer                 :: i,j,n

  do i=1,n
    call calculate_row(i)
    do j=1,i
      hm(j,i)=hr(j)
      hm(i,j)=hr(j)
      wm(j,i)=wr0(j)
      wm(i,j)=wr0(j)
      hm0(j,i)=hr0(j)
      hm0(i,j)=hr0(j)
      tm0(j,i)=tr0(j)
      tm0(i,j)=tr0(j)
      vm0(j,i)=vr0(j)
      vm0(i,j)=vr0(j)
      sm(j,i)=sr0(j)
      sm(i,j)=sr0(j)
      om(j,i)=or(j)
      om(i,j)=or(j)      
      om0(j,i)=or0(j)
      om0(i,j)=or0(j)      
      pm(j,i)=pr0(j)
      pm(i,j)=pr0(j)
      write(99,*)j,i
      write(99,*)hm(j,i)
      write(99,*)tm0(j,i)+vm0(j,i)+sm(j,i)+pm(j,i)+wm(j,i),pm(j,i)
    end do
  end do

end subroutine calculate_matrix_elements

subroutine svm_precal
! calculate the energy  on a given basis
  implicit none

  double precision               :: theta,ener
  integer                        :: i,j,k,n
  double precision               :: energies(n_basis_start,n_basis_start)
    
  n=100
  open(1,file='cs.dat')
  write(1,*)n,n_basis_start
  do k=0,n
    write(6,*)k
    theta=0.002*k
    cs=cos(theta)
    write(1,*)theta
    write(6,*)theta
    call calculate_matrix_elements(N_basis_start)
    do i=1,n_basis_start  
      do j=1,n_basis_start  
        write(1,*)hm(j,i),hm0(j,i),tm0(j,i),vm0(j,i),sm(j,i),om(j,i),pm(j,i),wm(j,i)
      end do
    end do
  end do
  close(1)

 call partial_diag(n_basis_start,n_basis_start,ener)
 write(6,*)ener
  
end subroutine svm_precal




   
END MODULE CORRELATED_GAUSSIANS


USE CORRELATED_GAUSSIANS
implicit none
integer :: ip,k

  call initialize
  call input_basis
  full_matrix=.true.

  call svm_precal




end
