use linalg
implicit none
  integer                       :: i,j,n,k
  double precision, allocatable :: h(:,:),o(:,:),e(:),ev(:,:)
  double precision, allocatable :: hm(:,:),om(:,:),em(:),evm(:,:),ener(:,:)
  
  read(55,*)n
  allocate(h(n,n),o(n,n),ev(n,n),e(n),ener(n,n))
  do i=1,n
    do j=1,n
      read(55,*)h(j,i),o(j,i)
    end do
  end do

  do k=1,n
    allocate(hm(k,k),om(k,k),evm(k,k),em(k))
    do i=1,k
      do j=1,k
        hm(j,i)=h(j,i)
        om(j,i)=o(j,i)
      end do
    end do

    if(k>50) then
      call diag(hm,om,k,em,evm)
      write(6,*)k,em(1)
      do i=1,50
        ener(k,i)=em(i)
      end do
    endif
    deallocate(hm,om,em,evm)
  end do
  open(1,file='spectrum.dat')
  do k=51,n
    write(1,'(i5,50d16.8)')k,(ener(k,i),i=1,50)
  end do
  close(1)



  end
