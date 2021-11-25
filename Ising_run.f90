program Ising_run
  use Ising_lib

  CHARACTER(100) :: num1char
  REAL(8) :: T, Mavg, mavg_
  type(IsingSystem) :: sys
  integer :: n = 50, j, k
  character*3 :: fname1, fname2
  logical :: writeAnimation = .false.
  
  CALL GET_COMMAND_ARGUMENT(1,num1char)
  READ(num1char,*) T
  
  sys = new_IsingSystem(n)
  call sys%randomize()
  call sys%set_Temperature(T)
  
  Mavg = 0.d0
  mavg_ = 0.d0
  
  if (n .le. 30) call sys%writeSystem()
  write(fname1,'(F3.1)') T
  fname1 = fname1(1:1)//"-"//fname1(3:3)
  open(unit=12, file='out'//fname1//'.csv')
  
  do j=1,10000
    write(12,*) sys%M()/n**2d0, sys%M()
    Mavg = Mavg + sys%M()
    mavg_ = mavg_ + sys%M()/n**2d0
    if (writeAnimation .and. mod(j,10) == 0) then
      call writeOutForAnimation(j/10, sys)
    endif
    
    call sys%doStep()
  enddo
  
  Mavg = Mavg/10000
  mavg_ = mavg_/10000
  
  close(12)
  if (n .le. 30) call sys%writeSystem()
  
  open(unit=14, file='averages.csv', POSITION='APPEND')
  write(14,*) T, Mavg, mavg_
  write(*,*) mavg_
  close(14)

contains

subroutine writeOutForAnimation(Nr, sys)
  integer :: Nr
  class(IsingSystem) :: sys
  if (Nr .lt. 10) then
    write(fname2,'(i1)') Nr
    fname2 = trim(fname2)
  else if (Nr .lt. 100) then
    write(fname2,'(i2)') Nr
  else if (Nr .lt. 1000) then
    write(fname2,'(i3)') Nr
  end if
  open(unit=13, file='./animation/system'//fname1//"_"//trim(fname2)//'.dat')
  do k=1,sys%N
    write(13,*) sys%Sys(k,:)
  enddo
  close(13)
end subroutine

end program
