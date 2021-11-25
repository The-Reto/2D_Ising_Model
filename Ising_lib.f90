module Ising_lib
  type IsingSystem
    integer :: N
    real(8) :: J, T, expm4J, expm2J
    integer(kind=1), dimension(:,:), allocatable :: sys

  contains
    procedure :: init=>init
    procedure :: randomize=>randomize
    procedure :: S => S
    procedure :: writeSystem=>writeSystem
    procedure :: E => E
    procedure :: M => magnet
    procedure :: dE => dE
    procedure :: flip => flip
    procedure :: doStep=>performStep
    procedure :: set_Temperature => setT
  end type

  interface IsingSystem
    procedure :: new_IsingSystem
  end interface
contains

  subroutine setT(this, T)
    class(IsingSystem) :: this
    real(8) :: T
    this%T = T
    this%expm2j = exp(-4d0*this%J/this%T)
    this%expm4j = exp(-8d0*this%J/this%T)
  end subroutine

  type(IsingSystem) function new_IsingSystem(N)
    integer :: N
    new_IsingSystem%N = N
    new_IsingSystem%J = 1d0
    new_IsingSystem%T = 1d0
    call new_IsingSystem%init()
    new_IsingSystem%Sys = 1
  end function

  subroutine init(this)
    class(IsingSystem) :: this
    if (allocated(this%sys)) deallocate(this%sys)
    allocate(this%sys(this%N, this%N))
  end subroutine

  subroutine randomize(this)
    class(IsingSystem) :: this
    integer :: i,j
    real, dimension(this%N, this%N) :: r
    call random_number(r)
    this%sys = sign(1d0, r - 5d-1)
  end subroutine

  integer*1 pure function S(this, i,j)
    class(IsingSystem), intent(in) :: this
    integer, intent(in) :: i,j
    S = this%sys(modulo(i-1, this%N) + 1, modulo(j-1, this%N) + 1)
  end function

  subroutine flip(this, i,j)
    class(IsingSystem) :: this
    integer, intent(in) :: i,j
    this%sys(modulo(i-1, this%N) + 1, modulo(j - 1, this%N) + 1) = -this%sys(modulo(i-1, this%N) + 1, modulo(j-1, this%N) + 1)
  end subroutine

  subroutine writeSystem(this)
    class(IsingSystem) :: this
    integer :: i
    write(*,fmt="(A)",advance='no') trim(" ╔")
    do i = 1, 5*this%N + 3
      write(*,fmt="(A)",advance='no') "═"
    enddo
    write(*,fmt="(A)") trim("╗")
    do i = 1,this%N
      write(*,*) trim("║"), this%sys(i,:), trim("  ║")
    enddo
    write(*,fmt="(A)",advance='no') trim(" ╚")
    do i = 1, 5*this%N + 3
      write(*,fmt="(A)",advance='no') "═"
    enddo
    write(*,fmt="(A)") trim("╝")
  end subroutine

  real function E(this)
    class(IsingSystem) :: this
    integer :: i,j
    E = 0.d0
    do i = 1,this%N
      do j = 1,this%N
        E = E + this%S(i,j) * (this%S(i+1, j) + this%S(i, j+1))
      enddo
    enddo
    E = -this%J*E
  end function

  real pure function dE(this, i, j)
    class(IsingSystem), intent(in) :: this
    integer, value :: i,j
    dE = -2.d0*this%J*this%S(i,j)*(this%S(i+1, j) + this%S(i-1, j) + this%S(i, j + 1) + this%S(i, j - 1))
  end function
  
  real function magnet(this)
    class(IsingSystem), intent(in) :: this
    integer :: i,j
    magnet = 0d0
    do i = 1, this%N
      do j = 1, this%N
        magnet = magnet + this%S(i,j)
      enddo
    enddo
  end function
  
  real(8) function p(sys, i,j)
    class(IsingSystem) :: sys
    integer :: i, j, s
    s = sys%S(i, j)*(sys%S(i+1, j) + sys%S(i-1, j) + sys%S(i, j + 1) + sys%S(i, j - 1))
    if (s .le. 0) then
      p = 1.d0
    else if (s == 2) then
      p = sys%expm2j
    else if (s == 4) then
      p = sys%expm4j
    else
      write(*,*) "This shouldn't happen", s
      stop
    endif
  end function

  subroutine performStep(this)
    class(IsingSystem) :: this
    integer :: i,j,k
    real :: i_, j_, r
    do k = 1,1000
      call random_number(r)
      call random_number(i_)
      call random_number(j_)
      i = int(i_ * this%N) + 1
      j = int(j_ * this%N) + 1
      if (r .lt. p(this, i, j)) call this%flip(i,j)
    enddo
  end subroutine
end module
