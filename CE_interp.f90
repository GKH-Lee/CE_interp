module CE_interp


! T-p and x data arrays from table
  integer :: np, nmol, u
  character(len=10), allocatable, dimension(:) :: m_name
  real*8, allocatable, dimension(:) :: P_t, T_t, lP_t, lT_t
  real*8, allocatable, dimension(:,:,:) :: x_t
  real*8, allocatable, dimension(:) :: x_out

! hard code binary file name to read
  character(len=100) :: x_table = 'interp_table.bin'



public :: CE_interpolate
private :: CE_interpolate_init, locate

contains

  subroutine CE_interpolate(P_in, T_in, m_in, m_size, x_out)
    implicit none

    integer, intent(in) :: m_size
    real*8, intent(in) :: P_in, T_in
    character(len=10), dimension(m_size), intent(in) :: m_in
    real*8, intent(out), dimension(m_size) :: x_out


    integer :: m
    integer :: i_pl, i_pu, i_tl, i_tu
    real*8 :: x0,x1,x2,y0,y1,y2,a11,a21,a12,a22,norm

    logical :: first_call = .True.


    ! Initialise - read VMR from table for T and p grid
    if (first_call .eqv. .True.) then
      call CE_interpolate_init()
      first_call = .False.
    end if

   ! Find upper and lower T and P indexes
   call locate(P_t, np, P_in, i_pl)
   i_pu = i_pl + 1
   call locate(T_t, np, T_in, i_tl)
   i_tu = i_tl + 1

   x0 = log10(T_in)
   x1 = lT_t(i_tl)
   x2 = lT_t(i_tu)
   y0 = log10(P_in)
   y1 = lP_t(i_pl)
   y2 = lP_t(i_pu)

   norm = 1.0d0 / (x2 - x1) / (y2 - y1)

   ! Perform bi-linear interpolation for each molecule
   do m = 1, m_size
     ! Find VMR from table
     a11 = x_t(m,i_pl,i_tl)
     a21 = x_t(m,i_pu,i_tl)
     a12 = x_t(m,i_pl,i_tu)
     a22 = x_t(m,i_pu,i_tu)
 
     ! Bi-linear equation
     x_out(m) = a11 * (x2 - x0) * (y2 - y0) * norm &
              & + a21 * (x0 - x1) * (y2 - y0) * norm &
              & + a12 * (x2 - x0) * (y0 - y1) * norm &
              & + a22 * (x0 - x1) * (y0 - y1) * norm
   end do

   ! Un-log10 the VMR on output
   x_out(:) = 10.0d0**x_out(:)

  end subroutine CE_interpolate


  subroutine CE_interpolate_init()
    implicit none


    integer :: i, j

    open(newunit=u, file=x_table,action='read',form='unformatted')

    read(u) np, nmol
 
    print*, np, nmol
 
    allocate(m_name(nmol))
    read(u) m_name(1:nmol)

    print*, m_name

    allocate(T_t(np))
    do i = 1, np
      read(u) T_t(i)
    end do

    print*, T_t

    allocate(P_t(np))
    do j = 1, np
      read(u) P_t(j)
    end do

    print*, P_t

    allocate(x_t(nmol,np,np))
    do j = 1, np
      do i = 1, np
        read(u) x_t(1:nmol,j,i)
      end do
    end do

    ! Log10 arrays of grid
    allocate(lT_t(np),lP_t(np))
    lT_t = log10(T_t)
    lP_t = log10(P_t)

  end subroutine CE_interpolate_init


  pure subroutine locate(arr, n, var, idx)
    implicit none

    integer, intent(in) :: n
    integer, intent(out) :: idx
    real*8, dimension(n), intent(in) :: arr
    real*8,intent(in) ::  var
    integer :: jl, jm, ju

    ! Search an array using bi-section (numerical methods)
    ! Then return array index that is lower than var in arr

    jl = 0
    ju = n+1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if ((arr(n) > arr(1)).eqv.(var > arr(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate

end module CE_interp



program CE
  use CE_interp
  implicit none

  real*8 :: T = 3000, P = 0.001
  integer :: lnames = 8
  character(len=10), dimension(8) :: names
  real*8, dimension(8) :: x

  names = (/'H2O','CO ','CH4','NH3','H2 ' ,'He ','Na ','K  '/)
 
  call CE_interpolate(P,T,names,lnames,x)

  print*, T,P,x

end program CE

