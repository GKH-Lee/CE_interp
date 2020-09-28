!!!
! G.K.H. Lee - Sep 2020
! A self contained module to interpolate mmr and VMR data from 2D GGchem data,
! as prepared by CE_prepare.py
! (NOTE: you have to modify GGChem to output mean molecular weight too!!!)
!
!
!! TODO: Add boundary checks for pressure range
!!!

module CE_interp_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  integer, parameter :: qp = REAL128


  ! T-p and x data arrays from table
  integer :: np, nmol, u
  character(len=4), allocatable, dimension(:) :: m_name
  real(kind=dp), allocatable, dimension(:) :: P_t, T_t, lP_t, lT_t
  real(kind=dp), allocatable, dimension(:,:) :: mu_t
  real(kind=dp), allocatable, dimension(:,:,:) :: x_t

  ! hard code file produced from CE_prepare.py
  character(len=20) :: x_table = 'interp_table.txt'

  ! First call to module flag
  logical :: first_call = .True.

  public :: CE_interpolate
  private :: CE_interpolate_init, locate, linear_interp, bilinear_interp

contains

  subroutine CE_interpolate(P_in, T_in, m_in, m_size, x_out, mu_out)
    implicit none

    !! Input:
    ! - m_size = size of molecule array
    ! - P_in = input pressure [bar]
    ! - T_in = input temperature [K]
    ! - m_in = name array of species to interpolate
    !! Output:
    ! - x_out = VMR of each species @ P_in and T_in
    ! - mu_out = mean molecular weight @ P_in and T_in
    integer, intent(in) :: m_size
    real(kind=dp), intent(in) :: P_in, T_in
    character(len=4), dimension(m_size), intent(in) :: m_in
    real(kind=dp), intent(out), dimension(m_size) :: x_out
    real(kind=dp), intent(out) :: mu_out

    !! Work variables
    integer :: m
    integer :: i_pl, i_pu, i_tl, i_tu
    real(kind=dp) :: xval,x1,x2,yval,y1,y2,a11,a21,a12,a22,aval

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

   ! print*, P_t(i_pl), P_t(i_pu)
   ! print*, T_t(i_tl), T_t(i_tu)
   ! print*, mu_t(i_pl,i_tu), mu_t(i_pu,i_tu)

   if (i_tl == np) then
     !! Input higher than T grid boundary, make it = minval(T)
     !! Perform mu linear interpolation
     xval = log10(P_in)
     x1 = lP_t(i_pl)
     x2 = lP_t(i_pu)
     y1 = mu_t(i_pl,np)
     y2 = mu_t(i_pu,np)
     call linear_interp(xval, x1, x2, y1, y2, yval)
     mu_out = yval

     do m = 1, m_size
       y1 = x_t(m,i_pl,np)
       y2 = x_t(m,i_pu,np)
       call linear_interp(xval, x1, x2, y1, y2, yval)
       x_out(m) = 10.0_dp**yval ! Unlog for output
     end do

   else if (i_tl == 0) then
     !! Input lower than T grid boundary, make it = minval(T)
     !! Perform mu linear interpolation
     xval = log10(P_in)
     x1 = lP_t(i_pl)
     x2 = lP_t(i_pu)
     y1 = mu_t(i_pl,1)
     y2 = mu_t(i_pu,1)
     call linear_interp(xval, x1, x2, y1, y2, yval)
     mu_out = yval

     do m = 1, m_size
       y1 = x_t(m,i_pl,1)
       y2 = x_t(m,i_pu,1)
       call linear_interp(xval, x1, x2, y1, y2, yval)
       x_out(m) = 10.0_dp**yval ! Unlog for output
     end do

   else
     !! Within T grid bounds
     !! Coordinates for T-p bi-linear interpolation
     xval = log10(T_in)
     x1 = lT_t(i_tl)
     x2 = lT_t(i_tu)
     yval = log10(P_in)
     y1 = lP_t(i_pl)
     y2 = lP_t(i_pu)

     !! Perform bi-linear interpolation for mu
     a11 = mu_t(i_pl,i_tl)
     a21 = mu_t(i_pu,i_tl)
     a12 = mu_t(i_pl,i_tu)
     a22 = mu_t(i_pu,i_tu)
     call bilinear_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
     mu_out = aval

     !! Perform bi-linear interpolation for each species in a loop
     do m = 1, m_size
       a11 = x_t(m,i_pl,i_tl)
       a21 = x_t(m,i_pu,i_tl)
       a12 = x_t(m,i_pl,i_tu)
       a22 = x_t(m,i_pu,i_tu)
       call bilinear_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
       x_out(m) = 10.0_dp**aval ! Unlog for output
     end do
   end if

  end subroutine CE_interpolate


  subroutine CE_interpolate_init()
    implicit none

    integer :: i, j

    !! Read T and P grid from file + mu and VMR data

    open(newunit=u, file=x_table,action='read',form='formatted')

    read(u,*) np, nmol

    !print*, np, nmol

    allocate(m_name(nmol))
    do i = 1, nmol
      read(u,*) m_name(i)
      !print*, i, m_name(i)
    end do

    allocate(T_t(np))
    do i = 1, np
      read(u,*) T_t(i)
      !print*, i,np,T_t(i)
    end do

    allocate(P_t(np))
    do j = 1, np
      read(u,*) P_t(j)
      !print*, j, np, P_t(j)
    end do

    allocate(mu_t(np,np), x_t(nmol,np,np))
    do j = 1, np
      do i = 1, np
        read(u,*) mu_t(j,i), x_t(1:nmol,j,i)
      end do
    end do

    ! Log10 arrays of T-p grid
    allocate(lT_t(np),lP_t(np))
    lT_t = log10(T_t)
    lP_t = log10(P_t)

  end subroutine CE_interpolate_init


  pure subroutine locate(arr, n, var, idx)
    implicit none

    integer, intent(in) :: n
    integer, intent(out) :: idx
    real(kind=dp), dimension(n), intent(in) :: arr
    real(kind=dp), intent(in) ::  var
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

  pure subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(kind=dp), intent(in) :: xval, y1, y2, x1, x2
    real(kind=dp), intent(out) :: yval
    real(kind=dp) :: norm

    norm = 1.0_dp / (x2 - x1)

    yval = (y1 * (x2 - xval) + y2 * (xval - x1)) * norm

  end subroutine linear_interp


  pure subroutine bilinear_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
    implicit none

    real(kind=dp), intent(in) :: xval, yval, x1, x2, y1, y2, a11, a21, a12, a22
    real(kind=dp), intent(out) :: aval
    real(kind=dp) :: norm

    norm = 1.0_dp / (x2 - x1) / (y2 - y1)

    aval = a11 * (x2 - xval) * (y2 - yval) * norm &
      & + a21 * (xval - x1) * (y2 - yval) * norm &
      & + a12 * (x2 - xval) * (yval - y1) * norm &
      & + a22 * (xval - x1) * (yval - y1) * norm

  end subroutine bilinear_interp

end module CE_interp_mod



program CE_interp_test
  use CE_interp_mod, only : CE_interpolate
  implicit none

  real*8 :: T = 1000.0, P = 1.0
  integer, parameter :: lnames = 3
  character(len=4), dimension(lnames) :: names
  real*8, dimension(lnames) :: x
  real*8 :: mu

  names = (/'H2  ','H   ', 'He  '/)

  call CE_interpolate(P,T,names,lnames,x,mu)

  print*, T, P, mu
  print*, x

end program CE_interp_test
