module geometry
  use kinds
  implicit none

  integer, parameter :: nx = 48, ny = 24
  integer, parameter :: n = nx * ny
  real(dp), parameter :: lx = 20.0_dp, ly = 10.0_dp
  real(dp), parameter :: dx = lx / (nx + 1), dy = ly / (ny + 1)
  real(dp), parameter :: hbar = 1.0_dp, mass = 1.0_dp

  type :: disk_t
    real(dp) :: x, y, r
  end type disk_t

  type(disk_t), allocatable :: obstacles(:)

contains

  subroutine set_obstacles(preset)
    integer, intent(in) :: preset
    select case (preset)
    case (0)
      allocate(obstacles(0))
    case (1)
      allocate(obstacles(1))
      obstacles(1) = disk_t(10.0_dp, 5.0_dp, 1.6_dp)
    case default
      allocate(obstacles(5))
      obstacles(1) = disk_t(6.0_dp, 3.0_dp, 0.9_dp)
      obstacles(2) = disk_t(14.0_dp, 7.0_dp, 0.9_dp)
      obstacles(3) = disk_t(10.0_dp, 5.0_dp, 1.1_dp)
      obstacles(4) = disk_t(5.0_dp, 7.5_dp, 0.8_dp)
      obstacles(5) = disk_t(15.0_dp, 2.5_dp, 0.8_dp)
    end select
  end subroutine set_obstacles

  pure integer function grid_index(i, j)
    integer, intent(in) :: i, j
    grid_index = (j - 1) * nx + i
  end function grid_index

  pure real(dp) function grid_x(i)
    integer, intent(in) :: i
    grid_x = i * dx
  end function grid_x

  pure real(dp) function grid_y(j)
    integer, intent(in) :: j
    grid_y = j * dy
  end function grid_y

  real(dp) function potential(x, y)
    real(dp), intent(in) :: x, y
    real(dp), parameter :: amp = 400.0_dp, sig = 0.45_dp
    integer :: k
    potential = 0.0_dp
    do k = 1, size(obstacles)
      potential = potential + amp * exp(-((x - obstacles(k)%x)**2 + &
                  (y - obstacles(k)%y)**2) / (2.0_dp * sig * sig))
    end do
  end function potential

  logical function inside_obstacle(x, y)
    real(dp), intent(in) :: x, y
    integer :: k
    inside_obstacle = .false.
    do k = 1, size(obstacles)
      if ((x - obstacles(k)%x)**2 + (y - obstacles(k)%y)**2 <= obstacles(k)%r**2) then
        inside_obstacle = .true.
        return
      end if
    end do
  end function inside_obstacle

end module geometry
