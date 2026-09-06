module evolve
  use kinds
  use geometry
  implicit none

  type :: wavepacket_t !lol im a comment
    real(dp) :: x0, y0, kx0, ky0, sigma
  end type wavepacket_t

contains

  subroutine init_wavepacket(wp, psi0)
    type(wavepacket_t), intent(in) :: wp
    complex(dp), intent(out) :: psi0(n)
    real(dp) :: x, y, norm
    integer :: i, j, idx

    norm = 0.0_dp
    do j = 1, ny
      do i = 1, nx
        idx = grid_index(i, j)
        x = grid_x(i)
        y = grid_y(j)
        psi0(idx) = exp(-((x - wp%x0)**2 + (y - wp%y0)**2) / (4.0_dp * wp%sigma**2)) * &
                    cmplx(cos(wp%kx0 * x + wp%ky0 * y), sin(wp%kx0 * x + wp%ky0 * y), dp)
        norm = norm + abs(psi0(idx))**2 * dx * dy
      end do
    end do
    psi0 = psi0 / sqrt(norm)
  end subroutine init_wavepacket

  subroutine project_onto_eigenbasis(eigvecs, psi0, coeffs)
    complex(dp), intent(in) :: eigvecs(n, n), psi0(n)
    complex(dp), intent(out) :: coeffs(n)
    call zgemv('C', n, n, (1.0_dp, 0.0_dp), eigvecs, n, psi0, 1, (0.0_dp, 0.0_dp), coeffs, 1)
  end subroutine project_onto_eigenbasis

  subroutine density_at_time(eigvecs, eigvals, coeffs, t, density)
    complex(dp), intent(in) :: eigvecs(n, n), coeffs(n)
    real(dp), intent(in) :: eigvals(n), t
    real(dp), intent(out) :: density(n)
    complex(dp) :: phased(n), psit(n)
    integer :: k
! bleh :3
    do k = 1, n
      phased(k) = coeffs(k) * cmplx(cos(-eigvals(k) * t), sin(-eigvals(k) * t), dp)
    end do
    call zgemv('N', n, n, (1.0_dp, 0.0_dp), eigvecs, n, phased, 1, (0.0_dp, 0.0_dp), psit, 1)
    density = abs(psit)**2
  end subroutine density_at_time

end module evolve
