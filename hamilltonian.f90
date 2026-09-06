module hamiltonian
  use kinds
  use geometry
  implicit none

contains

  subroutine build_hamiltonian(H)
    complex(dp), intent(out) :: H(n, n)
    real(dp) :: kin, x, y
    integer :: i, j, idx, nb

    H = (0.0_dp, 0.0_dp)
    kin = hbar * hbar / (2.0_dp * mass)

    do j = 1, ny
      do i = 1, nx
        idx = grid_index(i, j)
        x = grid_x(i)
        y = grid_y(j)
        H(idx, idx) = 2.0_dp * kin * (1.0_dp / dx**2 + 1.0_dp / dy**2) + potential(x, y)

        if (i > 1) then
          nb = grid_index(i - 1, j)
          H(idx, nb) = H(idx, nb) - kin / dx**2
        end if
        if (i < nx) then
          nb = grid_index(i + 1, j)
          H(idx, nb) = H(idx, nb) - kin / dx**2
        end if
        if (j > 1) then
          nb = grid_index(i, j - 1)
          H(idx, nb) = H(idx, nb) - kin / dy**2
        end if
        if (j < ny) then
          nb = grid_index(i, j + 1)
          H(idx, nb) = H(idx, nb) - kin / dy**2
        end if
      end do
    end do
  end subroutine build_hamiltonian

  subroutine diagonalize(H, eigvecs, eigvals)
    complex(dp), intent(inout) :: H(n, n)
    complex(dp), intent(out) :: eigvecs(n, n)
    real(dp), intent(out) :: eigvals(n)
    complex(dp), allocatable :: work(:)
    real(dp) :: rwork(max(1, 3 * n - 2))
    integer :: lwork, info

    eigvecs = H
    lwork = -1
    allocate(work(1))
    call zheev('V', 'U', n, eigvecs, n, eigvals, work, lwork, rwork, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))
    call zheev('V', 'U', n, eigvecs, n, eigvals, work, lwork, rwork, info)
    if (info /= 0) error stop 'zheev failed'
  end subroutine diagonalize

end module hamiltonian
