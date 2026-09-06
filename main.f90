program qbilliard_term
  use kinds
  use geometry
  use hamiltonian
  use evolve
  use render
  implicit none

  complex(dp), allocatable :: H(:,:), eigvecs(:,:), psi0(:), coeffs(:)
  real(dp), allocatable :: eigvals(:), density(:)
  type(wavepacket_t) :: wp
  integer :: preset, nframes, delay_ms, f
  real(dp) :: dt, t, maxd, frame_max
  character(len=1) :: again
  integer(kind=8) :: clock_start, clock_end, clock_rate

  allocate(H(n, n), eigvecs(n, n), eigvals(n))
  allocate(psi0(n), coeffs(n), density(n))

  print '(A)', 'quantum billiards -- single particle, exact eigenbasis propagation'
  print '(A,I0,A,I0,A,I0,A)', 'grid ', nx, 'x', ny, ' (', n, ' basis states)'
  print '(A)', ''

  do
    call read_int('obstacle layout  [0]=empty [1]=single bumper [2]=five bumpers', 2, preset)
    call set_obstacles(preset)
    print '(A)', ''

    call system_clock(clock_start, clock_rate)
    call build_hamiltonian(H)
    call diagonalize(H, eigvecs, eigvals)
    call system_clock(clock_end)
    print '(A,F6.2,A)', 'diagonalized in ', real(clock_end - clock_start, dp) / real(clock_rate, dp), ' s'
    print '(A,F10.5)', 'ground state energy: ', eigvals(1)
    print '(A)', ''

    call read_real('wavepacket start x  (0 - ' // trim(to_str(lx)) // ')', 4.0_dp, wp%x0)
    call read_real('wavepacket start y  (0 - ' // trim(to_str(ly)) // ')', 5.0_dp, wp%y0)
    call read_real('momentum kx', 6.0_dp, wp%kx0)
    call read_real('momentum ky', 0.5_dp, wp%ky0)
    call read_real('wavepacket width sigma', 0.7_dp, wp%sigma)
    call read_int('number of frames', 400, nframes)
    call read_real('time step dt', 0.006_dp, dt)
    call read_int('frame delay (ms)', 35, delay_ms)

    call init_wavepacket(wp, psi0)
    call project_onto_eigenbasis(eigvecs, psi0, coeffs)

    maxd = 0.0_dp
    do f = 1, nframes
      t = (f - 1) * dt
      call density_at_time(eigvecs, eigvals, coeffs, t, density)
      frame_max = maxval(density)
      if (frame_max > maxd) maxd = frame_max
    end do

    write(*, '(A)', advance='no') achar(27) // '[2J'
    do f = 1, nframes
      t = (f - 1) * dt
      call density_at_time(eigvecs, eigvals, coeffs, t, density)
      call draw_frame(density, maxd, f, nframes, t)
      call sleep_ms(delay_ms)
    end do

    print '(A)', ''
    call read_char('run again with new parameters? [y/n]', 'y', again)
    if (again == 'n' .or. again == 'N') exit
    print '(A)', ''
  end do

contains

  subroutine read_real(prompt, default, value)
    character(len=*), intent(in) :: prompt
    real(dp), intent(in) :: default
    real(dp), intent(out) :: value
    character(len=200) :: line
    integer :: ios
    write(*, '(A,A,F0.4,A)', advance='no') trim(prompt), ' [', default, ']: '
    read(*, '(A)', iostat=ios) line
    value = default
    if (ios == 0 .and. len_trim(line) > 0) then
      read(line, *, iostat=ios) value
      if (ios /= 0) value = default
    end if
  end subroutine read_real

  subroutine read_int(prompt, default, value)
    character(len=*), intent(in) :: prompt
    integer, intent(in) :: default
    integer, intent(out) :: value
    character(len=200) :: line
    integer :: ios
    write(*, '(A,A,I0,A)', advance='no') trim(prompt), ' [', default, ']: '
    read(*, '(A)', iostat=ios) line
    value = default
    if (ios == 0 .and. len_trim(line) > 0) then
      read(line, *, iostat=ios) value
      if (ios /= 0) value = default
    end if
  end subroutine read_int

  subroutine read_char(prompt, default, value)
    character(len=*), intent(in) :: prompt
    character(len=1), intent(in) :: default
    character(len=1), intent(out) :: value
    character(len=200) :: line
    integer :: ios
    write(*, '(A,A,A,A)', advance='no') trim(prompt), ' [', default, ']: '
    read(*, '(A)', iostat=ios) line
    value = default
    if (ios == 0 .and. len_trim(line) > 0) value = line(1:1)
  end subroutine read_char

  function to_str(v) result(s)
    real(dp), intent(in) :: v
    character(len=16) :: s
    write(s, '(F0.1)') v
  end function to_str

end program qbilliard_term
