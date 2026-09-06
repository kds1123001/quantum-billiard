module render
  use kinds
  use geometry
  use iso_c_binding, only: c_int32_t
  implicit none

  interface
    subroutine usleep(useconds) bind(C, name="usleep")
      import :: c_int32_t
      integer(c_int32_t), value :: useconds
    end subroutine usleep
  end interface

contains

  subroutine cursor_home()
    write(*, '(A)', advance='no') achar(27) // '[H'
  end subroutine cursor_home

  subroutine draw_frame(density, maxval_global, frame_idx, nframes, t)
    real(dp), intent(in) :: density(n), maxval_global, t
    integer, intent(in) :: frame_idx, nframes
    integer :: i, j, idx, r, g, b, pos, clen
    real(dp) :: v, tt, x, y
    character(len=4200) :: line
    character(len=32) :: cell
    character(len=1) :: esc

    esc = achar(27)

    call cursor_home()

    do j = 1, ny
      pos = 1
      do i = 1, nx
        idx = grid_index(i, j)
        x = grid_x(i)
        y = grid_y(j)
        v = density(idx)
        tt = min(v / maxval_global, 1.0_dp) ** 0.5_dp

        if (inside_obstacle(x, y) .and. tt < 0.08_dp) then
          cell = esc // '[38;2;90;90;110m' // achar(226) // achar(150) // achar(146)
        else
          r = int(20 + 235 * tt)
          g = int(20 + 130 * tt**1.6_dp)
          b = int(60 + 160 * (1.0_dp - tt) * tt + 40)
          write(cell, '(A,I0,A,I0,A,I0,A)') esc // '[38;2;', r, ';', g, ';', b, &
               'm' // achar(226) // achar(150) // achar(136)
        end if
        clen = len_trim(cell)
        line(pos:pos+clen-1) = cell(1:clen)
        pos = pos + clen
      end do
      line(pos:pos+3) = esc // '[0m'
      pos = pos + 4
      write(*, '(A)') line(1:pos-1)
    end do

    write(*, '(A,I5,A,I5,A,F8.3,A)') 'frame ', frame_idx, ' / ', nframes, '   t = ', t, '        '
  end subroutine draw_frame

  subroutine sleep_ms(ms)
    integer, intent(in) :: ms
    call usleep(int(ms, c_int32_t) * 1000_c_int32_t)
  end subroutine sleep_ms ! lol cat!

end module render
