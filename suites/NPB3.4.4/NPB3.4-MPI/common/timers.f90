!---------------------------------------------------------------------
! NPB-MPI Fortran timers — bind(C) to tacvar_npb_* (shared with C path).
!---------------------------------------------------------------------

      module tacvar_npb_c
        use, intrinsic :: iso_c_binding
        implicit none
        interface
          subroutine tacvar_npb_timer_clear(n) bind(C, name="tacvar_npb_timer_clear")
            import :: c_int
            integer(c_int), value :: n
          end subroutine
          subroutine tacvar_npb_timer_start(n) bind(C, name="tacvar_npb_timer_start")
            import :: c_int
            integer(c_int), value :: n
          end subroutine
          subroutine tacvar_npb_timer_stop(n) bind(C, name="tacvar_npb_timer_stop")
            import :: c_int
            integer(c_int), value :: n
          end subroutine
          function tacvar_npb_timer_read(n) bind(C, name="tacvar_npb_timer_read")
            import :: c_int, c_double
            integer(c_int), value :: n
            real(c_double) :: tacvar_npb_timer_read
          end function
          subroutine tacvar_npb_prepare(expected_steps)  &
     &         bind(C, name="tacvar_npb_prepare")
            import :: c_int64_t
            integer(c_int64_t), value :: expected_steps
          end subroutine
          subroutine tacvar_npb_step_start() bind(C, name="tacvar_npb_step_start")
          end subroutine
          subroutine tacvar_npb_step_stop() bind(C, name="tacvar_npb_step_stop")
          end subroutine
          subroutine tacvar_npb_kernel_start() bind(C, name="tacvar_npb_kernel_start")
          end subroutine
          subroutine tacvar_npb_kernel_stop() bind(C, name="tacvar_npb_kernel_stop")
          end subroutine
        end interface
      end module tacvar_npb_c

      subroutine timer_clear(n)
      use tacvar_npb_c
      use, intrinsic :: iso_c_binding, only: c_int
      implicit none
      integer n
      call tacvar_npb_timer_clear(int(n, kind=c_int))
      return
      end

      subroutine timer_start(n)
      use tacvar_npb_c
      use, intrinsic :: iso_c_binding, only: c_int
      implicit none
      integer n
      call tacvar_npb_timer_start(int(n, kind=c_int))
      return
      end

      subroutine timer_stop(n)
      use tacvar_npb_c
      use, intrinsic :: iso_c_binding, only: c_int
      implicit none
      integer n
      call tacvar_npb_timer_stop(int(n, kind=c_int))
      return
      end

      double precision function timer_read(n)
      use tacvar_npb_c
      use, intrinsic :: iso_c_binding, only: c_int
      implicit none
      integer n
      timer_read = tacvar_npb_timer_read(int(n, kind=c_int))
      return
      end

!---------------------------------------------------------------------
! TacVar step helpers (instrumentation only; no workload change)
!---------------------------------------------------------------------
      subroutine tacvar_prepare(expected_steps)
      use tacvar_npb_c
      use, intrinsic :: iso_c_binding, only: c_int64_t
      implicit none
      integer(kind=8) expected_steps
      call tacvar_npb_prepare(int(expected_steps, kind=c_int64_t))
      return
      end

      subroutine tacvar_step_start()
      use tacvar_npb_c
      implicit none
      call tacvar_npb_step_start()
      return
      end

      subroutine tacvar_step_stop()
      use tacvar_npb_c
      implicit none
      call tacvar_npb_step_stop()
      return
      end

      subroutine check_timer_flag(timeron)
      implicit none
      logical timeron
      integer nc, ios
      character(len=20) val

      timeron = .false.
      call get_environment_variable('NPB_TIMER_FLAG', val, nc, ios)
      if (ios .eq. 0) then
         if (nc .le. 0) then
            timeron = .true.
         else if (val(1:1) .ge. '1' .and. val(1:1) .le. '9') then
            timeron = .true.
         else if (val .eq. 'on' .or. val .eq. 'ON' .or.  &
     &            val .eq. 'yes' .or. val .eq. 'YES' .or.  &
     &            val .eq. 'true' .or. val .eq. 'TRUE') then
            timeron = .true.
         endif
      else
         open(unit=99, file='timer.flag', status='old', iostat=ios)
         if (ios .eq. 0) then
            close(99)
            timeron = .true.
         endif
      endif
      return
      end
