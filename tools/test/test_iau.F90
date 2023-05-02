!-----------------------------------------------------------------------
!
    program test_iau
!
!-----------------------------------------------------------------------
!***  This driver is designed to test reading a cubed sphere inc file
!-----------------------------------------------------------------------
!***
!***  Revision history
!***
!     Oct 2022:  M. Potts             - initial code for fv3 read
!
!---------------------------------------------------------------------------------
!
      use mpi, only: MPI_Init, MPI_Finalize
      use atmosphere_mod,     only: Atm, mygrid

      use module_get_cubed_sphere_inc,  only : read_netcdf_incs
      use CCPP_data,          only: GFS_control
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!
!
      integer              :: rc, grid_id, mype, mpi_comm
      character(len=128), dimension(3) :: filenames
!     type(iau_external_data_type)        :: IAU_Data ! number of blocks
!     type(GFS_init_type)  :: Init_parm
      logical :: tests_passed, use_parallel_netcdf, testing

     
      filenames = (/'temp0.nc4', 'temp1.nc4', 'temp2.nc4'/)
      use_parallel_netcdf=.false.
      testing=.true.
      mpi_comm = 0
      mype = 0
      grid_id = 1
!     call MPI_Init(rc)
      mpi_comm=0
      mype=0
      GFS_control%isc=1
      GFS_control%nx=96
      GFS_control%jsc=1
      GFS_control%ny=96
!     call iau_initialize_netcdf(GFS_control,IAU_data,Init_parm)
!     call iau_initialize(GFS_control)
      call read_netcdf_incs(filenames, Atm, mygrid, &
                       use_parallel_netcdf, &
                       testing, tests_passed=tests_passed, rc=rc)

!     call MPI_Finalize(rc)
      if(tests_passed) then
        call exit(0)
      else
        call exit(1)
      endif
     end program
