#define NC_ERR_STOP(status) \
    if (status /= nf90_noerr) write(0,*) "file: ", __FILE__, " line: ", __LINE__, trim(nf90_strerror(status)); \
    if (status /= nf90_noerr) call ESMF_Finalize(endflag=ESMF_END_ABORT)

module module_get_cubed_sphere_inc

  use netcdf
  use esmf
  use mpp_mod,            only: mpp_pe, mpp_get_current_pelist
  use tracer_manager_mod, only: get_tracer_index, get_number_tracers
  use field_manager_mod,  only: MODEL_ATMOS
  use CCPP_data,          only: GFS_control
  use time_manager_mod,   only: time_type, get_time, set_time
  use fv_arrays_mod,      only: fv_atmos_type
  use mpi

  implicit none
  type iau_internal_data_type
    real,allocatable :: ua_inc(:,:,:)
    real,allocatable :: va_inc(:,:,:)
    real,allocatable :: temp_inc(:,:,:)
    real,allocatable :: delp_inc(:,:,:)
    real,allocatable :: delz_inc(:,:,:)
    real,allocatable :: tracer_inc(:,:,:,:)
  end type iau_internal_data_type

  public read_netcdf, read_netcdf_incs, read_netcdf_inc, iau_internal_data_type

  logical :: par

  contains

!----------------------------------------------------------------------------------------
  logical function compvals(a, b)
    real, intent(in)     :: a
    real, intent(in)     :: b

    if(abs(a - b) < 1.0e-20) then
       compvals = .true.
    else
       compvals = .false.
    endif
    write(6,*) 'HEY!!! comp is ',abs(a - b),compvals

  end function
  subroutine read_netcdf_inc(filename, increment_data, Atm, mygrid, &
                          testing, im_ret, jm_ret, pf_ret, tileCount, tests_passed,rc)
    character(*), intent(in)                         :: filename
    type(iau_internal_data_type), intent(inout)      :: increment_data
    type (fv_atmos_type), allocatable, intent(inout) :: Atm(:)
    integer, intent(in)                              :: mygrid
    logical, intent(in)                              :: testing
    integer, optional, intent(out)                   :: im_ret
    integer, optional, intent(out)                   :: jm_ret
    integer, optional, intent(out)                   :: pf_ret
    integer, optional, intent(out)                   :: tileCount
    logical, optional,intent(out)                    :: tests_passed
    integer, optional,intent(out)                    :: rc

! local variables
    integer :: i,j
    integer :: mpi_comm
    integer :: mype

    real(ESMF_KIND_R8), dimension(:,:,:,:,:), pointer :: array_r8_3d_tiled
    real(ESMF_KIND_R8)  :: testval
    integer :: ncerr
    integer :: ncid
    integer :: im_dimid, jm_dimid, tile_dimid, pfull_dimid, phalf_dimid, time_dimid 
    integer :: tm,ph
    integer :: im_varid, jm_varid, tile_varid, lon_varid, lat_varid, timeiso_varid
    integer :: liqwat_varid, sphum_varid, o3mr_varid, icewat_varid
    integer :: ugrd_varid, vgrd_varid, dpres_varid, tmp_varid 

    integer :: par_access, nvar, xtype, ndims, nAtts 
    integer :: sphum_idx, liq_wat_idx, ice_wat_idx, o3mr_idx
    integer :: isc, iec, jsc, jec, TC
    integer :: im, jm, pf
    integer :: mytile
    integer, dimension(:), allocatable :: varids
    character(len=NF90_MAX_NAME) :: varname
!
    TC = 6
    write(6,*) 'in read_netcdf_inc'
    if(present(tileCount)) TC = tileCount
    if(present(im_ret)) im = im_ret
    if(present(jm_ret)) jm = jm_ret
    if(present(pf_ret)) pf = pf_ret
    if(present(tests_passed)) tests_passed = .true.
    if(testing) then
      testval = 0.0
      mytile = 1
      mype = 0
      mpi_comm = 0
    else 
      testval = 1.0
      mytile = Atm(mygrid)%tile_of_mosaic
      mype = mpp_pe()
      call mpp_get_current_pelist(Atm(mygrid)%pelist, commID=mpi_comm)
    endif
    par = .false.

    if (par) then
       ! not implemented yet
       ncerr = nf90_open(trim(filename),&
               mode=nf90_nowrite, &
               comm=mpi_comm, info = MPI_INFO_NULL, ncid=ncid); NC_ERR_STOP(ncerr)
    else
       ncerr = nf90_open(trim(filename),&
               mode=nf90_nowrite, &
               ncid=ncid); NC_ERR_STOP(ncerr)
    end if
    ncerr = nf90_inquire(ncid, nvariables = nvar); NC_ERR_STOP(ncerr)
    write(6,*) 'nvars is ',nvar
    allocate(varids(nvar))
    ncerr = nf90_inq_varids(ncid, nvar, varids); NC_ERR_STOP(ncerr)
    do i=1,nvar
      ncerr = nf90_inquire_variable(ncid, varids(i), name=varname, xtype=xtype, ndims=ndims, nAtts=nAtts)
      NC_ERR_STOP(ncerr)
      write(6,*) 'Name is ',trim(varname), varids(i),ndims
    enddo
    !get dimensions of fields in the file
      
    varname = "grid_yt"
    ncerr = nf90_inq_dimid(ncid, trim(varname), jm_dimid) ; NC_ERR_STOP(ncerr)
    ncerr = nf90_inq_varid(ncid,trim(varname),jm_varid); NC_ERR_STOP(ncerr)
    ncerr = nf90_inquire_dimension(ncid,jm_dimid,len=jm); NC_ERR_STOP(ncerr)
    write(6,*) 'got jm ',jm
    varname = "grid_xt"
    ncerr = nf90_inq_dimid(ncid, trim(varname), im_dimid) ; NC_ERR_STOP(ncerr)
    write(6,*) 'im_dimid', ncerr, im_dimid
    ncerr = nf90_inq_varid(ncid,trim(varname),im_varid); NC_ERR_STOP(ncerr)
    write(6,*) 'im_varid', ncerr, im_varid
    ncerr = nf90_inquire_dimension(ncid,im_dimid,len=im); NC_ERR_STOP(ncerr)
    write(6,*) 'ncerr', ncerr
    write(6,*) 'got im ',im
    varname = "time"
    ncerr = nf90_inq_dimid(ncid, trim(varname), time_dimid) ; NC_ERR_STOP(ncerr)
    ncerr = nf90_inquire_dimension(ncid,time_dimid,len=tm); NC_ERR_STOP(ncerr)
    varname = "tile"
    write(6,*) 'got tm ',tm
    ncerr = nf90_inq_dimid(ncid, trim(varname), tile_dimid) ; NC_ERR_STOP(ncerr)
    ncerr = nf90_inquire_dimension(ncid,tile_dimid,len=TC); NC_ERR_STOP(ncerr)
    varname = "pfull"
    write(6,*) 'got TC ',TC
    ncerr = nf90_inq_dimid(ncid, trim(varname), pfull_dimid) ; NC_ERR_STOP(ncerr)
    ncerr = nf90_inquire_dimension(ncid,pfull_dimid,len=pf); NC_ERR_STOP(ncerr)
    varname = "phalf"
    write(6,*) 'got pf ',pf
    ncerr = nf90_inq_dimid(ncid, trim(varname), phalf_dimid) ; NC_ERR_STOP(ncerr)
    ncerr = nf90_inquire_dimension(ncid,phalf_dimid,len=ph); NC_ERR_STOP(ncerr)

    if(present(im_ret)) im_ret = im
    if(present(jm_ret)) jm_ret = jm
    if(present(pf_ret)) pf_ret = pf
    if(present(tileCount)) tileCount = TC
    !get the variable id's we will need for each variable to be retrieved
    !TODO? put this in an indexed array?
    varname = "ice_wat"
    ncerr = nf90_inq_varid(ncid,trim(varname),icewat_varid); NC_ERR_STOP(ncerr); NC_ERR_STOP(ncerr)
    varname = "liq_wat"
    ncerr = nf90_inq_varid(ncid,trim(varname),liqwat_varid); NC_ERR_STOP(ncerr); NC_ERR_STOP(ncerr)
    varname = "spfh"
    ncerr = nf90_inq_varid(ncid,trim(varname),sphum_varid); NC_ERR_STOP(ncerr); NC_ERR_STOP(ncerr)
    varname = "o3mr"
    ncerr = nf90_inq_varid(ncid,trim(varname),o3mr_varid); NC_ERR_STOP(ncerr)
    varname = "ugrd"
    ncerr = nf90_inq_varid(ncid,trim(varname),ugrd_varid); NC_ERR_STOP(ncerr)
    varname = "vgrd"
    ncerr = nf90_inq_varid(ncid,trim(varname),vgrd_varid); NC_ERR_STOP(ncerr)
    varname = "dpres"
    ncerr = nf90_inq_varid(ncid,trim(varname),dpres_varid); NC_ERR_STOP(ncerr)
    varname = "tmp"
    ncerr = nf90_inq_varid(ncid,trim(varname),tmp_varid); NC_ERR_STOP(ncerr)
    write(6,*) 'at allocating step',im,jm,pf,TC,tm
    if(.not. allocated(increment_data%ua_inc)) then
      ! Allocate space in increment 
      allocate(increment_data%ua_inc(im,jm,pf))
      allocate(increment_data%va_inc(im,jm,pf))
      allocate(increment_data%temp_inc(im,jm,pf))
      allocate(increment_data%delp_inc(im,jm,pf))
      allocate(increment_data%tracer_inc(im,jm,pf,4)) ! currently, just need specific humdity, liq_wat, ice_wat, and o3mr
    endif
    if(testing) then
      ! assign dummy indices
      sphum_idx   = 1
      liq_wat_idx = 2
      ice_wat_idx = 3
      o3mr_idx    = 4
    else
      ! get the tracer index to update variables in Atm(t)%tracer_inc
      ! this in not available when on the write component
      write(6,*) 'getting indices'
      sphum_idx   = get_tracer_index(MODEL_ATMOS, 'sphum')
      liq_wat_idx = get_tracer_index(MODEL_ATMOS, 'liq_wat')
      ice_wat_idx = get_tracer_index(MODEL_ATMOS, 'ice_wat')
      o3mr_idx    = get_tracer_index(MODEL_ATMOS, 'o3mr')
    endif
    ! allocate temporary array to hold variables
    ! TODO read only what we need instead of the whole field
    allocate(array_r8_3d_tiled(im,jm,pf,TC,tm))
    isc = GFS_control%isc
    iec = GFS_control%isc+GFS_control%nx-1
    jsc = GFS_control%jsc
    jec = GFS_control%jsc+GFS_control%ny-1
    !Read u
    ncerr = nf90_get_var(ncid, ugrd_varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
    increment_data%ua_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)
    write(6,*) 'reading u',increment_data%ua_inc(isc,jsc,1) 
    
    !Read v
    ncerr = nf90_get_var(ncid, vgrd_varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
    increment_data%va_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)

    !Read potential temp
    ncerr = nf90_get_var(ncid, tmp_varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
    increment_data%temp_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)

    !Read delp
    ncerr = nf90_get_var(ncid, dpres_varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
    increment_data%delp_inc(isc:iec,jsc:jec,:) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)

    !Update sphum
    ncerr = nf90_get_var(ncid, sphum_varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
    increment_data%tracer_inc(isc:iec,jsc:jec,:,sphum_idx) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)

    !Update ice_wat
    ncerr = nf90_get_var(ncid, icewat_varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
    increment_data%tracer_inc(isc:iec,jsc:jec,:,ice_wat_idx) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)

    !Update liq_wat
    ncerr = nf90_get_var(ncid, liqwat_varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
    increment_data%tracer_inc(isc:iec,jsc:jec,:,liq_wat_idx) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)

    !Update o3mr 
    ncerr = nf90_get_var(ncid, o3mr_varid, array_r8_3d_tiled) ; NC_ERR_STOP(ncerr)
    increment_data%tracer_inc(isc:iec,jsc:jec,:,o3mr_idx) = array_r8_3d_tiled(isc:iec,jsc:jec,:,mytile,1)

    deallocate(array_r8_3d_tiled)
  
  end subroutine read_netcdf_inc
  subroutine read_netcdf(filename, Atm, mygrid, &
                          use_parallel_netcdf, &
                          testing,tests_passed,rc)
!
    character(*), intent(in)                         :: filename
    type (fv_atmos_type), allocatable, intent(inout) :: Atm(:)
    integer,          intent(in)                     :: mygrid
    logical, intent(in)                              :: use_parallel_netcdf
    logical, intent(in)                              :: testing
    logical, optional,intent(out)                    :: tests_passed
    integer, optional,intent(out)                    :: rc
!
!** local vars
    type(iau_internal_data_type)               :: increment_data
    integer :: i,j,k
    integer :: im, jm
    integer :: mpi_comm
    integer :: mype

    real(ESMF_KIND_R8), dimension(:,:,:,:,:), pointer :: array_r8_3d_tiled
    real(ESMF_KIND_R8)  :: testval

    integer :: ncerr
    integer :: ncid
    integer :: im_dimid, jm_dimid, tile_dimid, pfull_dimid, phalf_dimid, time_dimid 
    integer :: tm,pf,ph, tileCount
    integer :: im_varid, jm_varid, tile_varid, lon_varid, lat_varid, timeiso_varid
    integer :: liqwat_varid, sphum_varid, o3mr_varid, icewat_varid
    integer :: ugrd_varid, vgrd_varid, dpres_varid, tmp_varid 

    integer :: par_access, nvar, xtype, ndims, nAtts 
    integer :: sphum_idx, liq_wat_idx, ice_wat_idx, o3mr_idx
    integer :: isc, iec, jsc, jec
    integer :: mytile
    integer, dimension(:), allocatable :: varids
    character(len=NF90_MAX_NAME) :: varname

    if(testing) then
      testval = 0.0
      mytile = 1
      tileCount = 6
      mype = 0
      mpi_comm = 0
      allocate(increment_data%ua_inc(96,96,127))
      allocate(increment_data%va_inc(96,96,127))
      allocate(increment_data%temp_inc(96,96,127))
      allocate(increment_data%delp_inc(96,96,127))
      allocate(increment_data%delz_inc(96,96,127))
      allocate(increment_data%tracer_inc(96,96,127,6))
    else 
      testval = 1.0
      mytile = Atm(mygrid)%tile_of_mosaic
      mype = mpp_pe()
      call mpp_get_current_pelist(Atm(mygrid)%pelist, commID=mpi_comm)
    endif
    par = .false.
    call read_netcdf_inc(filename, increment_data,Atm, mygrid, testing,im_ret=im, jm_ret=jm, pf_ret=pf, tileCount=tileCount, tests_passed=tests_passed, rc=rc)
    if(testing) then
      ! allocate 6 tiles for Atm
      write(6,*) "im, jm, etc are",im,jm,pf
      allocate(Atm(6))
      ! assign dummy indices
      sphum_idx   = 1
      liq_wat_idx = 2
      ice_wat_idx = 3
      o3mr_idx    = 4
      ! Allocate space in Atm for testing 
      do i=1,tileCount
        allocate(Atm(i)%u(im,jm,pf))
        allocate(Atm(i)%v(im,jm,pf))
        allocate(Atm(i)%pt(im,jm,pf))
        allocate(Atm(i)%delp(im,jm,pf))
        allocate(Atm(i)%q(im,jm,pf,4)) ! currently, just need specific humdity, liq_wat, ice_wat, and o3mr
        Atm(i)%u(:,:,:) = 0.0
        Atm(i)%v(:,:,:) = 0.0
        Atm(i)%pt(:,:,:) = 0.0
        Atm(i)%delp(:,:,:) = 0.0
        Atm(i)%q(:,:,:,:) = 0.0
      enddo
    else
      ! get the tracer index to update variables in Atm(t)%tracer_inc
      ! this in not available when on the write component
      sphum_idx   = get_tracer_index(MODEL_ATMOS, 'sphum')
      liq_wat_idx = get_tracer_index(MODEL_ATMOS, 'liq_wat')
      ice_wat_idx = get_tracer_index(MODEL_ATMOS, 'ice_wat')
      o3mr_idx    = get_tracer_index(MODEL_ATMOS, 'o3mr')
    endif
    isc = GFS_control%isc
    iec = GFS_control%isc+GFS_control%nx-1
    jsc = GFS_control%jsc
    jec = GFS_control%jsc+GFS_control%ny-1
    !Update u
    Atm(mygrid)%u(isc:iec,jsc:jec,:) = Atm(mygrid)%u(isc:iec,jsc:jec,:) + &
      testval * increment_data%ua_inc(isc:iec,jsc:jec,:)
    
    !Update v
    Atm(mygrid)%v(isc:iec,jsc:jec,:) = Atm(mygrid)%v(isc:iec,jsc:jec,:) + &
      testval * increment_data%va_inc(isc:iec,jsc:jec,:)

    !Update potential temp
    Atm(mygrid)%pt(isc:iec,jsc:jec,:) = Atm(mygrid)%pt(isc:iec,jsc:jec,:) + &
      testval * increment_data%temp_inc(isc:iec,jsc:jec,:)

    !Update delp
    Atm(mygrid)%delp(isc:iec,jsc:jec,:) = Atm(mygrid)%delp(isc:iec,jsc:jec,:) + &
      testval * increment_data%delp_inc(isc:iec,jsc:jec,:)

    !Update sphum
    Atm(mygrid)%q(isc:iec,jsc:jec,:,sphum_idx) = Atm(mygrid)%q(isc:iec,jsc:jec,:,sphum_idx) + &
      testval * increment_data%tracer_inc(isc:iec,jsc:jec,:,sphum_idx)

    !Update ice_wat
    Atm(mygrid)%q(isc:iec,jsc:jec,:,ice_wat_idx) = Atm(mygrid)%q(isc:iec,jsc:jec,:,ice_wat_idx) + &
      testval * increment_data%tracer_inc(isc:iec,jsc:jec,:,ice_wat_idx)

    !Update liq_wat
    Atm(mygrid)%q(isc:iec,jsc:jec,:,liq_wat_idx) = Atm(mygrid)%q(isc:iec,jsc:jec,:,liq_wat_idx) + &
      testval * increment_data%tracer_inc(isc:iec,jsc:jec,:,liq_wat_idx)
    !Update o3mr 
    Atm(mygrid)%q(isc:iec,jsc:jec,:,o3mr_idx) = Atm(mygrid)%q(isc:iec,jsc:jec,:,o3mr_idx) + &
      testval * increment_data%tracer_inc(isc:iec,jsc:jec,:,o3mr_idx)

    if(testing) then
        tests_passed = tests_passed .and. compvals(increment_data%ua_inc(isc,jsc,1) ,-1.8169951587765354E-007)
        tests_passed = tests_passed .and. compvals(increment_data%va_inc(isc,jsc,1) , 3.0927015537418612E-007) 
        tests_passed = tests_passed .and. compvals(increment_data%delp_inc(isc,jsc,1), -1.8873791418627661E-015)
        tests_passed = tests_passed .and. compvals(increment_data%temp_inc(isc,jsc,1) , -6.8499218741635559E-008)
        tests_passed = tests_passed .and. compvals(increment_data%tracer_inc(isc,jsc,1,ice_wat_idx) , -4.4666691960233954E-019)
        tests_passed = tests_passed .and. compvals(increment_data%tracer_inc(isc,jsc,1,liq_wat_idx) , -7.3514147181582930E-022)
        tests_passed = tests_passed .and. compvals(increment_data%tracer_inc(isc,jsc,1,o3mr_idx) , -1.2690141736645521E-017)
        tests_passed = tests_passed .and. compvals(increment_data%tracer_inc(isc,jsc,1,sphum_idx) , -2.9569591390493574E-006)

      do i=1,tileCount
        deallocate(Atm(i)%u)
        deallocate(Atm(i)%v)
        deallocate(Atm(i)%pt)
        deallocate(Atm(i)%delp)
        deallocate(Atm(i)%q) ! currently, just need specific humdity, liq_wat, ice_wat, and o3mr
      enddo
      deallocate(Atm) ! currently, just need specific humdity, liq_wat, ice_wat, and o3mr
    endif 
      
  end subroutine read_netcdf
  subroutine read_netcdf_incs(filenames, &
                          Atm, mygrid, &
                          use_parallel_netcdf, &
                          testing, tests_passed, rc)
!
    character(len=128), dimension(3), intent(in) :: filenames
    type (fv_atmos_type), allocatable, intent(inout) :: Atm(:)
    integer,                             intent(in) :: mygrid
    logical, intent(in)                             :: use_parallel_netcdf
    logical, intent(in)                             :: testing
    logical, intent(out),optional                   :: tests_passed
    integer, intent(out),optional                   :: rc
!
!** local vars
    integer :: mpi_comm
    integer :: mype
    type(iau_internal_data_type)               :: increment_data(3)
    type(iau_internal_data_type)               :: tendency
    integer :: i,j,k
    integer :: im, jm

    real(ESMF_KIND_R8), dimension(:,:,:,:,:), pointer :: array_r8_3d_tiled
    real(ESMF_KIND_R8)  :: testval

    integer :: ncerr
    integer :: ncid
    integer :: im_dimid, jm_dimid, tile_dimid, pfull_dimid, phalf_dimid, time_dimid 
    integer :: tm,pf,ph, tileCount
    integer :: im_varid, jm_varid, tile_varid, lon_varid, lat_varid, timeiso_varid
    integer :: liqwat_varid, sphum_varid, o3mr_varid, icewat_varid
    integer :: ugrd_varid, vgrd_varid, dpres_varid, tmp_varid 

    integer :: par_access, nvar, xtype, ndims, nAtts 
    integer :: sphum_idx, liq_wat_idx, ice_wat_idx, o3mr_idx
    integer :: isc, iec, jsc, jec
    integer :: mytile
    real :: dt_atmos
    integer, dimension(:), allocatable :: varids
    character(len=NF90_MAX_NAME) :: varname

    dt_atmos = 250.
    if(testing) then
      testval = 0.0
      mytile = 1
      mype = 0
      mpi_comm = 0
    else 
      testval = 1.0
      mytile = Atm(mygrid)%tile_of_mosaic
      mype = mpp_pe()
      call mpp_get_current_pelist(Atm(mygrid)%pelist, commID=mpi_comm)
    endif
    par = .false.
    ! read in 3? increments
    do i=1,3
      call read_netcdf_inc(filenames(i), increment_data(i),Atm,mygrid,testing,im_ret=im, jm_ret=jm, pf_ret=pf, tileCount=tileCount, tests_passed=tests_passed, rc=rc)
    enddo
    ! allocate space for tendencies
    allocate(tendency%ua_inc(im,jm,pf))
    allocate(tendency%va_inc(im,jm,pf))
    allocate(tendency%temp_inc(im,jm,pf))
    allocate(tendency%delp_inc(im,jm,pf))
    allocate(tendency%tracer_inc(im,jm,pf,4)) ! currently, just need specific humdity, liq_wat, ice_wat, and o3mr
    ! Calculate the tendencies by subtracting the first increment from the second and dividing by dt_atmos
    tendency%ua_inc(:,:,:) = (increment_data(2)%ua_inc(:,:,:) - increment_data(1)%ua_inc(:,:,:))/dt_atmos
    tendency%va_inc(:,:,:) = (increment_data(2)%va_inc(:,:,:) - increment_data(1)%va_inc(:,:,:))/dt_atmos
    tendency%temp_inc(:,:,:) = (increment_data(2)%temp_inc(:,:,:) - increment_data(1)%temp_inc(:,:,:))/dt_atmos
    tendency%delp_inc(:,:,:) = (increment_data(2)%delp_inc(:,:,:) - increment_data(1)%delp_inc(:,:,:))/dt_atmos
    tendency%tracer_inc(:,:,:,:) = (increment_data(2)%tracer_inc(:,:,:,:) - increment_data(1)%tracer_inc(:,:,:,:))/dt_atmos

    write(6,*) 'dt_atmos is ', dt_atmos, mygrid
    write(6,*) increment_data(2)%ua_inc(1,1,1) 
    write(6,*) increment_data(2)%ua_inc(11,1,1) 
    write(6,*) tendency%ua_inc(1,1,1)*dt_atmos/increment_data(2)%ua_inc(1,1,1)
    write(6,*) tendency%ua_inc(11,1,1)*dt_atmos/increment_data(2)%ua_inc(11,1,1)
    write(6,*) tendency%ua_inc(20,11,100)*dt_atmos/increment_data(2)%ua_inc(20,11,100)
    write(6,*) tendency%ua_inc(81,91,15)*dt_atmos/increment_data(2)%ua_inc(81,91,15)
  end subroutine read_netcdf_incs
  
  
!----------------------------------------------------------------------------------------
end module module_get_cubed_sphere_inc
