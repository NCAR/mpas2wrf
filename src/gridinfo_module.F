!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE GRIDINFO_MODULE
!
! This module handles (i.e., acquires, stores, and makes available) all data
!   describing the model domains to be processed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module gridinfo_module

   use constants_module
   use misc_definitions_module
   use scan_input
 
   ! Variables
   integer :: iproj_type, io_form_output, dyn_opt
   real :: known_lat, known_lon, pole_lat, pole_lon, stand_lon, truelat1, truelat2, &
           known_x, known_y, dxkm, dykm, phi, lambda, ref_lat, ref_lon, ref_x, ref_y, &
           dlatdeg, dlondeg
   character (len=1) :: gridtype,uv_gridtype
   character (len=1024) :: output_prefix, geogrid_filename
   integer :: ixdim, jydim, frames_per_outfile, use_geogrid_file 
   contains
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! Name: get_grid_params
   !
   ! Purpose: This subroutine retrieves all parameters regarding the model domains
   !    to be processed by geogrid.exe. This includes map parameters, domain
   !    size and location, and nest information. 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine get_grid_params(namelist_file)
 
      implicit none
 

      character (len=*), intent(in) :: namelist_file
 
      ! Local variables
      integer :: i, j, funit, io_form_geogrid, interval_seconds
      real :: dx, dy 
      character (len=128) :: map_proj
      logical :: is_used 
      integer :: s_we, e_we, s_sn, e_sn
      type (input_map_attributes) :: mapatts
      integer, parameter :: FIELD_TYPE_UNSUPPORTED = -1, &
                            FIELD_TYPE_REAL        =  1, &
                            FIELD_TYPE_DOUBLE      =  2, &
                            FIELD_TYPE_INTEGER     =  3, &
                            FIELD_TYPE_CHARACTER   =  4

      namelist /geogrid/ s_we, e_we, s_sn, e_sn, &
                         map_proj, ref_x, ref_y, ref_lat, ref_lon, &
                         pole_lat, pole_lon, truelat1, truelat2, stand_lon, &
                         dx, dy, geogrid_filename
      namelist /share/ output_prefix, frames_per_outfile, use_geogrid_file, uv_gridtype
 
      ! Set defaults for namelist variables
!      io_form_geogrid = 2
      output_prefix='NEW'
      frames_per_outfile = 1
      use_geogrid_file = 0
      ref_x = NAN
      ref_y = NAN
      ref_lat = NAN
      ref_lon = NAN
      dx = NAN
      dy = NAN
      map_proj = 'Lambert'
      pole_lat = 90.0
      pole_lon = 0.0
      truelat1 = NAN
      truelat2 = NAN
      stand_lon = NAN
      s_we = 1
      e_we = INVALID
      s_sn = 1
      e_sn = INVALID
      geogrid_filename = 'geo_em.d01.nc'
      uv_gridtype = 'b' 
      ! Read parameters from Fortran namelist
      do funit=10,100
         inquire(unit=funit, opened=is_used)
         if (.not. is_used) exit
      end do
      open(funit,file=namelist_file,status='old',form='formatted',err=1000)
      read(funit,share)
      close(funit)
      open(funit,file=namelist_file,status='old',form='formatted',err=1000)
      read(funit,geogrid)
      close(funit)


! BUG: More properly handle debug_level in module_debug

      dxkm = dx
      dykm = dy

      known_lat = ref_lat
      known_lon = ref_lon
      known_x = ref_x
      known_y = ref_y

      ! Before doing anything else, we must have a valid grid type 
         gridtype = 'C'
         dyn_opt = 2

      ! Check for valid io_form_geogrid
      if ( &
#ifdef IO_BINARY
          io_form_geogrid /= BINARY .and. &
#endif
#ifdef IO_NETCDF
          io_form_geogrid /= NETCDF .and. &
#endif
#ifdef IO_GRIB1
          io_form_geogrid /= GRIB1 .and. &
#endif
          .true. ) then
      end if
      io_form_output = io_form_geogrid
  
      ! Convert map_proj to uppercase letters
      do i=1,len(map_proj)
         if (ichar(map_proj(i:i)) >= 97) map_proj(i:i) = char(ichar(map_proj(i:i))-32)
      end do
  
      ! Assign parameters to module variables
      if ((index(map_proj, 'LAMBERT') /= 0) .and. &
          (len_trim(map_proj) == len('LAMBERT'))) then
         iproj_type = PROJ_LC 
  
      else if ((index(map_proj, 'MERCATOR') /= 0) .and. &
               (len_trim(map_proj) == len('MERCATOR'))) then
         iproj_type = PROJ_MERC 
  
      else if ((index(map_proj, 'POLAR') /= 0) .and. &
               (len_trim(map_proj) == len('POLAR'))) then
         iproj_type = PROJ_PS 
  
      else if ((index(map_proj, 'ROTATED_LL') /= 0) .and. &
               (len_trim(map_proj) == len('ROTATED_LL'))) then
         iproj_type = PROJ_ROTLL 
  
      else if ((index(map_proj, 'LAT-LON') /= 0) .and. &
               (len_trim(map_proj) == len('LAT-LON'))) then
         iproj_type = PROJ_CASSINI 
  
      end if

      ! For Cassini / lat-lon projections
      if (iproj_type == PROJ_CASSINI) then

         ! If no dx,dy specified, assume global grid
         if (dx == NAN .and. dy == NAN) then
            dlondeg = 360. / (e_we-s_we)   ! Here, we really do not want e_we-s_we+1
            dlatdeg = 180. / (e_sn-s_sn)   ! Here, we really do not want e_we-s_we+1
            known_x = 1.
            known_y = 1.
            known_lon = stand_lon + dlondeg/2.
            known_lat = -90. + dlatdeg/2.
            dxkm = EARTH_RADIUS_M * PI * 2.0 / (e_we-s_we)
            dykm = EARTH_RADIUS_M * PI       / (e_sn-s_sn)

         ! If dx,dy specified, however, assume regional grid
         else
            dlatdeg = dy
            dlondeg = dx
            dxkm = dlondeg * EARTH_RADIUS_M * PI * 2.0 / 360.0
            dykm = dlatdeg * EARTH_RADIUS_M * PI * 2.0 / 360.0
            if (known_lat == NAN .or. known_lon == NAN) then
               print*,'For lat-lon projection, if dx/dy are specified, '// &
                      'a regional domain is assumed, and a ref_lat,ref_lon must also be specified'
            end if
         end if
      end if

      ! Manually set truelat2 = truelat1 if truelat2 not specified for Lambert
      if (iproj_type == PROJ_LC .and. truelat2 == NAN) then
         if (truelat1 == NAN) then
            print*,"No TRUELAT1 specified for Lambert conformal projection." 
            truelat2 = truelat1
         end if
      end if

  
      ! For C grid, let ixdim and jydim be the number of velocity points in 
      !    each direction; 
         ixdim = e_we - s_we + 1
         jydim = e_sn - s_sn + 1
  
      ! If the user hasn't supplied a known_x and known_y, assume the center of domain 1
         if (known_x == NAN .and. known_y == NAN) then
            known_x = ixdim / 2.
            known_y = jydim / 2.
         else if (known_x == NAN .or. known_y == NAN) then
            print*, 'In namelist file, neither or both of ref_x, ref_y must be specified.'
         end if 

         ! C grid does not support the rotated lat/lon projection
          if (iproj_type == PROJ_ROTLL) then
               print*,'Rotated lat/lon projection is not supported '// &
                      'Valid projecitons are "lambert", "mercator", "polar", and "lat-lon".'
          end if

      return
  
 1000 print*,'Error opening file namelist file'
 
   end subroutine get_grid_params

   subroutine create_map_attributes(stand_lon,e_sn,e_we,dx,cen_lat,cen_lon, &
                                   iproj_type,truelat1,truelat2,mapatts)
    
      IMPLICIT NONE

      ! Declare arguments
      REAL, INTENT(IN)        :: stand_lon
      INTEGER, INTENT(IN)     :: e_sn 
      INTEGER, INTENT(IN)     :: e_we 
      REAL, INTENT(IN)        :: dx
      REAL, INTENT(IN)        :: cen_lat
      REAL, INTENT(IN)        :: cen_lon
      INTEGER, INTENT(IN)     :: iproj_type 
      REAL, INTENT(IN)        :: truelat1
      REAL, INTENT(IN)        :: truelat2
      TYPE (input_map_attributes), INTENT(OUT)      :: mapatts 

      INTEGER, PARAMETER :: FIELD_TYPE_UNSUPPORTED = -1, &
                            FIELD_TYPE_REAL        =  1, &
                            FIELD_TYPE_DOUBLE      =  2, &
                            FIELD_TYPE_INTEGER     =  3, &
                            FIELD_TYPE_CHARACTER   =  4

 
      ! Create global attributes for map parameters
      ! ---- Set attribute names ----
      mapatts % num_atts = 9
      ! ---- Set attribute names ----
      mapatts % attnames(1)='STAND_LON'
      mapatts % attnames(2)='SOUTH-NORTH_GRID_DIMENSION'
      mapatts % attnames(3)='WEST-EAST_GRID_DIMENSION'
      mapatts % attnames(4)='DX'
      mapatts % attnames(5)='CEN_LAT'
      mapatts % attnames(6)='CEN_LON'
      mapatts % attnames(7)='MAP_PROJ'
      mapatts % attnames(8)='TRUELAT1'
      mapatts % attnames(9)='TRUELAT2'
      ! ---- Set attribute values ----
      mapatts % real_atts(1) = stand_lon
      mapatts % int_atts(2) = e_sn
      mapatts % int_atts(3) = e_we
      mapatts % real_atts(4) = dx
      mapatts % real_atts(5) = cen_lat
      mapatts % real_atts(6) = cen_lon
      mapatts % int_atts(7) = iproj_type
      mapatts % real_atts(8) = truelat1
      mapatts % real_atts(9) = truelat2
      ! ---- Set attribute types ----
      mapatts % att_type(1) = FIELD_TYPE_REAL
      mapatts % att_type(2) = FIELD_TYPE_INTEGER
      mapatts % att_type(3) = FIELD_TYPE_INTEGER
      mapatts % att_type(4) = FIELD_TYPE_REAL
      mapatts % att_type(5) = FIELD_TYPE_REAL
      mapatts % att_type(6) = FIELD_TYPE_REAL
      mapatts % att_type(7) = FIELD_TYPE_INTEGER
      mapatts % att_type(8) = FIELD_TYPE_REAL
      mapatts % att_type(9) = FIELD_TYPE_REAL

   end subroutine create_map_attributes
  
end module gridinfo_module
