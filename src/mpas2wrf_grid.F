program convert_mpas

    use scan_input
    use mpas_mesh
    use target_mesh
    use remapper
    use file_output
    use field_list
    use timer
    use map_utils
    use gridinfo_module
    use read_geogrid

    implicit none

    ! Timers
    type (timer_type) :: total_timer, &
                         read_timer, &
                         remap_timer, &
                         write_timer

    integer :: stat, n, i
    character (len=1024) :: mesh_filename, data_filename, namelist_filename
    character (len=1024) :: output_filename, input_filename
!
!
    logical :: path_flag
    type (mpas_mesh_type) :: source_mesh
    type (target_mesh_type) :: destination_mesh_m
    type (target_mesh_type) :: destination_mesh_u_cgrid
    type (target_mesh_type) :: destination_mesh_v_cgrid
    type (target_mesh_type) :: destination_mesh_uv_bgrid
    type (input_handle_type) :: handle
    type (input_field_type) :: field
    type (remap_info_type) :: remap_info_m
    type (remap_info_type) :: remap_info_u_cgrid
    type (remap_info_type) :: remap_info_v_cgrid
    type (remap_info_type) :: remap_info_uv_bgrid
    type (output_handle_type) :: output_handle
    type (target_field_type) :: target_field, mapfac_m, f
    type (field_list_type) :: include_field_list, exclude_field_list
    type (input_global_attributes) :: globatts
    type (input_map_attributes) :: mapatts

    integer :: iRec, iatt_rec
    integer :: nRecordsIn, nRecordsOut
    integer :: iFile
    integer :: fileArgStart, nArgs


    !---- PRISCILLA's READ GEOGRID FILE  mods ----------
    real, dimension(:,:), target, allocatable :: lat2dmf(:,:)
    real, dimension(:,:), target, allocatable :: lon2dmf(:,:)
    real, dimension(:,:), target, allocatable :: lat2duf_cgrid(:,:)
    real, dimension(:,:), target, allocatable :: lon2duf_cgrid(:,:)
    real, dimension(:,:), target, allocatable :: lat2dvf_cgrid(:,:)
    real, dimension(:,:), target, allocatable :: lon2dvf_cgrid(:,:)
    integer :: gf_iproj_type
    real :: gf_known_lat, gf_known_lon, gf_stand_lon, gf_truelat1, gf_truelat2, &
            gf_known_x, gf_known_y, gf_dxkm, gf_dykm
    !---- PAM's GEOGRID NAMELIST mods ----------
    type (proj_info) :: proj
    integer :: i_ll, j_ll, num_lats, num_lons
    real, dimension(:,:), allocatable :: lat2dm 
    real, dimension(:,:), allocatable :: lon2dm 
    real, dimension(:,:), allocatable :: lat2du_cgrid 
    real, dimension(:,:), allocatable :: lon2du_cgrid 
    real, dimension(:,:), allocatable :: lat2dv_cgrid 
    real, dimension(:,:), allocatable :: lon2dv_cgrid 
    real, dimension(:,:), allocatable :: lat2dm_deg
    real, dimension(:,:), allocatable :: lat2duv_bgrid
    real, dimension(:,:), allocatable :: lon2duv_bgrid
    real, dimension(:,:), allocatable :: lon2dm_deg 
    real :: cen_lat, cen_lon


    call timer_start(total_timer)

    if (command_argument_count() < 2) then
        write(0,*) ' '
        write(0,*) 'Usage: convert_mpas mesh-file [data-files]'
        write(0,*) ' '
        write(0,*) 'If only one file argument is given, both the MPAS mesh information and'
        write(0,*) 'the fields will be read from the specified file.'
        write(0,*) 'If two or more file arguments are given, the MPAS mesh information will'
        write(0,*) 'be read from the first file and fields to be remapped will be read from'
        write(0,*) 'the subsequent files.'
        write(0,*) 'All time records from input files will be processed and appended to'
        write(0,*) 'the output file.'
        stop 1
    end if

    nArgs = command_argument_count()
    
    call get_command_argument(2, mesh_filename)
    if (nArgs == 2) then
        fileArgStart = 2
    else 
        fileArgStart = 3
    end if


    ! Have the gridinfo module retrieve description of the grid setup
    call get_command_argument(1, namelist_filename)
    call get_grid_params(namelist_filename)


    if (use_geogrid_file == 0) then
       num_lats=jydim-1
       num_lons=ixdim-1
       call push_source_projection(iproj_type, stand_lon, truelat1, truelat2, &
                       dxkm, dykm, dykm, dxkm, known_x, known_y, &
                       known_lat, known_lon, EARTH_RADIUS_M)
       allocate(lat2dm(ixdim-1,jydim-1))
       allocate(lon2dm(ixdim-1,jydim-1))
       allocate(lat2du_cgrid(ixdim,jydim-1))
       allocate(lon2du_cgrid(ixdim,jydim-1))
       allocate(lat2dv_cgrid(ixdim-1,jydim))
       allocate(lon2dv_cgrid(ixdim-1,jydim))
       allocate(lat2dm_deg(ixdim-1,jydim-1))
       allocate(lon2dm_deg(ixdim-1,jydim-1))
       do i_ll=1,num_lons
          do j_ll=1,num_lats
              call ij_to_latlon(proj_stack,real(i_ll),real(j_ll),lat2dm(i_ll, j_ll),lon2dm(i_ll, j_ll))
              !
              ! Calling target_mesh_setup later changes lat2dm from degrees to radians. 
              ! So need to store lat-lons in degrees for calculation of map factor later. 
              !
              lat2dm_deg(i_ll,j_ll)=lat2dm(i_ll,j_ll)
              lon2dm_deg(i_ll,j_ll)=lon2dm(i_ll,j_ll)
          end do
       end do
       do i_ll=1,num_lons+1
          do j_ll=1,num_lats
              call ij_to_latlon(proj_stack,real(i_ll)-0.5,real(j_ll),lat2du_cgrid(i_ll, j_ll),lon2du_cgrid(i_ll, j_ll))
          end do
       end do
       do i_ll=1,num_lons
          do j_ll=1,num_lats+1
              call ij_to_latlon(proj_stack,real(i_ll),real(j_ll)-0.5,lat2dv_cgrid(i_ll, j_ll),lon2dv_cgrid(i_ll, j_ll))
          end do
       end do
       !
       ! Get map attributes
       !
       call ij_to_latlon(proj_stack,real(ixdim)/2.0,real(jydim)/2.0,cen_lat,cen_lon)
       call create_map_attributes(stand_lon,jydim,ixdim,dxkm,cen_lat,cen_lon, &
                                   iproj_type,truelat1,truelat2,mapatts) 
    elseif (use_geogrid_file == 1) then
       write(0,*) 'GEOGRID FILE is '''//trim(geogrid_filename)//''''
       if (read_geogrid_getlatlons(geogrid_filename, "XLAT_M", "XLONG_M",lat2dmf, lon2dmf, ixdim, jydim) /= 0) then
           write(0,*) 'Error: Problems setting up target mesh'
           stop 2
       end if
       allocate(lat2dm(ixdim,jydim))
       allocate(lon2dm(ixdim,jydim))
       if (read_geogrid_getlatlons(geogrid_filename, "XLAT_U", "XLONG_U",lat2duf_cgrid, lon2duf_cgrid, ixdim, jydim) /= 0) then
           write(0,*) 'Error: Problems setting up target mesh'
           stop 2
       end if
       allocate(lat2du_cgrid(ixdim,jydim))
       allocate(lon2du_cgrid(ixdim,jydim))
       if (read_geogrid_getlatlons(geogrid_filename, "XLAT_V", "XLONG_V",lat2dvf_cgrid, lon2dvf_cgrid, ixdim, jydim) /= 0) then
           write(0,*) 'Error: Problems setting up target mesh'
           stop 2
       end if
       allocate(lat2dv_cgrid(ixdim,jydim))
       allocate(lon2dv_cgrid(ixdim,jydim))
       lat2dm=lat2dmf
       lon2dm=lon2dmf
       lat2du_cgrid=lat2duf_cgrid
       lon2du_cgrid=lon2duf_cgrid
       lat2dv_cgrid=lat2dvf_cgrid
       lon2dv_cgrid=lon2dvf_cgrid
       stat = scan_input_get_map_attributes(geogrid_filename, mapatts)
       do iatt_rec=1,mapatts%num_atts
          if (mapatts%attnames(iatt_rec) == 'MAP_PROJ') then
            gf_iproj_type=mapatts%int_atts(iatt_rec)
          end if
          if (mapatts%attnames(iatt_rec) == 'STAND_LON') then
            gf_stand_lon=mapatts%real_atts(iatt_rec)
          end if
          if (mapatts%attnames(iatt_rec) == 'DX') then
            gf_dxkm=mapatts%real_atts(iatt_rec)
          end if
          if (mapatts%attnames(iatt_rec) == 'DY') then
            gf_dykm=mapatts%real_atts(iatt_rec)
          end if
          if (mapatts%attnames(iatt_rec) == 'CEN_LAT') then
            gf_known_lat=mapatts%real_atts(iatt_rec)
          end if
          if (mapatts%attnames(iatt_rec) == 'CEN_LON') then
            gf_known_lon=mapatts%real_atts(iatt_rec)
          end if
          if (mapatts%attnames(iatt_rec) == 'TRUELAT1') then
            gf_truelat1=mapatts%real_atts(iatt_rec)
          end if
          if (mapatts%attnames(iatt_rec) == 'TRUELAT2') then
            gf_truelat2=mapatts%real_atts(iatt_rec)
          end if
          if (mapatts%attnames(iatt_rec) == 'WEST-EAST_GRID_DIMENSION') then
            gf_known_x=mapatts%int_atts(iatt_rec)/2
            ixdim=mapatts%int_atts(iatt_rec)
            num_lons=ixdim-1
          end if
          if (mapatts%attnames(iatt_rec) == 'SOUTH-NORTH_GRID_DIMENSION') then
            gf_known_y=mapatts%int_atts(iatt_rec)/2
            jydim=mapatts%int_atts(iatt_rec)
            num_lats=jydim-1
          end if
       end do
       call push_source_projection(gf_iproj_type, gf_stand_lon, gf_truelat1, gf_truelat2, &
                       gf_dxkm, gf_dykm, gf_dykm, gf_dxkm, gf_known_x, gf_known_y, &
                       gf_known_lat, gf_known_lon, EARTH_RADIUS_M)
    end if
       allocate(lat2duv_bgrid(ixdim,jydim))
       allocate(lon2duv_bgrid(ixdim,jydim))
       do i_ll=1,num_lons+1
          do j_ll=1,num_lats+1
              call ij_to_latlon(proj_stack,real(i_ll)-0.5,real(j_ll)-0.5,lat2duv_bgrid(i_ll, j_ll),lon2duv_bgrid(i_ll, j_ll))
          end do
       end do
    write(0,*) 'Reading MPAS mesh information from file '''//trim(mesh_filename)//''''

    !
    ! Generate the target grid
    !
    if (target_mesh_setup(destination_mesh_m, lat2dm, lon2dm,'m') /= 0) then
        write(0,*) 'Error: Problems setting up target mesh'
        stop 2
    end if
    if (target_mesh_setup(destination_mesh_u_cgrid, lat2du_cgrid, lon2du_cgrid,'u') /= 0) then
        write(0,*) 'Error: Problems setting up target mesh'
        stop 2
    end if
    if (target_mesh_setup(destination_mesh_v_cgrid, lat2dv_cgrid, lon2dv_cgrid,'v') /= 0) then
        write(0,*) 'Error: Problems setting up target mesh'
        stop 2
    end if
    if (target_mesh_setup(destination_mesh_uv_bgrid, lat2duv_bgrid, lon2duv_bgrid,'b') /= 0) then
        write(0,*) 'Error: Problems setting up target mesh'
        stop 2
    end if

    !
    ! Get information defining the MPAS mesh
    !
    if (mpas_mesh_setup(mesh_filename, source_mesh) /= 0) then
        write(0,*) 'Error: Problems setting up MPAS mesh from file '//trim(mesh_filename)
        stat = target_mesh_free(destination_mesh_m)
        stat = target_mesh_free(destination_mesh_u_cgrid)
        stat = target_mesh_free(destination_mesh_v_cgrid)
        stat = target_mesh_free(destination_mesh_uv_bgrid)
        stop 3
    end if

    !
    ! Compute weights for mapping from MPAS mesh to target grid
    !
    write(0,*) ' '
    write(0,*) 'Computing remapping weights'
    call timer_start(remap_timer)
    if (remap_info_setup(source_mesh, destination_mesh_m, remap_info_m) /= 0) then
        write(0,*) 'Error: Problems setting up remapping'
        stat = mpas_mesh_free(source_mesh)
        stat = target_mesh_free(destination_mesh_m)
        stop 4
    end if
    if (remap_info_setup(source_mesh, destination_mesh_u_cgrid, remap_info_u_cgrid) /= 0) then
        write(0,*) 'Error: Problems setting up remapping'
        stat = mpas_mesh_free(source_mesh)
        stat = target_mesh_free(destination_mesh_u_cgrid)
        stop 4
    end if
    if (remap_info_setup(source_mesh, destination_mesh_v_cgrid, remap_info_v_cgrid) /= 0) then
        write(0,*) 'Error: Problems setting up remapping'
        stat = mpas_mesh_free(source_mesh)
        stat = target_mesh_free(destination_mesh_v_cgrid)
        stop 4
    end if
    if (remap_info_setup(source_mesh, destination_mesh_uv_bgrid, remap_info_uv_bgrid) /= 0) then
        write(0,*) 'Error: Problems setting up remapping'
        stat = mpas_mesh_free(source_mesh)
        stat = target_mesh_free(destination_mesh_uv_bgrid)
        stop 4
    end if
    call timer_stop(remap_timer)
    write(0,'(a,f10.6,a)') '    Time to compute remap weights: ', timer_time(remap_timer), ' s'


    !
    ! Loop over input data files
    !
    do iFile=fileArgStart,nArgs

        call get_command_argument(iFile, data_filename)
        n = len_trim(data_filename)
        do i=n,1,-1
           if (data_filename(i:i) == '/') then
              write(input_filename,'(a)') data_filename(i+1:n)
              path_flag=.TRUE.
              exit
! JGP mod
           else 
	      path_flag=.FALSE.	
! JGP mod
           end if
        end do
        if (.NOT. path_flag) then
           write (input_filename,'(a)') data_filename
        end if
        write(0,*) 'Remapping MPAS fields from file '''//trim(input_filename)//''''

        !
        ! Open input data file
        !
        if (scan_input_open(data_filename, handle, nRecords=nRecordsIn) /= 0) then
            write(0,*) 'Error: Problems opening input file '//trim(data_filename)
            write(0,*) '       This could result from an input file with no unlimited dimension.'
            stat = file_output_close(output_handle)
            stat = scan_input_close(handle)
            stat = mpas_mesh_free(source_mesh)
            stat = target_mesh_free(destination_mesh_m)
            stat = target_mesh_free(destination_mesh_u_cgrid)
            stat = target_mesh_free(destination_mesh_v_cgrid)
            stat = target_mesh_free(destination_mesh_uv_bgrid)
            stat = remap_info_free(remap_info_m)
            stat = remap_info_free(remap_info_u_cgrid)
            stat = remap_info_free(remap_info_v_cgrid)
            stat = remap_info_free(remap_info_uv_bgrid)
            stop 6
        end if
   

        if (frames_per_outfile == 1) then
          !
          ! Open output file
          !

          output_filename=trim(output_prefix) // trim(input_filename)
          write(0,*) 'output_filename = '//trim(output_filename)
          if (file_output_open(output_filename, output_handle, mode=FILE_MODE_APPEND, nRecords=nRecordsOut) /= 0) then
              write(0,*) 'Error: Problems opening output file'
              stat = mpas_mesh_free(source_mesh)
              stat = target_mesh_free(destination_mesh_m)
              stat = target_mesh_free(destination_mesh_u_cgrid)
              stat = target_mesh_free(destination_mesh_v_cgrid)
              stat = target_mesh_free(destination_mesh_uv_bgrid)
              stat = remap_info_free(remap_info_m)
              stat = remap_info_free(remap_info_u_cgrid)
              stat = remap_info_free(remap_info_v_cgrid)
              stat = remap_info_free(remap_info_uv_bgrid)
              stop 5
          end if

          if (nRecordsOut /= 0) then
              write(0,*) 'Existing output file has ', nRecordsOut, ' records'
          else
              write(0,*) 'Created a new output file'
          end if

          !
          ! Get list of fields to include or exclude from input file
          !
          stat = field_list_init(include_field_list, exclude_field_list)


        end if

 
        write(0,*) 'Input file has ', nRecordsIn, ' records'
    
        ! generally, we should make sure dimensions match in existing output files
        ! and in subsequent MPAS input files
    
        !
        ! Scan through input file, determine which fields will be remapped,
        ! and define those fields in the output file; this only needs to be done
        ! if there are no existing records in the output file (i.e., the output
        ! file is a new file)
        !
        write(0,*) ' '
        if (nRecordsOut == 0) then
            write(0,*) 'Defining fields in output file'

            ! Define 'lat' and 'lon' fields for target mesh
            stat = remap_get_target_latitudes(remap_info_m, target_field,'M')
            stat = file_output_register_field(output_handle, target_field)
            stat = free_target_field(target_field)

            stat = remap_get_target_longitudes(remap_info_m, target_field,'M')
            stat = file_output_register_field(output_handle, target_field)
            stat = free_target_field(target_field)

            ! Define 'lat_u' and 'lon_u' fields for target mesh
            stat = remap_get_target_latitudes(remap_info_u_cgrid, target_field, 'U')
            stat = file_output_register_field(output_handle, target_field)
            stat = free_target_field(target_field)

            stat = remap_get_target_longitudes(remap_info_u_cgrid, target_field, 'U')
            stat = file_output_register_field(output_handle, target_field)
            stat = free_target_field(target_field)

            ! Define 'lat_v' and 'lon_v' fields for target mesh
            stat = remap_get_target_latitudes(remap_info_v_cgrid, target_field, 'V')
            stat = file_output_register_field(output_handle, target_field)
            stat = free_target_field(target_field)

            stat = remap_get_target_longitudes(remap_info_v_cgrid, target_field, 'V')
            stat = file_output_register_field(output_handle, target_field)
            stat = free_target_field(target_field)

            ! Define 'lat_uv_Bgrid' and 'lon_uv_Bgrid' fields for target mesh
            stat = remap_get_target_latitudes(remap_info_uv_bgrid, target_field, 'B')
            stat = file_output_register_field(output_handle, target_field)
            stat = free_target_field(target_field)

            stat = remap_get_target_longitudes(remap_info_uv_bgrid, target_field, 'B')
            stat = file_output_register_field(output_handle, target_field)
            stat = free_target_field(target_field)

            ! Define 'mapfac_m' field for target mesh
            if (use_geogrid_file == 0) then
               call get_map_factor(lat2dm_deg, lon2dm_deg, target_field, ixdim-1, jydim-1 )
            elseif (use_geogrid_file == 1) then 
               stat = read_geogrid_getmapfac(geogrid_filename, target_field)
            end if
            stat = file_output_register_field(output_handle, target_field)
            stat = free_target_field(target_field)

            ! Define 'F' field for target mesh
            if (use_geogrid_file == 0) then
               call get_coriolis_parameters(lat2dm, target_field, ixdim-1, jydim-1 )
            elseif (use_geogrid_file == 1) then
               stat = read_geogrid_getcoriol(geogrid_filename, target_field) 
            end if
            stat = file_output_register_field(output_handle, target_field)
            stat = free_target_field(target_field)

            ! Define 'xtime' field for target mesh
            stat = remap_get_xtime(handle, target_field)
            write(0,*) 'xtime = ', target_field % array1c
            stat = file_output_register_field(output_handle, target_field)
            stat = free_target_field(target_field)

            do while (scan_input_next_field(handle, field) == 0) 
                if (can_remap_field(field) .and. &
                    should_remap_field(field, include_field_list, exclude_field_list)) then
                    if (field % name == 'uReconstructZonal') then
                       if (uv_gridtype == 'A' .OR. uv_gridtype == 'a') then
                          stat = remap_field_dryrun(remap_info_m, field, target_field)
                          print*,'variable = ',field % name
                          print*,'using Arakawa A grid type' 
                       elseif (uv_gridtype == 'B' .OR. uv_gridtype == 'b') then
                          stat = remap_field_dryrun(remap_info_uv_bgrid, field, target_field)
                          print*,'variable = ',field % name
                          print*,'using Arakawa B grid type' 
                       elseif (uv_gridtype == 'C' .OR. uv_gridtype == 'c') then
                          stat = remap_field_dryrun(remap_info_u_cgrid, field, target_field)
                          print*,'variable = ',field % name
                          print*,'using Arakawa C grid type'
                       end if
                    elseif (field % name == 'uReconstructMeridional') then
                       if (uv_gridtype == 'A' .OR. uv_gridtype == 'a') then
                          stat = remap_field_dryrun(remap_info_m, field, target_field)
                          print*,'variable = ',field % name
                          print*,'using Arakawa A grid type' 
                       elseif (uv_gridtype == 'B' .OR. uv_gridtype == 'b') then
                          stat = remap_field_dryrun(remap_info_uv_bgrid, field, target_field)
                          print*,'variable = ',field % name
                          print*,'using Arakawa B grid type' 
                       elseif (uv_gridtype == 'C' .OR. uv_gridtype == 'c') then
                          stat = remap_field_dryrun(remap_info_v_cgrid, field, target_field)
                          print*,'variable = ',field % name
                          print*,'using Arakawa C grid type'
                       end if
                    else 
                       stat = remap_field_dryrun(remap_info_m, field, target_field)
                       print*,'using m coords'
                       print*,'variable = ',field % name
                    end if 
                    stat = file_output_register_field(output_handle, target_field)
                    if (stat /= 0) then
                        stat = free_target_field(target_field)
                        stat = scan_input_free_field(field)
                        stat = scan_input_close(handle)
                        stat = file_output_close(output_handle)
                        stat = mpas_mesh_free(source_mesh)
                        stat = target_mesh_free(destination_mesh_m)
                        stat = target_mesh_free(destination_mesh_u_cgrid)
                        stat = target_mesh_free(destination_mesh_v_cgrid)
                        stat = target_mesh_free(destination_mesh_uv_bgrid)
                        stat = remap_info_free(remap_info_m)
                        stat = remap_info_free(remap_info_u_cgrid)
                        stat = remap_info_free(remap_info_v_cgrid)
                        stat = remap_info_free(remap_info_uv_bgrid)
                        stat = field_list_finalize(include_field_list, exclude_field_list)
                        stop 7
                    end if
                    stat = free_target_field(target_field)
                end if
                stat = scan_input_free_field(field)
            end do
        end if


        ! 
        ! Write 'lat' and 'lon' fields for target mesh
        ! And write 'mapfac_m' and 'F' fields for target mesh
        ! And write 'xtime' to output file 
        ! 
        if (nRecordsOut == 0) then
            stat = remap_get_target_latitudes(remap_info_m, target_field, 'M')
            stat = file_output_write_field(output_handle, target_field, frame=0)
            stat = free_target_field(target_field)

            stat = remap_get_target_longitudes(remap_info_m, target_field, 'M')
            stat = file_output_write_field(output_handle, target_field, frame=0)
            stat = free_target_field(target_field)
        
            stat = remap_get_target_latitudes(remap_info_u_cgrid, target_field, 'U')
            stat = file_output_write_field(output_handle, target_field, frame=0)
            stat = free_target_field(target_field)

            stat = remap_get_target_longitudes(remap_info_u_cgrid, target_field, 'U')
            stat = file_output_write_field(output_handle, target_field, frame=0)
            stat = free_target_field(target_field)
        
            stat = remap_get_target_latitudes(remap_info_v_cgrid, target_field, 'V')
            stat = file_output_write_field(output_handle, target_field, frame=0)
            stat = free_target_field(target_field)

            stat = remap_get_target_longitudes(remap_info_v_cgrid, target_field, 'V')
            stat = file_output_write_field(output_handle, target_field, frame=0)
            stat = free_target_field(target_field)
        
            stat = remap_get_target_latitudes(remap_info_uv_bgrid, target_field, 'B')
            stat = file_output_write_field(output_handle, target_field, frame=0)
            stat = free_target_field(target_field)

            stat = remap_get_target_longitudes(remap_info_uv_bgrid, target_field, 'B')
            stat = file_output_write_field(output_handle, target_field, frame=0)
            stat = free_target_field(target_field)

            if (use_geogrid_file == 0) then
               call get_map_factor(lat2dm_deg, lon2dm_deg, target_field, ixdim-1, jydim-1 )
            elseif (use_geogrid_file == 1) then
               stat = read_geogrid_getmapfac(geogrid_filename, target_field)
            end if
            stat = file_output_write_field(output_handle, target_field, frame=0)
            stat = free_target_field(target_field)
            call get_coriolis_parameters(lat2dm, target_field, ixdim-1, jydim-1 )

            if (use_geogrid_file == 0) then
               call get_coriolis_parameters(lat2dm, target_field, ixdim-1, jydim-1 )
            elseif (use_geogrid_file == 1) then
               stat = read_geogrid_getcoriol(geogrid_filename, target_field) 
            end if
            stat = file_output_write_field(output_handle, target_field, frame=0)
            stat = free_target_field(target_field)
        end if

    
        !
        ! Loop over all times in the input file
        !
        do iRec=1,nRecordsIn
            stat = scan_input_rewind(handle)
    
            !
            ! Scan through list of fields in the input file, remapping fields and writing
            ! them to the output file
            !
            do while (scan_input_next_field(handle, field) == 0) 
                if (can_remap_field(field) .and. &
                    should_remap_field(field, include_field_list, exclude_field_list)) then
                    write(0,*) 'Remapping field '//trim(field % name)//', frame ', irec
    
                    call timer_start(read_timer)
                    stat = scan_input_read_field(field, frame=iRec)
                    call timer_stop(read_timer)
                    write(0,'(a,f10.6,a)') '    read: ', timer_time(read_timer), ' s'
    
                    call timer_start(remap_timer)
                    if (field % name == 'uReconstructZonal') then
                       if (uv_gridtype == 'A' .OR. uv_gridtype == 'a') then
                          stat = remap_field(remap_info_m, field, target_field)
                          print*,'using Arakawa A grid type' 
                       elseif (uv_gridtype == 'B' .OR. uv_gridtype == 'b') then
                          stat = remap_field(remap_info_uv_bgrid, field, target_field)
                          print*,'using Arakawa B grid type' 
                       elseif (uv_gridtype == 'C' .OR. uv_gridtype == 'c') then
                          stat = remap_field(remap_info_u_cgrid, field, target_field)
                          print*,'using Arakawa C grid type'
                       end if
                    elseif (field % name == 'uReconstructMeridional') then
                       stat = remap_field(remap_info_v_cgrid, field, target_field)
                       print*,'using v coords'
                       if (uv_gridtype == 'B' .OR. uv_gridtype == 'a') then
                          stat = remap_field(remap_info_m, field, target_field)
                          print*,'using Arakawa A grid type' 
                       elseif (uv_gridtype == 'B' .OR. uv_gridtype == 'b') then
                          stat = remap_field(remap_info_uv_bgrid, field, target_field)
                          print*,'using Arakawa B grid type' 
                       elseif (uv_gridtype == 'C' .OR. uv_gridtype == 'c') then
                          stat = remap_field(remap_info_v_cgrid, field, target_field)
                          print*,'using Arakawa C grid type'
                       end if
                    else 
                       stat = remap_field(remap_info_m, field, target_field)
                       print*,'using m coords'
                    end if 
                    call timer_stop(remap_timer)
                    write(0,'(a,f10.6,a)') '    remap: ', timer_time(remap_timer), ' s'
    
                    call timer_start(write_timer)
                    stat = file_output_write_field(output_handle, target_field, frame=(nRecordsOut+iRec))
                    call timer_stop(write_timer)
                    write(0,'(a,f10.6,a)') '    write: ', timer_time(write_timer), ' s'
    
                    stat = free_target_field(target_field)
                end if
                stat = scan_input_free_field(field)
            end do
        stat = remap_get_xtime(handle, target_field)
        stat = file_output_write_field(output_handle, target_field,frame=(nRecordsOut+iRec))
        stat = free_target_field(target_field)
        end do
    
        nRecordsOut = nRecordsOut + nRecordsIn

        write(0,*) 'Writing out global attributes from MPAS'
        stat = scan_input_get_global_attributes(handle, globatts)
        if (use_geogrid_file == 1) then
           stat = scan_input_get_map_attributes(geogrid_filename, mapatts)
        end if

        stat = file_output_global_attributes(output_handle, globatts)
        stat = file_output_map_attributes(output_handle, mapatts)

        stat = scan_input_close(handle)
    end do


    !
    ! Cleanup
    !
    stat = file_output_close(output_handle)

    stat = mpas_mesh_free(source_mesh)
    stat = target_mesh_free(destination_mesh_m)
    stat = target_mesh_free(destination_mesh_u_cgrid)
    stat = target_mesh_free(destination_mesh_v_cgrid)
    stat = target_mesh_free(destination_mesh_uv_bgrid)
    stat = remap_info_free(remap_info_m)
    stat = remap_info_free(remap_info_u_cgrid)
    stat = remap_info_free(remap_info_v_cgrid)
    stat = remap_info_free(remap_info_uv_bgrid)
    stat = field_list_finalize(include_field_list, exclude_field_list)

    call timer_stop(total_timer)

    write(0,*) ' '
    write(0,'(a,f10.6)') 'Total runtime: ', timer_time(total_timer)
    write(0,*) ' '

    stop

end program convert_mpas
