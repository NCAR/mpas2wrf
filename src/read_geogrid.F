module read_geogrid

   use netcdf
   use remapper

   contains

   integer function read_geogrid_getlatlons(filename, lat_name, lon_name, lat2d_geo, lon2d_geo, numLatLons, numLatLats) result(stat)
   
   implicit none
   
   
   ! This is the name of the data file we will read. 
   character (len = *), intent(in) :: filename 
   character (len = *), intent(in) :: lat_name 
   character (len = *), intent(in) :: lon_name 

   ! We are reading 3D data
   real, dimension(:,:,:), allocatable :: lat2d_time
   real, dimension(:,:,:), allocatable :: lon2d_time
   
   real, dimension(:,:), allocatable, intent(out) :: lat2d_geo
   real, dimension(:,:), allocatable, intent(out) :: lon2d_geo

   ! This will be the netCDF ID for the file and data variable.
   integer :: ncid, lonid, latid, lat_type, lat_ndims
   integer :: numLatLons, numLatLats, numLatTimes
   integer :: numLonLons, numLonLats, numLonTimes
   integer, dimension(nf90_max_var_dims) :: LondimIDs, LatdimIDs

   ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
   ! the file.
   stat = nf90_open(filename, NF90_NOWRITE, ncid) 
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if


   ! Get the varid of the XLAT_M variable, based on its name.
   stat = nf90_inq_varid(ncid, lat_name, latid)
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if 
   ! Get the varid of the XLONG_M variable, based on its name.
   stat = nf90_inq_varid(ncid, lon_name, lonid)
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if
 
   ! Get the varid of the XLAT_M variable, based on its name.
   stat = nf90_inquire_variable(ncid, latid, xtype=lat_type, ndims=lat_ndims)
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if 
   ! --- GET DIMENSION IDS and size for LATS ---
   stat = nf90_inquire_variable(ncid, latId, dimids = latdimIDs)
   if(stat /= nf90_NoErr) then
      stat = 1
       return
   end if
   stat = nf90_inquire_dimension(ncid, latdimIDs(1), len = numLatLons)
   if(stat /= nf90_NoErr) then
       stat = 1
       return
   end if 
   stat = nf90_inquire_dimension(ncid, latdimIDs(2), len = numLatLats)
   if(stat /= nf90_NoErr) then
      stat = 1
       return
   end if 
   stat = nf90_inquire_dimension(ncid, latdimIDs(3), len = numLatTimes)
   if(stat /= nf90_NoErr) then
       stat = 1
       return
   end if 
   allocate(lat2d_time(numLatLons, numLatLats, numLatTimes))
   allocate(lat2d_geo(numLatLons, numLatLats))

   ! --- FINISH GETTING DIMENSION IDS and size for LATS ---

   ! Read the LATITUDES.
   stat = nf90_get_var(ncid, latid, lat2d_time) 
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if
   ! STORE the LATITUDES in a 2D ARRAY i.e. remove the time dimension.
   lat2d_geo=lat2d_time(:,:,1)

   ! --- GET DIMENSION IDS and size for LATS ---
   stat = nf90_inquire_variable(ncid, lonId, dimids = londimIDs)
   if(stat /= nf90_NoErr) then
      stat = 1
       return
   end if
   stat = nf90_inquire_dimension(ncid, londimIDs(1), len = numLonLons)
   if(stat /= nf90_NoErr) then
       stat = 1
       return
   end if 
   stat = nf90_inquire_dimension(ncid, londimIDs(2), len = numLonLats)
   if(stat /= nf90_NoErr) then
      stat = 1
       return
   end if 
   stat = nf90_inquire_dimension(ncid, londimIDs(3), len = numLonTimes)
   if(stat /= nf90_NoErr) then
       stat = 1
       return
   end if 
   allocate(lon2d_time(numLonLons, numLonLats, numLonTimes))
   allocate(lon2d_geo(numLonLons, numLonLats))

   ! --- FINISH GETTING DIMENSION IDS and size for LATS ---

   ! Read the LONGITUDES.
   stat = nf90_get_var(ncid, lonid, lon2d_time)
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if
   ! STORE the LONGITUDES in a 2D ARRAY i.e. remove the time dimension.
   lon2d_geo=lon2d_time(:,:,1)

  ! Close the file, freeing all resources.
   stat =  nf90_close(ncid) 
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if

  end function read_geogrid_getlatlons


 integer function read_geogrid_getmapfac(filename, mapfac_field) result(stat)
   
   implicit none
  
   type (target_field_type), intent(out) :: mapfac_field
   ! This is the name of the data file we will read. 
   character (len = *), intent(in) :: filename 

   ! We are reading 3D data
   real, dimension(:,:,:), allocatable :: mapfac_time 
   
   real, dimension(:,:), allocatable :: mapfac_geo 

   ! This will be the netCDF ID for the file and data variable.
   integer :: ncid, mapfacid, mapfac_type, mapfac_ndims
   integer :: nummapfacLons, nummapfacLats, nummapfacTimes
   integer, dimension(nf90_max_var_dims) :: mapfacdimIDs

   ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
   ! the file.
   stat = nf90_open(filename, NF90_NOWRITE, ncid) 
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if

   ! Get the varid of the MAPFAC_M variable, based on its name.
   stat = nf90_inq_varid(ncid, "MAPFAC_M", mapfacid)
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if 
 
   ! Get the varid of the XLAT_M variable, based on its name.
   stat = nf90_inquire_variable(ncid, mapfacid, xtype=mapfac_type, ndims=mapfac_ndims)
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if 
        if (mapfac_type == NF90_FLOAT) then
            mapfac_type = FIELD_TYPE_REAL
        else if (mapfac_type == NF90_DOUBLE) then
            mapfac_type = FIELD_TYPE_DOUBLE
        else if (mapfac_type == NF90_INT) then
            mapfac_type = FIELD_TYPE_INTEGER
        else if (mapfac_type == NF90_CHAR) then
            mapfac_type = FIELD_TYPE_CHARACTER
        else
            mapfac_type = FIELD_TYPE_UNSUPPORTED
        end if
   ! --- GET DIMENSION IDS and size for LATS ---
   stat = nf90_inquire_variable(ncid, mapfacid, dimids = mapfacdimIDs)
   if(stat /= nf90_NoErr) then
      stat = 1
       return
   end if
   stat = nf90_inquire_dimension(ncid, mapfacdimIDs(1), len = nummapfacLons)
   if(stat /= nf90_NoErr) then
       stat = 1
       return
   end if 
   stat = nf90_inquire_dimension(ncid, mapfacdimIDs(2), len = nummapfacLats)
   if(stat /= nf90_NoErr) then
      stat = 1
       return
   end if 
   stat = nf90_inquire_dimension(ncid, mapfacdimIDs(3), len = nummapfacTimes)
   if(stat /= nf90_NoErr) then
       stat = 1
       return
   end if 
   allocate(mapfac_time(nummapfacLons, nummapfacLats, nummapfacTimes))
   allocate(mapfac_geo(nummapfacLons, nummapfacLats))
   

   ! --- FINISH GETTING DIMENSION IDS and size for LATS ---

   ! Read the MAP FACTOR. 
   stat = nf90_get_var(ncid, mapfacid, mapfac_time) 
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if
   ! STORE the LATITUDES in a 2D ARRAY i.e. remove the time dimension.
   mapfac_geo=mapfac_time(:,:,1)
   ! STORE the mapfac variable in a field 
   mapfac_field % name = 'MAPFAC_M'
   mapfac_field % ndims = 2 
   mapfac_field % xtype = mapfac_type
   mapfac_field % isTimeDependent = .false.
   allocate(mapfac_field % dimnames(mapfac_field % ndims))
   allocate(mapfac_field % dimlens(mapfac_field % ndims))
   mapfac_field % dimnames(1) = 'nLon_m'
   mapfac_field % dimnames(2) = 'nLat_m'
   mapfac_field % dimlens(1) = nummapfacLons
   mapfac_field % dimlens(2) = nummapfacLats 
   allocate(mapfac_field % array2r(mapfac_field % dimlens(1),mapfac_field % dimlens(2))) 
   mapfac_field % array2r(:,:) = mapfac_geo


  ! Close the file, freeing all resources.
   stat =  nf90_close(ncid) 
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if

  end function read_geogrid_getmapfac
 integer function read_geogrid_getcoriol(filename, cor_field) result(stat)
   
   implicit none
  
   type (target_field_type), intent(out) :: cor_field 
   ! This is the name of the data file we will read. 
   character (len = *), intent(in) :: filename 

   ! We are reading 3D data
   real, dimension(:,:,:), allocatable :: cor_param_time 
   
   real, dimension(:,:), allocatable :: cor_param_geo 

   ! This will be the netCDF ID for the file and data variable.
   integer :: ncid, corid, cor_type, cor_ndims
   integer :: numcorLons, numcorLats, numcorTimes
   integer, dimension(nf90_max_var_dims) :: cordimIDs

   ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
   ! the file.
   stat = nf90_open(filename, NF90_NOWRITE, ncid) 
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if


   ! Get the varid of the F variable, based on its name.
   stat = nf90_inq_varid(ncid, "F", corid)
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if
 
   ! Get the varid of the F variable, based on its name.
   stat = nf90_inquire_variable(ncid, corid, xtype=cor_type, ndims=cor_ndims)
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if 
        if (cor_type == NF90_FLOAT) then
            cor_type = FIELD_TYPE_REAL
        else if (cor_type == NF90_DOUBLE) then
            cor_type = FIELD_TYPE_DOUBLE
        else if (cor_type == NF90_INT) then
            cor_type = FIELD_TYPE_INTEGER
        else if (cor_type == NF90_CHAR) then
            cor_type = FIELD_TYPE_CHARACTER
        else
            cor_type = FIELD_TYPE_UNSUPPORTED
        end if
   ! --- GET DIMENSION IDS and size for LATS ---
   stat = nf90_inquire_variable(ncid, corId, dimids = cordimIDs)
   if(stat /= nf90_NoErr) then
      stat = 1
       return
   end if
   stat = nf90_inquire_dimension(ncid, cordimIDs(1), len = numcorLons)
   if(stat /= nf90_NoErr) then
       stat = 1
       return
   end if 
   stat = nf90_inquire_dimension(ncid, cordimIDs(2), len = numcorLats)
   if(stat /= nf90_NoErr) then
      stat = 1
       return
   end if 
   stat = nf90_inquire_dimension(ncid, cordimIDs(3), len = numcorTimes)
   if(stat /= nf90_NoErr) then
       stat = 1
       return
   end if 
   allocate(cor_param_time(numcorLons, numcorLats, numcorTimes))
   allocate(cor_param_geo(numcorLons, numcorLats))

   ! --- FINISH GETTING DIMENSION IDS and size for LATS ---

   ! Read the LONGITUDES.
   stat = nf90_get_var(ncid, corid, cor_param_time)
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if
   ! STORE the LONGITUDES in a 2D ARRAY i.e. remove the time dimension.
   cor_param_geo=cor_param_time(:,:,1)
   ! STORE the variable F in a field type
   cor_field % name = 'F'
   cor_field % ndims = 2 
   cor_field % xtype = cor_type
   cor_field % isTimeDependent = .false.
   allocate(cor_field % dimnames(cor_field % ndims))
   allocate(cor_field % dimlens(cor_field % ndims))
   cor_field % dimnames(1) = 'nLon_m'
   cor_field % dimnames(2) = 'nLat_m'
   cor_field % dimlens(1) = numcorLons
   cor_field % dimlens(2) = numcorLats 
   allocate(cor_field % array2r(cor_field % dimlens(1),cor_field % dimlens(2))) 
   cor_field % array2r(:,:)=cor_param_geo
  ! Close the file, freeing all resources.
   stat =  nf90_close(ncid) 
   if (stat /= NF90_NOERR) then
       stat = 1
       return
   end if

  end function read_geogrid_getcoriol
end module read_geogrid 

