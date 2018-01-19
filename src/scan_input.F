module scan_input

    use netcdf

    type input_handle_type
        integer :: ncid
        integer :: num_vars = 0
        integer :: current_var = 0
        integer, dimension(:), pointer :: varids => null()
        integer :: unlimited_dimid
    end type input_handle_type

    type input_field_type
        character (len=64) :: name
        logical :: isTimeDependent = .false.
        integer :: varid = -1
        integer :: xtype = -1
        integer :: ndims = -1
        character (len=64), dimension(:), pointer :: dimnames
        integer, dimension(:), pointer :: dimlens
        integer, dimension(:), pointer :: dimids
        type (input_handle_type), pointer :: file_handle

        !  Members to store field data
        real :: array0r
        real, dimension(:), pointer :: array1r => null()
        real, dimension(:,:), pointer :: array2r => null()
        real, dimension(:,:,:), pointer :: array3r => null()
        real, dimension(:,:,:,:), pointer :: array4r => null()
        double precision :: array0d
        double precision, dimension(:), pointer :: array1d => null()
        double precision, dimension(:,:), pointer :: array2d => null()
        double precision, dimension(:,:,:), pointer :: array3d => null()
        double precision, dimension(:,:,:,:), pointer :: array4d => null()
        integer :: array0i
        integer, dimension(:), pointer :: array1i => null()
        integer, dimension(:,:), pointer :: array2i => null()
        integer, dimension(:,:,:), pointer :: array3i => null()
        character (len=64) :: array0c
        character (len=64), dimension(:), pointer :: array1c => null()
        character (len=64), dimension(:,:), pointer :: array2c => null()
    end type input_field_type

    type input_global_attributes
        integer  :: num_atts = 1000
        character (len=64), dimension(1:1000) :: attnames
        real, dimension(1:1000) :: real_atts
        double precision, dimension(1:1000) :: double_atts
        integer, dimension(1:1000) :: int_atts
        character (len=64), dimension(1:1000) :: char_atts
        integer, dimension(1:1000) :: att_type
    end type input_global_attributes

    type input_map_attributes
        integer :: num_atts = 1000
        character (len=64), dimension(1:1000) :: attnames
        real, dimension(1:1000) :: real_atts
        double precision, dimension(1:1000) :: double_atts
        integer, dimension(1:1000) :: int_atts
        character (len=64), dimension(1:1000) :: char_atts
        integer, dimension(1:1000) :: att_type
    end type input_map_attributes

    integer, parameter :: FIELD_TYPE_UNSUPPORTED = -1, &
                          FIELD_TYPE_REAL        =  1, &
                          FIELD_TYPE_DOUBLE      =  2, &
                          FIELD_TYPE_INTEGER     =  3, &
                          FIELD_TYPE_CHARACTER   =  4


    contains


    integer function scan_input_open(filename, handle, nRecords) result(stat)

        implicit none

        character (len=*), intent(in) :: filename
        type (input_handle_type), intent(out) :: handle
        integer, intent(out), optional :: nRecords

 
        stat = 0

        stat = nf90_open(trim(filename), NF90_NOWRITE, handle % ncid)
        if (stat /= NF90_NOERR) then
            stat = 1
            return
        end if

        stat = nf90_inquire(handle % ncid, nVariables=handle % num_vars)
        if (stat /= NF90_NOERR) then
            stat = 1
            return
        end if

        allocate(handle % varids(handle % num_vars))

        stat = nf90_inq_varids(handle % ncid, handle % num_vars, handle % varids)
        if (stat /= NF90_NOERR) then
            stat = 1
            return
        end if

        stat = nf90_inquire(handle % ncid, unlimitedDimId=handle % unlimited_dimid)
        if (stat /= NF90_NOERR) then
            stat = 1
            return
        end if

        if (present(nRecords)) then
            stat = nf90_inquire_dimension(handle % ncid, handle % unlimited_dimid, len=nRecords)
            if (stat /= NF90_NOERR) then
                stat = 1
                return
            end if

            ! In case we have an input file that no time-varying records but
            ! does have time-invariant fields, set nRecords = 1 so that we can
            ! at least extract these fields
            if ((nRecords == 0) .and. (handle % num_vars > 0)) then
                nRecords = 1
            end if
        end if

        handle % current_var = 1

    end function scan_input_open


    integer function scan_input_close(handle) result(stat)

        implicit none

        type (input_handle_type), intent(inout) :: handle


        stat = 0

        stat = nf90_close(handle % ncid)
        if (stat /= NF90_NOERR) then
            stat = 1
        end if

        if (associated(handle % varids)) then
            deallocate(handle % varids)
        end if
        handle % current_var = 0

    end function scan_input_close


    integer function scan_input_rewind(handle) result(stat)

        implicit none

        type (input_handle_type), intent(inout) :: handle


        stat = 0

        handle % current_var = 1

    end function scan_input_rewind


    integer function scan_input_next_field(handle, field) result(stat)

        implicit none

        type (input_handle_type), intent(inout), target :: handle
        type (input_field_type), intent(out) :: field

        integer :: idim


        stat = 0

        if (handle % current_var < 1 .or. handle % current_var > handle % num_vars) then
            stat = 1
            return
        end if

        field % varid = handle % varids(handle % current_var)
        stat = nf90_inquire_variable(handle % ncid, field % varid, &
                                     name=field % name, &
                                     xtype=field % xtype, &
                                     ndims=field % ndims)
        if (stat /= NF90_NOERR) then
            stat = 1
            return
        end if

        if (field % xtype == NF90_FLOAT) then
            field % xtype = FIELD_TYPE_REAL
        else if (field % xtype == NF90_DOUBLE) then
            field % xtype = FIELD_TYPE_DOUBLE
        else if (field % xtype == NF90_INT) then
            field % xtype = FIELD_TYPE_INTEGER
        else if (field % xtype == NF90_CHAR) then
            field % xtype = FIELD_TYPE_CHARACTER
        else
            field % xtype = FIELD_TYPE_UNSUPPORTED
        end if

        allocate(field % dimids(field % ndims))

        stat = nf90_inquire_variable(handle % ncid, field % varid, &
                                     dimids=field % dimids)
        if (stat /= NF90_NOERR) then
            stat = 1
            deallocate(field % dimids)
            return
        end if

        allocate(field % dimlens(field % ndims))
        allocate(field % dimnames(field % ndims))

        do idim=1,field % ndims
            stat = nf90_inquire_dimension(handle % ncid, field % dimids(idim), &
                                          name=field % dimnames(idim), &
                                          len=field % dimlens(idim))
            if (field % dimids(idim) == handle % unlimited_dimid) then
                field % isTimeDependent = .true.
            end if
        end do

        field % file_handle => handle

        handle % current_var = handle % current_var + 1

    end function scan_input_next_field


    integer function scan_input_for_field(handle, fieldname, field) result(stat)

        implicit none

        type (input_handle_type), intent(inout), target :: handle
        character (len=*), intent(in) :: fieldname
        type (input_field_type), intent(out) :: field

        integer :: idim

        stat = 0

        stat = nf90_inq_varid(handle % ncid, trim(fieldname), field % varid)
        if (stat /= NF90_NOERR) then
            stat = 1
            return
        end if

        stat = nf90_inquire_variable(handle % ncid, field % varid, &
                                     name=field % name, &
                                     xtype=field % xtype, &
                                     ndims=field % ndims)
        if (stat /= NF90_NOERR) then
            stat = 1
            return
        end if

        if (field % xtype == NF90_FLOAT) then
            field % xtype = FIELD_TYPE_REAL
        else if (field % xtype == NF90_DOUBLE) then
            field % xtype = FIELD_TYPE_DOUBLE
        else if (field % xtype == NF90_INT) then
            field % xtype = FIELD_TYPE_INTEGER
        else if (field % xtype == NF90_CHAR) then
            field % xtype = FIELD_TYPE_CHARACTER
        else
            field % xtype = FIELD_TYPE_UNSUPPORTED
        end if

        allocate(field % dimids(field % ndims))

        stat = nf90_inquire_variable(handle % ncid, field % varid, &
                                     dimids=field % dimids)
        if (stat /= NF90_NOERR) then
            stat = 1
            deallocate(field % dimids)
            return
        end if

        allocate(field % dimlens(field % ndims))
        allocate(field % dimnames(field % ndims))

        do idim=1,field % ndims
            stat = nf90_inquire_dimension(handle % ncid, field % dimids(idim), &
                                          name=field % dimnames(idim), &
                                          len=field % dimlens(idim))
            if (field % dimids(idim) == handle % unlimited_dimid) then
                field % isTimeDependent = .true.
            end if
        end do

        field % file_handle => handle

    end function scan_input_for_field


    integer function scan_input_read_field(field, frame) result(stat)

        implicit none

        type (input_field_type), intent(inout) :: field
        integer, intent(in), optional :: frame

        integer :: local_frame
        integer, dimension(5) :: start, count
        real, dimension(1) :: temp1r
        double precision, dimension(1) :: temp1d
        integer, dimension(1) :: temp1i
        character (len=64), dimension(1) :: temp1c 

        stat = 0

        local_frame = 1
        if (present(frame)) then
            local_frame = frame
        end if

        if (field % xtype == FIELD_TYPE_REAL) then
            if (field % ndims == 0 .or. (field % ndims == 1 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = local_frame
                    count(1) = 1
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, temp1r, &
                                        start=start(1:1), count=count(1:1))
                    field % array0r = temp1r(1)
                else
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array0r)
                end if
            else if (field % ndims == 1 .or. (field % ndims == 2 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = local_frame
                    count(2) = 1
                    allocate(field % array1r(count(1)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array1r, &
                                        start=start(1:2), count=count(1:2))
                else
                    allocate(field % array1r(field%dimlens(1)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array1r)
                end if
            else if (field % ndims == 2 .or. (field % ndims == 3 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = 1
                    count(2) = field % dimlens(2)
                    start(3) = local_frame
                    count(3) = 1
                    allocate(field % array2r(count(1),count(2)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array2r, &
                                        start=start(1:3), count=count(1:3))
                else
                    allocate(field % array2r(field%dimlens(1),field%dimlens(2)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array2r)
                end if
            else if (field % ndims == 3 .or. (field % ndims == 4 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = 1
                    count(2) = field % dimlens(2)
                    start(3) = 1
                    count(3) = field % dimlens(3)
                    start(4) = local_frame
                    count(4) = 1
                    allocate(field % array3r(count(1),count(2),count(3)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array3r, &
                                        start=start(1:4), count=count(1:4))
                else
                    allocate(field % array3r(field%dimlens(1),field%dimlens(2),field%dimlens(3)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array3r)
                end if
            else if (field % ndims == 4 .or. (field % ndims == 5 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = 1
                    count(2) = field % dimlens(2)
                    start(3) = 1
                    count(3) = field % dimlens(3)
                    start(4) = 1
                    count(4) = field % dimlens(4)
                    start(5) = local_frame
                    count(5) = 1
                    allocate(field % array4r(count(1),count(2),count(3),count(4)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array4r, &
                                        start=start(1:5), count=count(1:5))
                else
                    allocate(field % array4r(field%dimlens(1),field%dimlens(2),field%dimlens(3),field%dimlens(4)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array4r)
                end if
            end if
        else if (field % xtype == FIELD_TYPE_DOUBLE) then
            if (field % ndims == 0 .or. (field % ndims == 1 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = local_frame
                    count(1) = 1
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, temp1d, &
                                        start=start(1:1), count=count(1:1))
                    field % array0d = temp1d(1)
                else
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array0d)
                end if
            else if (field % ndims == 1 .or. (field % ndims == 2 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = local_frame
                    count(2) = 1
                    allocate(field % array1d(count(1)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array1d, &
                                        start=start(1:2), count=count(1:2))
                else
                    allocate(field % array1d(field%dimlens(1)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array1d)
                end if
            else if (field % ndims == 2 .or. (field % ndims == 3 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = 1
                    count(2) = field % dimlens(2)
                    start(3) = local_frame
                    count(3) = 1
                    allocate(field % array2d(count(1),count(2)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array2d, &
                                        start=start(1:3), count=count(1:3))
                else
                    allocate(field % array2d(field%dimlens(1),field%dimlens(2)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array2d)
                end if
            else if (field % ndims == 3 .or. (field % ndims == 4 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = 1
                    count(2) = field % dimlens(2)
                    start(3) = 1
                    count(3) = field % dimlens(3)
                    start(4) = local_frame
                    count(4) = 1
                    allocate(field % array3d(count(1),count(2),count(3)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array3d, &
                                        start=start(1:4), count=count(1:4))
                else
                    allocate(field % array3d(field%dimlens(1),field%dimlens(2),field%dimlens(3)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array3d)
                end if
            else if (field % ndims == 4 .or. (field % ndims == 5 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = 1
                    count(2) = field % dimlens(2)
                    start(3) = 1
                    count(3) = field % dimlens(3)
                    start(4) = 1
                    count(4) = field % dimlens(4)
                    start(5) = local_frame
                    count(5) = 1
                    allocate(field % array4d(count(1),count(2),count(3),count(4)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array4d, &
                                        start=start(1:5), count=count(1:5))
                else
                    allocate(field % array4d(field%dimlens(1),field%dimlens(2),field%dimlens(3),field%dimlens(4)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array4d)
                end if
            end if

        else if (field % xtype == FIELD_TYPE_INTEGER) then
            if (field % ndims == 0 .or. (field % ndims == 1 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = local_frame
                    count(1) = 1
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, temp1i, &
                                        start=start(1:1), count=count(1:1))
                    field % array0i = temp1i(1)
                else
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array0i)
                end if
            else if (field % ndims == 1 .or. (field % ndims == 2 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = local_frame
                    count(2) = 1
                    allocate(field % array1i(count(1)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array1i, &
                                        start=start(1:2), count=count(1:2))
                else
                    allocate(field % array1i(field%dimlens(1)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array1i)
                end if
            else if (field % ndims == 2 .or. (field % ndims == 3 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = 1
                    count(2) = field % dimlens(2)
                    start(3) = local_frame
                    count(3) = 1
                    allocate(field % array2i(count(1),count(2)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array2i, &
                                        start=start(1:3), count=count(1:3))
                else
                    allocate(field % array2i(field%dimlens(1),field%dimlens(2)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array2i)
                end if
            else if (field % ndims == 3 .or. (field % ndims == 4 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = 1
                    count(2) = field % dimlens(2)
                    start(3) = 1
                    count(3) = field % dimlens(3)
                    start(4) = local_frame
                    count(4) = 1
                    allocate(field % array3i(count(1),count(2),count(3)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array3i, &
                                        start=start(1:4), count=count(1:4))
                else
                    allocate(field % array3i(field%dimlens(1),field%dimlens(2),field%dimlens(3)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array3i)
                end if
            end if

        else if (field % xtype == FIELD_TYPE_CHARACTER) then
            if (field % ndims == 0 .or. (field % ndims == 1 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = local_frame
                    count(1) = 1
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, temp1c, &
                                        start=start(1:1), count=count(1:1))
                    field % array0c = temp1c(1)
                else
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array0c)
                end if
            else if (field % ndims == 1 .or. (field % ndims == 2 .and. field % isTimeDependent)) then
                if (field % isTimeDependent) then
                    start(1) = 1
                    count(1) = field % dimlens(1)
                    start(2) = local_frame
                    count(2) = 1
                    allocate(field % array1c(count(1)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array1c, &
                                        start=start(1:2), count=count(1:2))
                else
                    allocate(field % array1c(field%dimlens(1)))
                    stat = nf90_get_var(field % file_handle % ncid, field % varid, field % array1c)
                end if
            end if

        end if

        if (stat /= NF90_NOERR) then
            write(0,*) ' '
            write(0,*) 'NetCDF error: reading '//trim(field % name)//' returned ', stat
            write(0,*) ' '
            stat = 1
        else
            stat = 0
        end if

    end function scan_input_read_field


    integer function scan_input_free_field(field) result(stat)

        implicit none

        type (input_field_type), intent(inout) :: field


        stat = 0

        if (associated(field % dimids)) then
            deallocate(field % dimids)
        end if
        if (associated(field % dimlens)) then
            deallocate(field % dimlens)
        end if
        if (associated(field % dimnames)) then
            deallocate(field % dimnames)
        end if

        if (associated(field % array1r)) then
            deallocate(field % array1r)
        end if
        if (associated(field % array2r)) then
            deallocate(field % array2r)
        end if
        if (associated(field % array3r)) then
            deallocate(field % array3r)
        end if
        if (associated(field % array4r)) then
            deallocate(field % array4r)
        end if

        if (associated(field % array1d)) then
            deallocate(field % array1d)
        end if
        if (associated(field % array2d)) then
            deallocate(field % array2d)
        end if
        if (associated(field % array3d)) then
            deallocate(field % array3d)
        end if
        if (associated(field % array4d)) then
            deallocate(field % array4d)
        end if

        if (associated(field % array1i)) then
            deallocate(field % array1i)
        end if
        if (associated(field % array2i)) then
            deallocate(field % array2i)
        end if
        if (associated(field % array3i)) then
            deallocate(field % array3i)
        end if

        nullify(field % file_handle)

    end function scan_input_free_field
   
    integer function scan_input_get_global_attributes(handle, globatts) result(stat)
 
        implicit none
 
        type (input_handle_type), intent(in) :: handle
        type (input_global_attributes), intent(out) :: globatts
 
        integer :: iatts
        integer, parameter :: FIELD_TYPE_UNSUPPORTED = -1, &
                              FIELD_TYPE_REAL        =  1, &
                              FIELD_TYPE_DOUBLE      =  2, &
                              FIELD_TYPE_INTEGER     =  3, &
                              FIELD_TYPE_CHARACTER   =  4
        stat = 0
 
        stat = nf90_inquire(handle % ncid, nAttributes=globatts % num_atts)
        if (stat /= NF90_NOERR) then
            stat = 1
            return
        end if
 
        do iatts=1,globatts % num_atts
            stat=nf90_inq_attname(handle % ncid, nf90_global, iatts, name=globatts % attnames(iatts))
            if (stat /= NF90_NOERR) then
                stat = 1
                return
            end if
 
            stat=nf90_inquire_attribute(handle % ncid, nf90_global, globatts % attnames(iatts), xtype=globatts % att_type(iatts))
            if (stat /= NF90_NOERR) then
                stat = 1
                return
            end if
 
 
            if (globatts % att_type(iatts) == NF90_FLOAT) then
                globatts % att_type(iatts) = FIELD_TYPE_REAL
                stat=nf90_get_att(handle % ncid, nf90_global, globatts % attnames(iatts), globatts % real_atts(iatts))
                if (stat /= NF90_NOERR) then
                    stat = 1
                    return
                end if
            else if (globatts % att_type(iatts) == NF90_DOUBLE) then
                globatts % att_type(iatts) = FIELD_TYPE_DOUBLE
                stat=nf90_get_att(handle % ncid, nf90_global, globatts % attnames(iatts), globatts % double_atts(iatts))
                if (stat /= NF90_NOERR) then
                    stat = 1
                    return
                end if
            else if (globatts % att_type(iatts) == NF90_INT) then
                globatts % att_type(iatts) = FIELD_TYPE_INTEGER
                stat=nf90_get_att(handle % ncid, nf90_global, globatts % attnames(iatts), globatts % int_atts(iatts))
                if (stat /= NF90_NOERR) then
                    stat = 1
                    return
                end if
            else if (globatts % att_type(iatts) == NF90_CHAR) then
                globatts % att_type(iatts) = FIELD_TYPE_CHARACTER
                stat=nf90_get_att(handle % ncid, nf90_global, globatts % attnames(iatts), globatts % char_atts(iatts))
                if (stat /= NF90_NOERR) then
                    stat = 1
                    return
                end if
            else
                globatts % att_type(iatts) = FIELD_TYPE_UNSUPPORTED
            end if
 
        end do

    end function scan_input_get_global_attributes


    integer function scan_input_get_map_attributes(filename, mapatts) result(stat)
 
        implicit none
 
        character (len=*), intent(in) :: filename
        type (input_map_attributes), intent(out) :: mapatts
 
        integer :: mapiatts, ncid
        integer, parameter :: FIELD_TYPE_UNSUPPORTED = -1, &
                              FIELD_TYPE_REAL        =  1, &
                              FIELD_TYPE_DOUBLE      =  2, &
                              FIELD_TYPE_INTEGER     =  3, &
                              FIELD_TYPE_CHARACTER   =  4

        integer :: attlen

        stat = 0
        stat = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        if (stat /= NF90_NOERR) then
            stat = 1
            return
        end if

        stat = nf90_inquire(ncid, nAttributes=mapatts % num_atts)
        if (stat /= NF90_NOERR) then
            stat = 1
            return
        end if

        do mapiatts=1,mapatts % num_atts

            stat=nf90_inq_attname(ncid, nf90_global, mapiatts, name=mapatts % attnames(mapiatts))
            if (stat /= NF90_NOERR) then
                stat = 1
                return
            end if
 
            stat=nf90_inquire_attribute(ncid, nf90_global, mapatts % attnames(mapiatts), xtype=mapatts % att_type(mapiatts), len=attlen)
            if (stat /= NF90_NOERR) then
                stat = 1
                return
            end if
 
            if (mapatts % att_type(mapiatts) == NF90_FLOAT .and. attlen == 1) then
                mapatts % att_type(mapiatts) = FIELD_TYPE_REAL
                stat=nf90_get_att(ncid, nf90_global, mapatts % attnames(mapiatts), mapatts % real_atts(mapiatts))
                if (stat /= NF90_NOERR) then
                    stat = 1 
                    return
                end if
            else if (mapatts % att_type(mapiatts) == NF90_DOUBLE .and. attlen == 1) then
                mapatts % att_type(mapiatts) = FIELD_TYPE_DOUBLE
                stat=nf90_get_att(ncid, nf90_global, mapatts % attnames(mapiatts), mapatts % double_atts(mapiatts))
                if (stat /= NF90_NOERR) then
                    stat = 1 
                    return
                end if
            else if (mapatts % att_type(mapiatts) == NF90_INT .and. attlen == 1) then
                mapatts % att_type(mapiatts) = FIELD_TYPE_INTEGER
                stat=nf90_get_att(ncid, nf90_global, mapatts % attnames(mapiatts), mapatts % int_atts(mapiatts))
                if (stat /= NF90_NOERR) then
                    stat = 1 
                    return
                end if
            else if (mapatts % att_type(mapiatts) == NF90_CHAR .and. attlen == 1) then
                mapatts % att_type(mapiatts) = FIELD_TYPE_CHARACTER
                stat=nf90_get_att(ncid, nf90_global, mapatts % attnames(mapiatts), mapatts % char_atts(mapiatts))
                if (stat /= NF90_NOERR) then
                    stat = 1 
                    return
                end if
            else
                mapatts % att_type(mapiatts) = FIELD_TYPE_UNSUPPORTED
            end if

        end do
    end function scan_input_get_map_attributes

end module scan_input
