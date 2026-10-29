module profit_output
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use profit_status
    implicit none
    private
    public :: write_study_csv
contains
    subroutine write_study_csv(path, metadata_keys, metadata_values, &
            column_names, values, status)
        character(*), intent(in) :: path, metadata_keys(:), metadata_values(:)
        character(*), intent(in) :: column_names(:)
        real(dp), intent(in) :: values(:, :)
        integer, intent(out) :: status
        character(:), allocatable :: line
        character(32) :: number
        integer :: unit, ios, close_status, row, col
        status = profit_invalid_shape
        if (size(metadata_keys) /= size(metadata_values)) return
        if (size(column_names) /= size(values, 1)) return
        if (size(column_names) == 0) return
        status = profit_invalid_input
        if (len_trim(path) == 0) return
        if (any(len_trim(column_names) == 0)) return
        if (any(len_trim(metadata_keys) == 0)) return
        if (.not. all(ieee_is_finite(values))) return
        status = profit_io_error
        open (newunit=unit, file=path, status="replace", action="write", iostat=ios)
        if (ios /= 0) return
        line = "record_type,metadata_key,metadata_value"
        do col = 1, size(column_names)
            line = line//","//quote(trim(column_names(col)))
        end do
        write (unit, '(a)', iostat=ios) line
        if (ios == 0) then
            do row = 1, size(metadata_keys)
                line = 'metadata,'//quote(trim(metadata_keys(row)))// &
                    ','//quote(trim(metadata_values(row)))
                do col = 1, size(column_names)
                    line = line//","
                end do
                write (unit, '(a)', iostat=ios) line
                if (ios /= 0) exit
            end do
        end if
        if (ios == 0) then
            do row = 1, size(values, 2)
                line = 'data,,'
                do col = 1, size(values, 1)
                    write (number, '(es25.17e3)') values(col, row)
                    line = line//","//trim(adjustl(number))
                end do
                write (unit, '(a)', iostat=ios) line
                if (ios /= 0) exit
            end do
        end if
        close (unit, iostat=close_status)
        if (ios /= 0 .or. close_status /= 0) return
        status = profit_success
    end subroutine write_study_csv

    function quote(value) result(escaped)
        character(*), intent(in) :: value
        character(:), allocatable :: escaped
        integer :: i
        escaped = '"'
        do i = 1, len(value)
            escaped = escaped//value(i:i)
            if (value(i:i) == '"') escaped = escaped//'"'
        end do
        escaped = escaped//'"'
    end function quote
end module profit_output
