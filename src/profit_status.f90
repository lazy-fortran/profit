module profit_status
    implicit none
    private
    integer, parameter, public :: profit_success = 0
    integer, parameter, public :: profit_invalid_shape = 1
    integer, parameter, public :: profit_invalid_input = 2
    integer, parameter, public :: profit_range_error = 3
    integer, parameter, public :: profit_io_error = 4
end module profit_status
