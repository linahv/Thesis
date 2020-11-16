program alpha_program

use bmad
use alpha_modules

implicit none

type (lat_struct), target :: lat   ! This structure holds the lattice info
type (ele_struct), pointer :: ele
type (ele_struct) :: ele2

! Printing parameters
integer i, j, nmax
REAL(rp) start_s, end_s, delta_s ,s
INTEGER n_slices
logical err

!Calculations parameters
integer N 
type(coord_struct), target, allocatable :: orb(:)
REAL(rp), ALLOCATABLE :: F0(:)
REAL(rp), ALLOCATABLE :: F1(:)
REAL(rp), ALLOCATABLE :: F2(:)
REAL(rp) alpha0, nux, alpha1, alpha2

! Programs should always implement "intelligent bookkeeping".
bmad_com%auto_bookkeeper = .false.

N = 5000
n_slices = 100

ALLOCATE(F0(0:N-1))
ALLOCATE(F1(0:N-1))
ALLOCATE(F2(0:N-1))

call bmad_parser ("lat.lat", lat)  ! Read in a lattice.

! Propagate Twiss parameters
if (lat%param%geometry == closed$) call twiss_at_start (lat)
call twiss_propagate_all (lat)  
    
nmax = lat%n_ele_max
nux = lat%ele(lat%n_ele_track)%a%phi /2./pi

call FNdisp(lat, N, n_slices, nux, F0, F1)
call F2N(lat, N, n_slices, nux, F0, F1, F2)

call alpha0calc(lat, F0, F1, F2, N, n_slices, nux, alpha0, alpha1, alpha2)

!call getalpha(lat, N, n_slices, nux, alpha0, alpha1, alpha2)

print *, alpha0
print *, alpha1
print *, alpha2

endprogram
