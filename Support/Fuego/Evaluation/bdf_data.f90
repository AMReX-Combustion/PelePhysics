module bdf_data
  use bdf, only : bdf_ts
  implicit none
  type(bdf_ts), save :: ts
  logical, save :: reuse_jac = .true.
  !$omp threadprivate(ts)
end module bdf_data

