module Gurobi_Interface_mod
    use iso_c_binding
    implicit none

    interface
        ! C prototype:
        ! void solve_set_cover_wrapper(int* n_cand, int* n_targ, int* row_ptr, int* col_idx, int* sol, int* stat)
        subroutine solve_set_cover_wrapper(n_cand, n_targ, row_ptr, col_idx, sol, stat) &
            bind(C, name="solve_set_cover_wrapper")
            
            import :: C_INT
            implicit none
            
            integer(C_INT), intent(in)    :: n_cand
            integer(C_INT), intent(in)    :: n_targ
            integer(C_INT), intent(in)    :: row_ptr(*) ! Array pointer
            integer(C_INT), intent(in)    :: col_idx(*) ! Array pointer
            integer(C_INT), intent(out)   :: sol(*)     ! Array pointer
            integer(C_INT), intent(out)   :: stat
        end subroutine solve_set_cover_wrapper
    end interface

end module Gurobi_Interface_mod