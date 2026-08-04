#include <stdlib.h>
#include "gurobi_c.h"

extern "C" {
    // Fortran interface wrapper
    void solve_set_cover_wrapper(int* n_cand, int* n_targ,
        int* row_ptr, int* col_idx,
        int* solution, int* status) {

        // ==========================================
        // 1. Variable Declaration
        // ==========================================
        GRBenv* env = NULL;
        GRBmodel* model = NULL;
        int       error = 0;
        int       cols = *n_cand;
        int       rows = *n_targ;
        int       optimstatus = 0;

        // Initialize pointers to NULL for safe cleanup via goto
        double* obj = NULL;
        char* vtype = NULL;
        int* cind = NULL;
        double* cval = NULL;
        double* x = NULL;

        // ==========================================
        // 2. Logic Start
        // ==========================================

        // --- Step A: Create Environment Silently ---
        // 1. Create an empty environment container first
        error = GRBemptyenv(&env);
        if (error) goto QUIT;

        // 2. Disable console output BEFORE starting the environment
        //    This suppresses the "Set parameter LicenseID..." message.
        if (!error) error = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);

        // 3. Disable gurobi.log file generation
        if (!error) error = GRBsetstrparam(env, GRB_STR_PAR_LOGFILE, "");

        if (error) goto QUIT;

        // 4. Now start the environment (License check happens here, but silently)
        error = GRBstartenv(env);
        if (error) goto QUIT;


        // --- Step B: Create Model ---
        error = GRBnewmodel(env, &model, "SetCover", 0, NULL, NULL, NULL, NULL, NULL);
        if (error) goto QUIT;

        // --- Step C: Add Variables ---
        obj = (double*)malloc(cols * sizeof(double));
        vtype = (char*)malloc(cols * sizeof(char));

        if (!obj || !vtype) { error = 10001; goto QUIT; } // Allocation error

        for (int j = 0; j < cols; j++) {
            obj[j] = 1.0;       // Objective coefficient (Minimize number of sensors)
            vtype[j] = GRB_BINARY;
        }

        error = GRBaddvars(model, cols, 0, NULL, NULL, NULL, obj, NULL, NULL, vtype, NULL);

        // Free temporary memory immediately
        free(obj); obj = NULL;
        free(vtype); vtype = NULL;

        if (error) goto QUIT;

        // --- Step D: Add Constraints (CSR Format) ---
        for (int i = 0; i < rows; i++) {
            // Convert Fortran 1-based indices to C 0-based indices
            int start = row_ptr[i] - 1;
            int end = row_ptr[i + 1] - 1;
            int len = end - start;

            if (len > 0) {
                cind = (int*)malloc(len * sizeof(int));
                cval = (double*)malloc(len * sizeof(double));

                if (!cind || !cval) { error = 10002; goto QUIT; }

                for (int k = 0; k < len; k++) {
                    cind[k] = col_idx[start + k] - 1; // 1-based to 0-based
                    cval[k] = 1.0;
                }

                // Constraint: Sum of selected sensors covering target i >= 1
                error = GRBaddconstr(model, len, cind, cval, GRB_GREATER_EQUAL, 1.0, NULL);

                free(cind); cind = NULL;
                free(cval); cval = NULL;

                if (error) goto QUIT;
            }
        }

        // --- Step E: Optimize ---
        error = GRBoptimize(model);
        if (error) goto QUIT;

        // --- Step F: Get Solution ---
        error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);

        if (optimstatus == GRB_OPTIMAL) {
            x = (double*)malloc(cols * sizeof(double));
            if (!x) { error = 10003; goto QUIT; }

            error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, cols, x);
            if (!error) {
                // Convert double (0.0/1.0) to int (0/1)
                for (int j = 0; j < cols; j++) solution[j] = (x[j] > 0.5) ? 1 : 0;
                *status = 0; // Success
            }
            free(x); x = NULL;
        }
        else {
            *status = 1; // Infeasible or other issue
        }

    QUIT:
        // Set return status code if an error occurred
        if (error) *status = error;

        // Unified resource cleanup
        // (Safe to call free on NULL or already freed pointers if set to NULL)
        if (obj) free(obj);
        if (vtype) free(vtype);
        if (cind) free(cind);
        if (cval) free(cval);
        if (x) free(x);

        if (model) GRBfreemodel(model);
        if (env)   GRBfreeenv(env);
    }
}