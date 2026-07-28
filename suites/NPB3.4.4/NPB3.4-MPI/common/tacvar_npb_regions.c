/**
 * @file tacvar_npb_regions.c
 * @brief Logical region registry for every NPB-MPI named timer + step.
 *
 * Legacy detail regions are documented for NPB_TIMER_FLAG reporting but are
 * not written to TacVar event CSVs (recorded_by_tacvar=0).
 */
#include "tacvar_npb_regions.h"
#include "tacvar_measure.h"
#include <string.h>
#include <strings.h>

#ifndef TACVAR_REGION_STEP
#define TACVAR_REGION_STEP 1000
#endif

#define R_ALWAYS      "always"
#define R_FLAG        "NPB_TIMER_FLAG"
#define R_STEP        "TACVAR_ENABLE_PER_STEP_TIMING=1"

#define REG(id, name, loc, desc, when, rec) \
    { (id), (name), (loc), (desc), (when), (rec) }

static void register_list(const tacvar_region_info_t *list, size_t n)
{
    tacvar_region_info_clear();
    tacvar_region_info_register(list, n);
}

static void register_mg(void)
{
    static const tacvar_region_info_t r[] = {
        REG(1, "T_bench", "MG/mg.f90:mg_mpi:timed-section",
            "Full timed MG section (pre-resid + all V-cycles + post-norm)",
            R_ALWAYS, 1),
        REG(2, "T_init", "MG/mg.f90:mg_mpi:init",
            "Initialization and setup before timed section",
            R_ALWAYS, 0),
        REG(3, "t_psinv", "MG/mg.f90:psinv",
            "Approximate inverse smoother", R_FLAG, 0),
        REG(4, "t_resid", "MG/mg.f90:resid",
            "Residual computation", R_FLAG, 0),
        REG(5, "t_rprj3", "MG/mg.f90:rprj3",
            "Restriction of residual to coarser grid", R_FLAG, 0),
        REG(6, "t_interp", "MG/mg.f90:interp",
            "Prolongation / interpolation of correction", R_FLAG, 0),
        REG(7, "t_norm2u3", "MG/mg.f90:norm2u3",
            "L2 / infinity norm reduction", R_FLAG, 0),
        REG(8, "t_comm3", "MG/mg.f90:comm3",
            "Halo exchange communication", R_FLAG, 0),
        REG(9, "t_rcomm", "MG/mg.f90:norm2u3",
            "Norm reduction communication", R_FLAG, 0),
        REG(TACVAR_REGION_STEP, "step", "MG/mg.f90:mg_mpi:it-loop",
            "One MG outer iteration: mg3P V-cycle + resid",
            R_STEP, 1),
    };
    register_list(r, sizeof(r) / sizeof(r[0]));
}

static void register_cg(void)
{
    static const tacvar_region_info_t r[] = {
        REG(1, "t_total", "CG/cg.f90:main:timed-section",
            "All inverse-power iterations", R_ALWAYS, 1),
        REG(2, "t_conjg", "CG/cg.f90:conj_grad",
            "Conjugate-gradient solve", R_FLAG, 0),
        REG(3, "t_rcomm", "CG/cg.f90:conj_grad",
            "CG reduction communication", R_FLAG, 0),
        REG(4, "t_ncomm", "CG/cg.f90:main",
            "Norm / zeta reduction communication", R_FLAG, 0),
        REG(TACVAR_REGION_STEP, "step", "CG/cg.f90:main:it-loop",
            "One inverse-power iteration: conj_grad + normalize",
            R_STEP, 1),
    };
    register_list(r, sizeof(r) / sizeof(r[0]));
}

static void register_sp(void)
{
    static const tacvar_region_info_t r[] = {
        REG(1, "t_total", "SP/sp.f90:main:timed-section",
            "All ADI time steps", R_ALWAYS, 1),
        REG(2, "t_rhs", "SP/rhs.f90:rhs",
            "Right-hand side computation", R_FLAG, 0),
        REG(3, "t_xsolve", "SP/x_solve.f90:x_solve",
            "X-direction approximate factorization", R_FLAG, 0),
        REG(4, "t_ysolve", "SP/y_solve.f90:y_solve",
            "Y-direction approximate factorization", R_FLAG, 0),
        REG(5, "t_zsolve", "SP/z_solve.f90:z_solve",
            "Z-direction approximate factorization", R_FLAG, 0),
        REG(6, "t_bpack", "SP/copy_faces.f90:copy_faces",
            "Boundary pack / unpack", R_FLAG, 0),
        REG(7, "t_exch", "SP/copy_faces.f90:copy_faces",
            "Face exchange communication", R_FLAG, 0),
        REG(8, "t_xcomm", "SP/x_solve.f90:x_solve",
            "X-solve communication", R_FLAG, 0),
        REG(9, "t_ycomm", "SP/y_solve.f90:y_solve",
            "Y-solve communication", R_FLAG, 0),
        REG(10, "t_zcomm", "SP/z_solve.f90:z_solve",
            "Z-solve communication", R_FLAG, 0),
        REG(TACVAR_REGION_STEP, "step", "SP/sp.f90:main:step-loop",
            "One ADI time step (copy_faces + solves + add)",
            R_STEP, 1),
    };
    register_list(r, sizeof(r) / sizeof(r[0]));
}

static void register_bt(void)
{
    static const tacvar_region_info_t r[] = {
        REG(1, "t_total", "BT/bt.f90:main:timed-section",
            "All ADI time steps (includes optional I/O)", R_ALWAYS, 1),
        REG(2, "t_io", "BT/bt.f90:main",
            "Checkpoint / timestep I/O", R_FLAG, 0),
        REG(3, "t_rhs", "BT/rhs.f90:rhs",
            "Right-hand side computation", R_FLAG, 0),
        REG(4, "t_xsolve", "BT/x_solve*.f90:x_solve",
            "X-direction solve", R_FLAG, 0),
        REG(5, "t_ysolve", "BT/y_solve*.f90:y_solve",
            "Y-direction solve", R_FLAG, 0),
        REG(6, "t_zsolve", "BT/z_solve*.f90:z_solve",
            "Z-direction solve", R_FLAG, 0),
        REG(7, "t_bpack", "BT/copy_faces.f90",
            "Boundary pack / unpack", R_FLAG, 0),
        REG(8, "t_exch", "BT/copy_faces.f90",
            "Face exchange communication", R_FLAG, 0),
        REG(9, "t_xcomm", "BT/x_solve*.f90",
            "X-solve communication", R_FLAG, 0),
        REG(10, "t_ycomm", "BT/y_solve*.f90",
            "Y-solve communication", R_FLAG, 0),
        REG(11, "t_zcomm", "BT/z_solve*.f90",
            "Z-solve communication", R_FLAG, 0),
        REG(12, "t_enorm", "BT/error.f90",
            "Error-norm computation / communication", R_FLAG, 0),
        REG(13, "t_iov", "BT/verify.f90",
            "Verification I/O", R_FLAG, 0),
        REG(TACVAR_REGION_STEP, "step", "BT/bt.f90:main:step-loop",
            "One BT outer time step including optional I/O",
            R_STEP, 1),
    };
    register_list(r, sizeof(r) / sizeof(r[0]));
}

static void register_lu(void)
{
    static const tacvar_region_info_t r[] = {
        REG(1, "t_total", "LU/ssor.f90:ssor:timed-section",
            "All SSOR time steps", R_ALWAYS, 1),
        REG(2, "t_rhs", "LU/rhs.f90:rhs",
            "Right-hand side computation", R_FLAG, 0),
        REG(3, "t_blts", "LU/blts.f90:blts;LU/ssor.f90",
            "Lower triangular solution", R_FLAG, 0),
        REG(4, "t_buts", "LU/buts.f90:buts;LU/ssor.f90",
            "Upper triangular solution", R_FLAG, 0),
        REG(5, "t_jacld", "LU/jacld.f90:jacld",
            "Lower triangular Jacobian formation", R_FLAG, 0),
        REG(6, "t_jacu", "LU/jacu.f90:jacu",
            "Upper triangular Jacobian formation", R_FLAG, 0),
        REG(7, "t_exch", "LU/rhs.f90:rhs",
            "RHS face exchange", R_FLAG, 0),
        REG(8, "t_lcomm", "LU/ssor.f90;LU/blts*.f90",
            "Lower-sweep communication", R_FLAG, 0),
        REG(9, "t_ucomm", "LU/ssor.f90;LU/buts*.f90",
            "Upper-sweep communication", R_FLAG, 0),
        REG(10, "t_rcomm", "LU/l2norm.f90:l2norm",
            "Residual norm communication", R_FLAG, 0),
        REG(TACVAR_REGION_STEP, "step", "LU/ssor.f90:ssor:istep-loop",
            "One complete SSOR timestep",
            R_STEP, 1),
    };
    register_list(r, sizeof(r) / sizeof(r[0]));
}

static void register_ft(void)
{
    static const tacvar_region_info_t r[] = {
        REG(1, "T_total", "FT/ft.f90:main:timed-section",
            "Setup + forward FFT + all iterations + verify", R_ALWAYS, 1),
        REG(2, "T_setup", "FT/ft.f90:main",
            "Index map / initial conditions / FFT init", R_FLAG, 0),
        REG(3, "T_fft", "FT/ft.f90:fft",
            "3D FFT (including transpose)", R_FLAG, 0),
        REG(4, "T_evolve", "FT/ft.f90:evolve",
            "Time-evolution multiply by twiddle", R_FLAG, 0),
        REG(5, "T_checksum", "FT/ft.f90:checksum",
            "Checksum reduction", R_FLAG, 0),
        REG(6, "T_fftlow", "FT/ft.f90:cfft*",
            "Low-level 1D FFT", R_FLAG, 0),
        REG(7, "T_fftcopy", "FT/ft.f90:cfft*",
            "FFT local copy", R_FLAG, 0),
        REG(8, "T_transpose", "FT/ft.f90:transpose*",
            "Global transpose communication", R_FLAG, 0),
        REG(9, "T_transxzloc", "FT/ft.f90",
            "XZ transpose local phase", R_FLAG, 0),
        REG(10, "T_transxzglo", "FT/ft.f90",
            "XZ transpose global MPI phase", R_FLAG, 0),
        REG(11, "T_transxzfin", "FT/ft.f90",
            "XZ transpose finalize phase", R_FLAG, 0),
        REG(12, "T_transxyloc", "FT/ft.f90",
            "XY transpose local phase", R_FLAG, 0),
        REG(13, "T_transxyglo", "FT/ft.f90",
            "XY transpose global MPI phase", R_FLAG, 0),
        REG(14, "T_transxyfin", "FT/ft.f90",
            "XY transpose finalize phase", R_FLAG, 0),
        REG(15, "T_synch", "FT/ft.f90",
            "Explicit synchronization", R_FLAG, 0),
        REG(16, "T_init", "FT/ft.f90:main",
            "Untimed initialization before T_total", R_ALWAYS, 0),
        REG(TACVAR_REGION_STEP, "step", "FT/ft.f90:main:iter-loop",
            "One FT iteration: evolve + inverse FFT + checksum",
            R_STEP, 1),
    };
    register_list(r, sizeof(r) / sizeof(r[0]));
}

static void register_is(void)
{
    static const tacvar_region_info_t r[] = {
        REG(0, "T_TOTAL", "IS/is.c:main:timed-section",
            "All ranking iterations", R_ALWAYS, 1),
        REG(1, "T_RANK", "IS/is.c:rank",
            "Key ranking / bucket work", R_FLAG, 0),
        REG(2, "T_RCOMM", "IS/is.c:rank",
            "Ranking communication", R_FLAG, 0),
        REG(3, "T_VERIFY", "IS/is.c",
            "Verification helpers", R_FLAG, 0),
        REG(TACVAR_REGION_STEP, "step", "IS/is.c:main:iteration-loop",
            "One rank(iteration) pass",
            R_STEP, 1),
    };
    register_list(r, sizeof(r) / sizeof(r[0]));
}

static void register_ep(void)
{
    static const tacvar_region_info_t r[] = {
        REG(1, "t_total", "EP/ep.f90:main:timed-section",
            "Setup + all local batches + final Allreduce", R_ALWAYS, 1),
        REG(2, "t_gpairs", "EP/ep.f90:main",
            "Gaussian acceptance-rejection pairs", R_FLAG, 0),
        REG(3, "t_randn", "EP/ep.f90:main",
            "Uniform random number generation (vranlc)", R_FLAG, 0),
        REG(4, "t_rcomm", "EP/ep.f90:main",
            "Final sx/sy/q Allreduce", R_FLAG, 0),
        REG(TACVAR_REGION_STEP, "step", "EP/ep.f90:main:k-loop",
            "One local batch: seed + vranlc + Gaussian pairs",
            R_STEP, 1),
    };
    register_list(r, sizeof(r) / sizeof(r[0]));
}

static void register_dt(void)
{
    static const tacvar_region_info_t r[] = {
        REG(0, "T_TOTAL", "DT/dt.c:main;DT/dt.c:ProcessNodes",
            "All-rank ProcessNodes total (TacVar); NPB report uses rank-0 only",
            R_ALWAYS, 1),
        REG(TACVAR_REGION_STEP, "step", "DT/dt.c:ProcessNodes:owned-node",
            "One owned Source/Comparator/Sink node body",
            R_STEP, 1),
    };
    register_list(r, sizeof(r) / sizeof(r[0]));
}

void tacvar_npb_register_regions(const char *benchmark)
{
    if (!benchmark)
        return;
    if (!strcasecmp(benchmark, "mg"))
        register_mg();
    else if (!strcasecmp(benchmark, "cg"))
        register_cg();
    else if (!strcasecmp(benchmark, "sp"))
        register_sp();
    else if (!strcasecmp(benchmark, "bt"))
        register_bt();
    else if (!strcasecmp(benchmark, "lu"))
        register_lu();
    else if (!strcasecmp(benchmark, "ft"))
        register_ft();
    else if (!strcasecmp(benchmark, "is"))
        register_is();
    else if (!strcasecmp(benchmark, "ep"))
        register_ep();
    else if (!strcasecmp(benchmark, "dt"))
        register_dt();
}
