#include "udf.h"

/* ========================================================================== */
/* USER CONSTANTS                                                             */
/* ========================================================================== */
#define TARGET_MDOT    0.89017
#define MAX_DELTA      (0.050478 * TARGET_MDOT)
#define TOL_RATIO      0.05
#define RECOVERY_STEPS 20

#define OUTLET_ID      79
#define INLET_ID       82

/* ========================================================================== */
/* GLOBAL STATE                                                               */
/* ========================================================================== */
static int recovery_counter = 0;
static real current_mdot = TARGET_MDOT;
static int time_step_counter = 0;

/* ========================================================================== */
/* SAFE HELPER                                                                */
/* ========================================================================== */

static void
get_outlet_conditions(Domain *domain, real *mdot, real *vf)
{
    Thread *t = NULL;
    face_t f;

    real sum_mdot = 0.0;
    real sum_vf = 0.0;
    int count = 0;

    t = Lookup_Thread(domain, OUTLET_ID);

    if (!t)
    {
        *mdot = 0.0;
        *vf = 0.0;
        return;
    }

    begin_f_loop(f, t)
    {
        real flux = F_FLUX(f, t);
        real vof = 0.0;

        cell_t c0 = F_C0(f, t);
        Thread *t0 = THREAD_T0(t);

        if (c0 >= 0 && t0)
        {
            if (THREAD_STORAGE(t0, SV_VOF) != NULL)
                vof = C_VOF(c0, t0);
        }

        if (flux == flux)
            sum_mdot += fabs(flux);

        if (vof == vof)
            sum_vf += vof;

        count++;
    }
    end_f_loop(f, t)

#if RP_NODE
    sum_mdot = PRF_GRSUM1(sum_mdot);
    sum_vf   = PRF_GRSUM1(sum_vf);
    count    = PRF_GISUM1(count);
#endif

    if (count > 0)
    {
        *mdot = sum_mdot;
        *vf   = sum_vf / count;
    }
    else
    {
        *mdot = 0.0;
        *vf   = 0.0;
    }
}

/* ========================================================================== */
/* EXECUTE AT END: TEMPERATURE STATISTICS                                    */
/* ========================================================================== */

DEFINE_EXECUTE_AT_END(volume_temp_stddev_mass_csv)
{
    Domain *domain = Get_Domain(1);
    Thread *t;
    cell_t c;

    real sumT_mass = 0.0;
    real sumT2_mass = 0.0;
    real totalMass = 0.0;

    thread_loop_c(t, domain)
    {
        if (!FLUID_THREAD_P(t))
            continue;

        begin_c_loop(c, t)
        {
            real T = C_T(c, t);
            real rho = C_R(c, t);
            real vol = C_VOLUME(c, t);

            if (!(T == T && rho == rho && vol == vol))
                continue;

            real mass = rho * vol;

            if (mass <= 0.0 || !(mass == mass))
                continue;

            sumT_mass  += T * mass;
            sumT2_mass += T * T * mass;
            totalMass  += mass;
        }
        end_c_loop(c, t)
    }

#if RP_NODE
    sumT_mass  = PRF_GRSUM1(sumT_mass);
    sumT2_mass = PRF_GRSUM1(sumT2_mass);
    totalMass  = PRF_GRSUM1(totalMass);
#endif

    if (totalMass <= 0.0)
        return;

    real meanT  = sumT_mass / totalMass;
    real meanT2 = sumT2_mass / totalMass;

    real stddev = 0.0;
    if (meanT2 > meanT * meanT)
        stddev = sqrt(meanT2 - meanT * meanT);

    if (I_AM_NODE_ZERO_P)
    {
        FILE *fp = fopen("temp_std_mass.csv", (N_ITER == 1) ? "w" : "a");
        if (!fp) return;

        if (N_ITER == 1)
            fprintf(fp, "Iteration,MeanMassTemperature,StdDevMassTemperature\n");

        fprintf(fp, "%d,%g,%g\n", N_ITER, meanT, stddev);
        fclose(fp);
    }
}

/* ========================================================================== */
/* EXECUTE AT END: FLOW CONTROL (WITH STARTUP CONDITION)                     */
/* ========================================================================== */

DEFINE_EXECUTE_AT_END(flow_control_logic)
{
    Domain *domain = Get_Domain(1);

    real mdot_out = 0.0;
    real vf_out = 0.0;

    time_step_counter++;

    /* FIRST 3 STEPS */
    if (time_step_counter <= 3)
    {
        current_mdot = TARGET_MDOT;
        return;
    }

    get_outlet_conditions(domain, &mdot_out, &vf_out);

    if (!(mdot_out == mdot_out && vf_out == vf_out))
        return;

    if (mdot_out < 0.0)
        mdot_out = fabs(mdot_out);

    if (vf_out < 0.0) vf_out = 0.0;
    if (vf_out > 1.0) vf_out = 1.0;

    real gain = 0.02;
    real candidate;
    real max_allowed;

    candidate = mdot_out * (1.0 + gain);

    max_allowed = mdot_out + MAX_DELTA;

    if (candidate > max_allowed)
        candidate = max_allowed;

    current_mdot = 0.8 * current_mdot + 0.2 * candidate;

    if (current_mdot > max_allowed)
        current_mdot = max_allowed;

    if (current_mdot > TARGET_MDOT)
        current_mdot = TARGET_MDOT;

    if (!(vf_out >= 0.9999 && fabs(mdot_out - TARGET_MDOT) <= TOL_RATIO * TARGET_MDOT))
    {
        recovery_counter = RECOVERY_STEPS;
    }

    if (current_mdot < 1e-4 || !(current_mdot == current_mdot))
        current_mdot = TARGET_MDOT;
}

/* ========================================================================== */
/* PROFILE (FIXED STARTUP BUG)                                               */
/* ========================================================================== */

DEFINE_PROFILE(inlet2_massflow_profile, t, i)
{
    face_t f;
    real val;

    /* CRITICAL FIX */
    if (N_TIME < 3)
    {
        val = TARGET_MDOT;
    }
    else
    {
        val = current_mdot;
    }

    if (!(val == val) || val < 1e-4)
        val = TARGET_MDOT;

    begin_f_loop(f, t)
    {
        F_PROFILE(f, t, i) = val;
    }
    end_f_loop(f, t)
}