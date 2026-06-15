#include "udf.h"

DEFINE_EXECUTE_AT_END(volume_temp_stddev_mass_csv)
{
    Domain *domain = Get_Domain(1);
    Thread *t;
    cell_t c;

    real T, rho, vol, mass;
    real sumT_mass = 0.0;
    real sumT2_mass = 0.0;
    real totalMass = 0.0;
    real meanT, meanT2, stddev;

    /* ---- MASS-WEIGHTED MEAN AND STDDEV CALCULATIONS ---- */
    thread_loop_c(t, domain)
    {
        begin_c_loop(c, t)
        {
            T = C_T(c, t);
            rho = C_R(c, t);
            vol = C_VOLUME(c, t);
            mass = rho * vol;

            sumT_mass  += T * mass;
            sumT2_mass += T * T * mass;
            totalMass  += mass;
        }
        end_c_loop(c, t)
    }

#if RP_NODE
    /* ---- SYNCHRONIZATION IN PARALLEL MODE ---- */
    sumT_mass  = PRF_GRSUM1(sumT_mass);
    sumT2_mass = PRF_GRSUM1(sumT2_mass);
    totalMass  = PRF_GRSUM1(totalMass);
#endif

    if (totalMass <= 0.0)
        return;

    meanT  = sumT_mass  / totalMass;
    meanT2 = sumT2_mass / totalMass;
    stddev = sqrt(meanT2 - meanT * meanT);

    /* ---- WRITE TO CSV ONLY ON NODE 0 ---- */
    if (I_AM_NODE_ZERO_P)
    {
        static FILE *fp = NULL;
        static int header_written = 0;

        if (fp == NULL)
        {
            fp = fopen("temp_std_mass.csv", "w");
            if (fp == NULL)
            {
                Message("ERROR: cannot create temp_std_mass.csv\n");
                return;
            }
        }

        /* ---- HEADER ONLY ONCE ---- */
        if (!header_written)
        {
            fprintf(fp, "Iteration,MeanMassTemperature,StdDevMassTemperature\n");
            header_written = 1;
        }

        /* ---- WRITE TO CSV ---- */
        fprintf(fp, "%d,%g,%g\n", N_ITER, meanT, stddev);
        fflush(fp);
    }

    /* ---- CONSOLE MESSAGE ---- */
    Message("Iter %d | Mean_m = %g K | StdDev_m = %g K\n", N_ITER, meanT, stddev);
}
