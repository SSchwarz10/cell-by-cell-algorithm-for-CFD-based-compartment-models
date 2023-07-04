/*This UDF is used to calculate the different source terms for the mean age moments. Every UDS has its own source term and needs the previous UDS
The source terms must be loaded individually for every moment/UDS. Since the solver in FLUENT considers the density, the source terms must be multiplied with the density*/
#include "udf.h"
#include "mem.h"
#include "models.h"
 

DEFINE_SOURCE(Moment_1_Source, cell, thread, dS, eqn)
{
	/*Source for moment 1 (UDS Index 0)*/
	real source;
	source = C_R(cell, thread) * 1;
	dS[eqn] = 0;
	return source;
}

DEFINE_SOURCE(Moment_2_Source, cell, thread, dS, eqn)
{
	/*Source for moment 2 (UDS Index 1)*/
	real source;
	source = C_R(cell, thread) * 2 * C_UDSI(cell, thread, 0);
	dS[eqn] = 0;
	return source;
}

DEFINE_SOURCE(Moment_3_Source, cell, thread, dS, eqn)
{
	/*Source for moment 3 (UDS Index 2)*/
	real source;
	source = C_R(cell, thread) * 3 * C_UDSI(cell, thread, 1);
	dS[eqn] = 0;
	return source;
}

DEFINE_SOURCE(Moment_4_Source, cell, thread, dS, eqn)
{
	/*Source for moment 4 (UDS Index 3)*/
	real source;
	source = C_R(cell, thread) * 4 * C_UDSI(cell, thread, 2);
	dS[eqn] = 0;
	return source;
}

DEFINE_SOURCE(Moment_5_Source, cell, thread, dS, eqn)
{
	/*Source for moment 5 (UDS Index 4)*/
	real source;
	source = C_R(cell, thread) * 5 * C_UDSI(cell, thread, 3);
	dS[eqn] = 0;
	return source;
}

DEFINE_SOURCE(Variance_Source, cell, thread, dS, eqn)
{
	/*Source for Variance (UDS Index 5)*/
	real source;
	source = 2 * C_UDSI_DIFF(cell, thread, 5) * NV_MAG2(C_UDSI_G(cell, thread, 0));
	dS[eqn] = 0;
	return source;
}

/*diffusion coefficient for every moment/UDS*/
DEFINE_DIFFUSIVITY(Diff_Moments, cell, thread, i)
{
	/*Diffusivity coefficient in FLUENT is defined with the density (kg/m^3/s)*/
	real Diff_l = 1e-09; /*m^2/s*/
	return C_R(cell, thread) * Diff_l;
}

