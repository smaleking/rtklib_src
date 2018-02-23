/* =============================================================================
* function name : pf.c
*        author : 
*          date : 
*
*   description : this file defines the function of position fix manager (pf)
* =============================================================================*/

#include "rtklib.h"
#include "pf.h"
/* constants and macros ------------------------------------------------------*/

#define SQR(x)   ((x)*(x))

#define RE_GLO   6378136.0        /* radius of earth (m)            ref [2] */
#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */
#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */

#define OMGE_GLO 7.292115E-5      /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */

#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */

#define ERREPH_GLO 5.0            /* error of glonass ephemeris (m) */
#define TSTEP    60.0             /* integration step glonass ephemeris (s) */
#define RTOL_KEPLER 1E-14         /* relative tolerance for Kepler equation */

#define DEFURASSR 0.15            /* default accurary of ssr corr (m) */
#define MAXECORSSR 10.0           /* max orbit correction of ssr (m) */
#define MAXCCORSSR (1E-6*CLIGHT)  /* max clock correction of ssr (m) */
#define MAXAGESSR 70.0            /* max age of ssr orbit and clock (s) */
#define MAXAGESSR_HRCLK 10.0      /* max age of ssr high-rate clock (s) */
#define STD_BRDCCLK 30.0          /* error of broadcast clock (m) */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_CBIAS   0.3         /* code bias error std (m) */

/*validate measurements based on 
** 1) Elevation,
** 2) CN0, e.g., 15.6 dB-Hz
** 3) Time in track, e.g., 0.2 s,
** 4) FLI (frequency lock indicator) or PLI,
** 5) Multipath indicator, e.g., slope, symmetric,
** 6) PR is longer than 80 ms
*/
static void pf_meas_validation(const obsd_t *obs)
{

}

/*select measurements based on residuals, etc.*/
static void pf_meas_selection(const obsd_t *obs)
{

}

/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
void pf_ecef2pos(const double *r, double *pos)
{
	double e2 = FE_WGS84*(2.0 - FE_WGS84), r2 = dot(r, r, 2), z, zk, v = RE_WGS84, sinp;

	for (z = r[2], zk = 0.0; fabs(z - zk) >= 1E-4;) {
		zk = z;
		sinp = z / sqrt(r2 + z*z);
		v = RE_WGS84 / sqrt(1.0 - e2*sinp*sinp);
		z = r[2] + v*e2*sinp;
	}
	pos[0] = r2>1E-12 ? atan(z / sqrt(r2)) : (r[2]>0.0 ? PI / 2.0 : -PI / 2.0);
	pos[1] = r2>1E-12 ? atan2(r[1], r[0]) : 0.0;
	pos[2] = sqrt(r2 + z*z) - v;
}

/* construct PR meas based on ToT ------------------------------------------
* 
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
void pf_construct_meas(int iter, const obsd_t *obs)
{
}

static int pf_calculate_res(int iter, const obsd_t *obs, int n, const double *rs,
	                         const double *dts, const double *vare, const int *svh,
	                         const nav_t *nav, const double *x, const prcopt_t *opt,
	                         double *v, double *H, double *var, double *azel, int *vsat,
	                         double *resp, int *nx)
{
	double r, dion, dtrp, vmeas, vion, vtrp, rr[3], pos[3], dtr, e[3], P;
	int i, j, nv = 0, ns[2] = { 0 }, sys;

	for (i = 0; i<3; i++) 
		rr[i] = x[i]; 
	
	 dtr = x[3];

	ecef2pos(rr, pos);

	for (i = 0; i < n&&i < MAXOBS; i++) 
	{
		/* geometric distance/azimuth/elevation angle */
		if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 ||
			satazel(pos, e, azel + i * 2)<0.1) continue;

		/* psudorange with code bias correction */
		//if ((P = prange(obs + i, nav, azel + i * 2, iter, opt, &vmeas)) == 0.0) continue;
		P = obs[i].P[0];
		
		/* ionospheric corrections */
		if (!ionocorr(obs[i].time, nav, obs[i].sat, pos, azel + i * 2,
			iter > 0 ? 0 : IONOOPT_BRDC, &dion, &vion)) {
			continue;
		}

		/* tropospheric corrections */
		if (!tropcorr(obs[i].time, nav, pos, azel + i * 2,
			iter > 0 ? 0 : TROPOPT_SAAS, &dtrp, &vtrp)) {
			continue;
		};

		/* pseudorange residual */
		v[nv] = P - (r + dtr - CLIGHT*dts[i * 2] + dion + dtrp);

		/* design matrix */
		for (j = 0; j<4; j++) H[j + nv * 4] = j<3 ? -e[j] : 1.0;

		vsat[i] = 1; resp[i] = v[nv];

		/* error variance */
		var[nv++] = 1.0;

		//for (i = 0; i<nv; i++) for (j = 0; j<4; j++) H[j + i * 4] = H[j + i * 5];

	}

	return nv;
}

/* Estimate receiver position - 3D ------------------------------------------------*/
static int pf_pos_lsq_3d(const double *A, const double *y, int n, int m, double *x,
	                     double *Q)
{
	double *Ay;
	int info;

	if (m<n) return -1;
	Ay = mat(n, 1);
	matmul("NN", n, 1, m, 1.0, A, y, 0.0, Ay); /* Ay=A*y */
	matmul("NT", n, n, m, 1.0, A, A, 0.0, Q);  /* Q=A*A' */
	if (!(info = matinv(Q, n))) matmul("NN", n, 1, n, 1.0, Q, Ay, 0.0, x); /* x=Q^-1*Ay */
	free(Ay);
	return info;
}

/* Estimate receiver velocity - 3D ------------------------------------------------*/
static int pf_vel_lsq_3d(const obsd_t *obs, int n, const double *rs, const double *dts,
	                     const double *vare, const int *svh, const nav_t *nav,
	                     const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,
	                     double *resp, char *msg)
{

}

/* estimate receiver position ------------------------------------------------*/
int pf_estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
	const double *vare, const int *svh, const nav_t *nav,
	const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,
	double *resp, char *msg)
{
	double x[4] = { 0 }, dx[4], Q[16], *v, *H, *var, sig;
	int i, j, k, info, stat, nx, nv;

	nx = 4;

	v = mat(n, 1); H = mat(4, n); var = mat(n, 1);

	/**/

	x[0] = -1200000;
	x[1] = -4700000;
	x[2] = 4070000;
	x[3] = 0;
	
	//for (i = 0; i<3; i++) x[i] = sol->rr[i];

	for (i = 0; i<MAXITR; i++) {

		/* pseudorange residuals */
		nv = pf_calculate_res(i, obs, n, rs, dts, vare, svh, nav, x, opt, v, H, var, azel, vsat, resp, &nx);

		if (nv<nx) 
		{
			sprintf(msg, "lack of valid sats ns=%d", nv);
			break;
		}

		/* weight by variance */
		for (j = 0; j<n; j++) 
		{
			//sig = sqrt(var[j]);
			sig = 1.0;
			v[j] /= sig;
			for (k = 0; k<n; k++) 
				H[k + j*n] /= sig;
		}
		/* least square estimation */
		if ((info = pf_pos_lsq_3d(H, v, 4, nv, dx, Q))) {
			sprintf(msg, "lsq error info=%d", info);
			break;
		}
		for (j = 0; j<4; j++) x[j] += dx[j];

		if (norm(dx, n)<1E-4) {
			sol->type = 0;
			sol->time = timeadd(obs[0].time, -x[3] / CLIGHT);
			sol->dtr[0] = x[3] / CLIGHT; /* receiver clock bias (s) */
			sol->dtr[1] = x[4] / CLIGHT; /* glonass-gps time offset (s) */
			for (j = 0; j<6; j++) sol->rr[j] = j<3 ? x[j] : 0.0;
			for (j = 0; j<3; j++) sol->qr[j] = (float)Q[j + j*nx];
			sol->qr[3] = (float)Q[1];    /* cov xy */
			sol->qr[4] = (float)Q[2 + nx]; /* cov yz */
			sol->qr[5] = (float)Q[2];    /* cov zx */
			sol->ns = (unsigned char)nv;
			sol->age = sol->ratio = 0.0;

			/* validate solution */
			//if ((stat = valsol(azel, vsat, n, opt, v, nv, nx, msg))) {
			//	sol->stat = opt->sateph == EPHOPT_SBAS ? SOLQ_SBAS : SOLQ_SINGLE;
			//}
			free(v); free(H); free(var);
			stat = 1;
			return stat;
		}
	}
	if (i >= MAXITR) sprintf(msg, "iteration divergent i=%d", i);

	free(v); free(H); free(var);

	return 0;
}

int pf_main(const obsd_t *obs, int n, const nav_t *navs)
{
	navState_type navState = {0};

	int i, stat, vsat[MAXOBS] = { 0 }, svh[MAXOBS];
	gtime_t time;
	sol_t *sol = {0};
	const nav_t *nav;
	double *rs, *dts, *var, *azel_, *resp;
	char msg[128] = "";

	/*initialization*/
	navState.pos_lla.lat = 0;
	navState.pos_lla.lon = 0;
	navState.pos_lla.alt = 0;
	//pos2ecef(navState.pos_lla, navState.pos_ecef);

	rs = mat(6, n); dts = mat(2, n); var = mat(1, n); azel_ = zeros(2, n); resp = mat(1, n);
	//sol->time = obs[0].time;
	time = obs[0].time;

	/* satellite positons, velocities and clocks */
	satposs(time, obs, n, navs, 0, rs, dts, var, svh);
	
	/* estimate receiver po sition with pseudorange */
	stat = pf_estpos(obs, n, rs, dts, var, svh, navs, 0, sol, azel_, vsat, resp, msg);
}