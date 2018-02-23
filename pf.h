#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "rtklib.h"

typedef signed char        S8;                 /* Signed 8 bits integer    */
typedef unsigned char      U8;                 /* Unsigned 8 bits integer  */
typedef signed short       S16;                /* Signed 16 bits integer   */
typedef unsigned short     U16;                /* Unsigned 16 bits integer */
typedef signed long        S32;                /* Signed 32 bits integer   */
typedef unsigned long      U32;                /* Unsigned 32 bits integer */
typedef float              FLT;                /* 4 bytes floating point   */
typedef double             DBL;                /* 8 bytes floating point   */

typedef struct
{
	DBL    lat;
	DBL    lon;
	DBL    alt;
} LLA;

typedef struct
{
	DBL    X;
	DBL    Y;
	DBL    Z;
} ECEF_XYZ;

typedef struct
{
	U32         second;
	DBL         seconds;            // 

	DBL         pos[3];
	DBL         pos_sd[3];      // NED position standard deviations
	DBL         vel[3];         // NED velocities
	DBL         vel_sd[3];      // NED velocity standard deviations
	DBL         clkBias;
	DBL         clkDrift;
	DBL         PDOP;
	DBL         NDOP;           // North DOP
	DBL         EDOP;           // East DOP
	DBL         VDOP;           // Down DOP

	LLA         pos_lla;            // GPS position (lat,lon,alt)
	LLA         vel_lla;
	ECEF_XYZ    pos_ecef;       // ECEF position
	ECEF_XYZ    pos_ecef_sd;        // ECEF position standard deviations
	ECEF_XYZ    vel_ecef;       // ECEF velocities
	ECEF_XYZ    vel_ecef_sd;        // ECEF velocity standard deviations

	U8          inTunnel;
	U8          highWay;
	U8          pedastrian;

	U8          posBias;
	U8          velBias;

	//U8          pos_valid : 1;
	//U8          vel_valid : 1;
	//U8          valid : 1;
	//U8          data_valid : 1;
	U8          fix_source;

} navState_type;

typedef struct
{
	DBL    time;         // Receive time
	DBL    nSV;          // Number of Satellites
	DBL    GDOP;         // Dilutions of precision
	DBL    PDOP;         // Position dilution of precision
	DBL    HDOP;         // Horizontal dilution of precision
	DBL    NDOP;         // North dilution of precision
	DBL    EDOP;         // East dilution of precision
	DBL    VDOP;         // Vertical dilution of precision
	DBL    TDOP;         // Time dilution of precision
} navDOP;

int pf_main(const obsd_t *obs, int n, const nav_t *navs);
int pf_estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
	const double *vare, const int *svh, const nav_t *nav,
	const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,
	double *resp, char *msg);