#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define PI          3.1415926535897932  /* pi */
#define D2R         (PI/180.0)          /* deg to rad */
#define R2D         (180.0/PI)          /* rad to deg */
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define SC2RAD      3.1415926535898     /* semi-circle to radian (IS-GPS) */
#define AU          149597870691.0      /* 1 AU (m) */
#define AS2R        (D2R/3600.0)        /* arc sec to radian */

#define OMGE        7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */

#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */

#define HION        350000.0            /* ionosphere height (m) */

#define MAXFREQ     7                   /* max NFREQ */

#define FREQ1       1.57542E9           /* L1/E1  frequency (Hz) */
#define FREQ2       1.22760E9           /* L2     frequency (Hz) */
#define FREQ5       1.17645E9           /* L5/E5a frequency (Hz) */
#define FREQ6       1.27875E9           /* E6/LEX frequency (Hz) */
#define FREQ7       1.20714E9           /* E5b    frequency (Hz) */
#define FREQ8       1.191795E9          /* E5a+b  frequency (Hz) */
#define FREQ1_GLO   1.60200E9           /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO   0.56250E6           /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO   1.24600E9           /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO   0.43750E6           /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO   1.202025E9          /* GLONASS G3 frequency (Hz) */
#define FREQ2_CMP   1.561098E9          /* BeiDou B1 frequency (Hz) */
#define FREQ7_CMP   1.20714E9           /* BeiDou B2 frequency (Hz) */
#define FREQ6_CMP   1.26852E9           /* BeiDou B3 frequency (Hz) */

#define EFACT_GPS   1.0                 /* error factor: GPS */
#define EFACT_GLO   1.5                 /* error factor: GLONASS */
#define EFACT_GAL   1.0                 /* error factor: Galileo */
#define EFACT_QZS   1.0                 /* error factor: QZSS */
#define EFACT_CMP   1.0                 /* error factor: BeiDou */
#define EFACT_SBS   3.0                 /* error factor: SBAS */

#define SYS_NONE    0x00                /* navigation system: none */
#define SYS_GPS     0x01                /* navigation system: GPS */
#define SYS_SBS     0x02                /* navigation system: SBAS */
#define SYS_GLO     0x04                /* navigation system: GLONASS */
#define SYS_GAL     0x08                /* navigation system: Galileo */
#define SYS_QZS     0x10                /* navigation system: QZSS */
#define SYS_CMP     0x20                /* navigation system: BeiDou */
#define SYS_ALL     0xFF                /* navigation system: all */

#define TSYS_GPS    0                   /* time system: GPS time */
#define TSYS_UTC    1                   /* time system: UTC */
#define TSYS_GLO    2                   /* time system: GLONASS time */
#define TSYS_GAL    3                   /* time system: Galileo time */
#define TSYS_QZS    4                   /* time system: QZSS time */
#define TSYS_CMP    5                   /* time system: BeiDou time */
typedef struct
{                                 /* navigation data type */
    int n, nmax;                  /* number of broadcast ephemeris */
    int ng, ngmax;                /* number of glonass ephemeris */
    int ns, nsmax;                /* number of sbas ephemeris */
    int ne, nemax;                /* number of precise ephemeris */
    int nc, ncmax;                /* number of precise clock */
    int na, namax;                /* number of almanac data */
    int nt, ntmax;                /* number of tec grid data */
    int nn, nnmax;                /* number of stec grid data */
    eph_t *eph;                   /* GPS/QZS/GAL ephemeris */
    geph_t *geph;                 /* GLONASS ephemeris */
    seph_t *seph;                 /* SBAS ephemeris */
    peph_t *peph;                 /* precise ephemeris */
    pclk_t *pclk;                 /* precise clock */
    alm_t *alm;                   /* almanac data */
    tec_t *tec;                   /* tec grid data */
    stec_t *stec;                 /* stec grid data */
    erp_t erp;                    /* earth rotation parameters */
    double utc_gps[4];            /* GPS delta-UTC parameters {A0,A1,T,W} */
    double utc_glo[4];            /* GLONASS UTC GPS time parameters */
    double utc_gal[4];            /* Galileo UTC GPS time parameters */
    double utc_qzs[4];            /* QZS UTC GPS time parameters */
    double utc_cmp[4];            /* BeiDou UTC parameters */
    double utc_sbs[4];            /* SBAS UTC parameters */
    double ion_gps[8];            /* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_gal[4];            /* Galileo iono model parameters {ai0,ai1,ai2,0} */
    double ion_qzs[8];            /* QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_cmp[8];            /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    int leaps;                    /* leap seconds (s) */
    double lam[MAXSAT][NFREQ];    /* carrier wave lengths (m) */
    double cbias[MAXSAT][3];      /* code bias (0:p1-p2,1:p1-c1,2:p2-c2) (m) */
    double wlbias[MAXSAT];        /* wide-lane bias (cycle) */
    double glo_cpbias[4];         /* glonass code-phase bias {1C,1P,2C,2P} (m) */
    char glo_fcn[MAXPRNGLO + 1];  /* glonass frequency channel number + 8 */
    pcv_t pcvs[MAXSAT];           /* satellite antenna pcv */
    sbssat_t sbssat;              /* SBAS satellite corrections */
    sbsion_t sbsion[MAXBAND + 1]; /* SBAS ionosphere corrections */
    dgps_t dgps[MAXSAT];          /* DGPS corrections */
    ssr_t ssr[MAXSAT];            /* SSR corrections */
    lexeph_t lexeph[MAXSAT];      /* LEX ephemeris */
    lexion_t lexion;              /* LEX ionosphere correction */
} nav_t;