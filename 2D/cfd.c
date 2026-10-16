#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ---------- grid / physics ---------- */
#define NX            500
#define NY            250
#define NSTEPS        8000
#define SAVE_EVERY    500
#define U0            0.10f
#define TAU           0.60f
#define VMAX_SCALE    1.6f
#define M_PI 3.14159265358979323846

/* ---------- airfoil definition ---------- */
/* NACA 4-digit: MPXX  e.g. 2412 -> M=0.02 P=0.40 T=0.12 */
#define NACA_M           0.02f
#define NACA_P           0.40f
#define NACA_T           0.12f

/* --- resolution knobs --- */
#define AIRFOIL_CHORD    140.0f    /* chord length in grid cells  */
#define AIRFOIL_SAMPLES  300       /* polygon surface resolution  */
#define AIRFOIL_AOA      5.0f      /* angle of attack, degrees    */

/* pivot = 1/4-chord point, in grid coordinates */
#define AIRFOIL_X        175.0f
#define AIRFOIL_Y        (NY/2.0f)

/* ----------------------------------------- */

static const int    cx[9] = {0, 1, 0,-1, 0, 1,-1,-1, 1};
static const int    cy[9] = {0, 0, 1, 0,-1, 1, 1,-1,-1};
static const int    op[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
static const double w [9] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9,
                             1.0/36,1.0/36,1.0/36,1.0/36};

static float  *f, *fnew, *rho, *ux, *uy;
static int    *solid;

/* ---- airfoil polygon ---- */
static float poly_x[2 * AIRFOIL_SAMPLES + 4];
static float poly_y[2 * AIRFOIL_SAMPLES + 4];
static int   poly_n = 0;

/* NACA 4-digit surface points at chord station xc (0..1) */
static void naca_surface(float xc,
                         float *xu, float *yu,
                         float *xl, float *yl)
{
    float yt = 5.0f * NACA_T *
               (0.2969f * sqrtf(xc)
              - 0.1260f * xc
              - 0.3516f * xc * xc
              + 0.2843f * xc * xc * xc
              - 0.1015f * xc * xc * xc * xc);

    float yc, dyc;
    if (xc < NACA_P) {
        yc  =  NACA_M / (NACA_P * NACA_P) * (2.0f*NACA_P*xc - xc*xc);
        dyc =  2.0f * NACA_M / (NACA_P * NACA_P) * (NACA_P - xc);
    } else {
        float q = 1.0f - NACA_P;
        yc  =  NACA_M / (q*q) *
               ((1.0f - 2.0f*NACA_P) + 2.0f*NACA_P*xc - xc*xc);
        dyc =  2.0f * NACA_M / (q*q) * (NACA_P - xc);
    }

    float theta = atanf(dyc);
    *xu = xc - yt * sinf(theta);
    *yu = yc + yt * cosf(theta);
    *xl = xc + yt * sinf(theta);
    *yl = yc - yt * cosf(theta);
}

/* Build the airfoil outline (upper LE->TE, lower TE->LE), rotated and placed */
static void build_airfoil_poly(void)
{
    float alpha = AIRFOIL_AOA * (float)M_PI / 180.0f;
    float ca = cosf(alpha), sa = sinf(alpha);
    const float pivot_lx = 0.25f, pivot_ly = 0.0f;

    poly_n = 0;

    /* upper surface, LE -> TE (cosine-spaced for good leading-edge detail) */
    for (int i = 0; i <= AIRFOIL_SAMPLES; ++i) {
        float beta = (float)M_PI * i / (float)AIRFOIL_SAMPLES;
        float xc   = 0.5f * (1.0f - cosf(beta));
        float xu, yu, xl, yl;
        naca_surface(xc, &xu, &yu, &xl, &yl);

        float dx = (xu - pivot_lx) * AIRFOIL_CHORD;
        float dy = (yu - pivot_ly) * AIRFOIL_CHORD;
        float rx = dx * ca - dy * sa;
        float ry = dx * sa + dy * ca;

        poly_x[poly_n] = AIRFOIL_X + rx;
        poly_y[poly_n] = AIRFOIL_Y + ry;
        poly_n++;
    }

    /* lower surface, TE -> LE, skipping shared TE and LE points */
    for (int i = AIRFOIL_SAMPLES - 1; i >= 1; --i) {
        float beta = (float)M_PI * i / (float)AIRFOIL_SAMPLES;
        float xc   = 0.5f * (1.0f - cosf(beta));
        float xu, yu, xl, yl;
        naca_surface(xc, &xu, &yu, &xl, &yl);

        float dx = (xl - pivot_lx) * AIRFOIL_CHORD;
        float dy = (yl - pivot_ly) * AIRFOIL_CHORD;
        float rx = dx * ca - dy * sa;
        float ry = dx * sa + dy * ca;

        poly_x[poly_n] = AIRFOIL_X + rx;
        poly_y[poly_n] = AIRFOIL_Y + ry;
        poly_n++;
    }
}

/* even-odd point-in-polygon test */
static int point_in_poly(float px, float py)
{
    int inside = 0;
    for (int i = 0, j = poly_n - 1; i < poly_n; j = i++) {
        float yi = poly_y[i], yj = poly_y[j];
        if (((yi > py) != (yj > py)) &&
            (px < (poly_x[j] - poly_x[i]) *
                   (py - yi) / (yj - yi) + poly_x[i]))
            inside = !inside;
    }
    return inside;
}

/* ---------- BMP writer ---------- */
static void write_bmp(const char *path, const unsigned char *img,
                      int W, int H)
{
    int row_bytes = (W * 3 + 3) & ~3;
    int data_size = row_bytes * H;
    int file_size = 54 + data_size;

    unsigned char hdr[54];
    memset(hdr, 0, 54);
    hdr[0]='B'; hdr[1]='M';
    hdr[2]= file_size        &0xFF;
    hdr[3]=(file_size  >>  8)&0xFF;
    hdr[4]=(file_size  >> 16)&0xFF;
    hdr[5]=(file_size  >> 24)&0xFF;
    hdr[10]=54; hdr[14]=40;
    hdr[18]= W       &0xFF; hdr[19]=(W       >> 8)&0xFF;
    hdr[20]=(W>>16)  &0xFF; hdr[21]=(W       >>24)&0xFF;
    hdr[22]= H       &0xFF; hdr[23]=(H       >> 8)&0xFF;
    hdr[24]=(H>>16)  &0xFF; hdr[25]=(H       >>24)&0xFF;
    hdr[26]=1; hdr[28]=24;
    hdr[34]= data_size       &0xFF;
    hdr[35]=(data_size >>  8)&0xFF;
    hdr[36]=(data_size >> 16)&0xFF;
    hdr[37]=(data_size >> 24)&0xFF;

    FILE *fp = fopen(path, "wb");
    if (!fp) { perror(path); return; }
    fwrite(hdr, 1, 54, fp);

    unsigned char *row = malloc(row_bytes);
    for (int y = 0; y < H; ++y) {
        memset(row, 0, row_bytes);
        int sy = H - 1 - y;
        for (int x = 0; x < W; ++x) {
            row[x*3+0] = img[(sy*W + x)*3 + 2];
            row[x*3+1] = img[(sy*W + x)*3 + 1];
            row[x*3+2] = img[(sy*W + x)*3 + 0];
        }
        fwrite(row, 1, row_bytes, fp);
    }
    free(row);
    fclose(fp);
}

/* jet colormap */
static void colormap(double t, unsigned char *r,
                     unsigned char *g, unsigned char *b)
{
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    double k = 4.0 * t;
    double rr = fmin(fmax(fmin(k-1.5, -k+4.5), 0.0), 1.0);
    double gg = fmin(fmax(fmin(k-0.5, -k+3.5), 0.0), 1.0);
    double bb = fmin(fmax(fmin(k+0.5, -k+2.5), 0.0), 1.0);
    *r = (unsigned char)(255.0 * rr);
    *g = (unsigned char)(255.0 * gg);
    *b = (unsigned char)(255.0 * bb);
}

static inline double feq(int k, double r, double u, double v)
{
    double cu = cx[k]*u + cy[k]*v;
    double u2 = u*u + v*v;
    return w[k] * r * (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*u2);
}

int main(void)
{
    const int N = NX * NY;
    f    = malloc(sizeof(float) * 9 * N);
    fnew = malloc(sizeof(float) * 9 * N);
    rho  = malloc(sizeof(float) * N);
    ux   = malloc(sizeof(float) * N);
    uy   = malloc(sizeof(float) * N);
    solid= calloc(N, sizeof(int));

    /* --- geometry --- */
    build_airfoil_poly();
    printf("Airfoil: NACA %.0f%.0f%.0f, chord=%.0f cells, %d polygon points, AOA=%.1f deg\n",
           NACA_M*100.0f, NACA_P*10.0f, NACA_T*100.0f,
           AIRFOIL_CHORD, poly_n, AIRFOIL_AOA);

    for (int y = 0; y < NY; ++y)
        for (int x = 0; x < NX; ++x)
            if (point_in_poly((float)x + 0.5f, (float)y + 0.5f))
                solid[y*NX + x] = 1;

    /* --- initialize: uniform flow + tiny perturbation --- */
    for (int y = 0; y < NY; ++y)
        for (int x = 0; x < NX; ++x) {
            int i = y*NX + x;
            rho[i] = 1.0f;
            ux[i]  = U0;
            uy[i]  = 1e-3f * sinf(0.13f * x) * (y > NY/2 ? 1.0f : -1.0f);
            if (solid[i]) { ux[i] = 0.0f; uy[i] = 0.0f; }
            for (int k = 0; k < 9; ++k)
                f[k*N + i] = (float)feq(k, rho[i], ux[i], uy[i]);
        }

    unsigned char *img = malloc(NX * NY * 3);
    char path[256];
    int frame = 0;

    printf("Simulating %d x %d grid, %d steps\n", NX, NY, NSTEPS);

    for (int step = 0; step <= NSTEPS; ++step) {

        /* 1. macroscopic moments */
        for (int i = 0; i < N; ++i) {
            if (solid[i]) continue;
            float r = 0.0f, mx = 0.0f, my = 0.0f;
            for (int k = 0; k < 9; ++k) {
                float fk = f[k*N + i];
                r  += fk;
                mx += fk * (float)cx[k];
                my += fk * (float)cy[k];
            }
            if (r > 1e-12f) {
                rho[i] = r;
                ux [i] = mx / r;
                uy [i] = my / r;
            }
        }

        /* 2. BGK collision */
        for (int i = 0; i < N; ++i) {
            if (solid[i]) continue;
            float u = ux[i], v = uy[i], r = rho[i];
            for (int k = 0; k < 9; ++k) {
                float fe = (float)feq(k, r, u, v);
                fnew[k*N + i] = f[k*N + i]
                              - (f[k*N + i] - fe) / (float)TAU;
            }
        }

        /* 3. streaming + boundaries */
        for (int y = 0; y < NY; ++y) {
            for (int x = 0; x < NX; ++x) {
                int i = y*NX + x;
                if (solid[i]) continue;
                for (int k = 0; k < 9; ++k) {
                    int xs = x - cx[k];
                    int ys = y - cy[k];
                    if (xs < 0) {
                        f[k*N + i] = (float)feq(k, 1.0, U0, 0.0);
                    } else if (xs >= NX) {
                        f[k*N + i] = fnew[k*N + i];
                    } else if (ys < 0 || ys >= NY) {
                        f[k*N + i] = fnew[op[k]*N + i];
                    } else {
                        int si = ys*NX + xs;
                        if (solid[si])
                            f[k*N + i] = fnew[op[k]*N + i];
                        else
                            f[k*N + i] = fnew[k*N + si];
                    }
                }
            }
        }

        /* 4. dump frame */
        if (step % SAVE_EVERY == 0) {
            double umax = VMAX_SCALE * U0;
            for (int i = 0; i < N; ++i) {
                unsigned char r, g, b;
                if (solid[i]) {
                    r = 25; g = 25; b = 30;
                } else {
                    double vm = sqrt((double)ux[i]*ux[i]
                                   + (double)uy[i]*uy[i]);
                    colormap(vm / umax, &r, &g, &b);
                }
                img[i*3+0] = r;
                img[i*3+1] = g;
                img[i*3+2] = b;
            }
            snprintf(path, sizeof(path), "frame_%05d.bmp", frame++);
            write_bmp(path, img, NX, NY);
            if (frame % 10 == 0)
                printf("  frame %d (step %d)\n", frame, step);
        }
    }

    printf("Done. %d frames written.\n", frame);
    printf("Encode:\n"
           "  ffmpeg -framerate 30 -i frame_%%05d.bmp "
           "-pix_fmt yuv420p out.mp4\n");

    free(f); free(fnew); free(rho); free(ux); free(uy);
    free(solid); free(img);
    return 0;
}