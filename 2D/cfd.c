#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ---------- tunable parameters ---------- */
#define NX           400          /* grid width  */
#define NY           150          /* grid height */
#define NSTEPS       4000         /* total iterations */
#define SAVE_EVERY   50           /* dump a frame every N steps */
#define U0           0.10         /* inlet velocity (keep < 0.2) */
#define TAU          0.60         /* relaxation time -> viscosity */
#define CYL_R        (NY/8.0)     /* obstacle radius */
#define CYL_X        (NX/4.0)     /* obstacle center x */
#define CYL_Y        (NY/2.0)     /* obstacle center y */
#define VMAX_SCALE   1.6          /* color scale: vmax = VMAX_SCALE * U0 */
/* ----------------------------------------- */

/* D2Q9 lattice: directions, weights, opposites */
static const int    cx[9] = {0, 1, 0,-1, 0, 1,-1,-1, 1};
static const int    cy[9] = {0, 0, 1, 0,-1, 1, 1,-1,-1};
static const int    op[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
static const double w [9] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9,
                             1.0/36,1.0/36,1.0/36,1.0/36};

static double *f, *fnew, *rho, *ux, *uy;
static int    *solid;

/* ---------- minimal 24-bit BMP writer (no deps) ---------- */
static void write_bmp(const char *path, const unsigned char *img,
                      int W, int H)
{
    int row_bytes  = (W * 3 + 3) & ~3;
    int data_size  = row_bytes * H;
    int file_size  = 54 + data_size;

    unsigned char hdr[54];
    memset(hdr, 0, 54);
    hdr[0]='B'; hdr[1]='M';
    hdr[2]= file_size        & 0xFF;
    hdr[3]=(file_size  >>  8)& 0xFF;
    hdr[4]=(file_size  >> 16)& 0xFF;
    hdr[5]=(file_size  >> 24)& 0xFF;
    hdr[10]=54;
    hdr[14]=40;
    hdr[18]= W        & 0xFF;  hdr[19]=(W        >> 8)&0xFF;
    hdr[20]=(W>>16)   & 0xFF;  hdr[21]=(W        >>24)&0xFF;
    hdr[22]= H        & 0xFF;  hdr[23]=(H        >> 8)&0xFF;
    hdr[24]=(H>>16)   & 0xFF;  hdr[25]=(H        >>24)&0xFF;
    hdr[26]=1;
    hdr[28]=24;
    hdr[34]= data_size       & 0xFF;
    hdr[35]=(data_size >>  8)& 0xFF;
    hdr[36]=(data_size >> 16)& 0xFF;
    hdr[37]=(data_size >> 24)& 0xFF;

    FILE *fp = fopen(path, "wb");
    if (!fp) { perror(path); return; }
    fwrite(hdr, 1, 54, fp);

    unsigned char *row = malloc(row_bytes);
    for (int y = 0; y < H; ++y) {
        memset(row, 0, row_bytes);
        int sy = H - 1 - y;              /* BMP rows are bottom-up i think */
        for (int x = 0; x < W; ++x) {
            row[x*3+0] = img[(sy*W + x)*3 + 2];   /* B */
            row[x*3+1] = img[(sy*W + x)*3 + 1];   /* G */
            row[x*3+2] = img[(sy*W + x)*3 + 0];   /* R */
        }
        fwrite(row, 1, row_bytes, fp);
    }
    free(row);
    fclose(fp);
}

/* ---------- classic 'jet' CFD colormap ---------- */
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

/* ---------- equilibrium distribution ---------- */
static inline double feq(int k, double r, double u, double v)
{
    double cu = cx[k]*u + cy[k]*v;
    double u2 = u*u + v*v;
    return w[k] * r * (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*u2);
}

int main(void)
{
    const int N = NX * NY;
    f    = malloc(sizeof(double) * 9 * N);
    fnew = malloc(sizeof(double) * 9 * N);
    rho  = malloc(sizeof(double) * N);
    ux   = malloc(sizeof(double) * N);
    uy   = malloc(sizeof(double) * N);
    solid= calloc(N, sizeof(int));

    /* --- geometry: circular obstacle --- */
    for (int y = 0; y < NY; ++y)
        for (int x = 0; x < NX; ++x) {
            double dx = x - CYL_X, dy = y - CYL_Y;
            if (dx*dx + dy*dy < CYL_R*CYL_R)
                solid[y*NX + x] = 1;
        }

    /* --- initialize: uniform flow + tiny asymmetric perturbation
       to trigger vortex shedding faster --- */
    for (int y = 0; y < NY; ++y)
        for (int x = 0; x < NX; ++x) {
            int i = y*NX + x;
            rho[i] = 1.0;
            ux[i]  = U0;
            uy[i]  = 1e-3 * sin(0.15 * x) * (y > NY/2 ? 1.0 : -1.0);
            if (solid[i]) { ux[i] = 0.0; uy[i] = 0.0; }
            for (int k = 0; k < 9; ++k)
                f[k*N + i] = feq(k, rho[i], ux[i], uy[i]);
        }

    unsigned char *img = malloc(NX * NY * 3);
    char path[256];
    int frame = 0;

    printf("Simulating %d x %d grid, %d steps, tau=%.2f, U0=%.2f\n",
           NX, NY, NSTEPS, TAU, U0);

    for (int step = 0; step <= NSTEPS; ++step) {

        /* 1. macroscopic moments */
        for (int i = 0; i < N; ++i) {
            if (solid[i]) continue;
            double r = 0.0, mx = 0.0, my = 0.0;
            for (int k = 0; k < 9; ++k) {
                double fk = f[k*N + i];
                r  += fk;
                mx += fk * cx[k];
                my += fk * cy[k];
            }
            if (r > 1e-12) {
                rho[i] = r;
                ux [i] = mx / r;
                uy [i] = my / r;
            }
        }

        /* 2. collision (BGK) */
        for (int i = 0; i < N; ++i) {
            if (solid[i]) continue;
            double u = ux[i], v = uy[i];
            for (int k = 0; k < 9; ++k) {
                double fe = feq(k, rho[i], u, v);
                fnew[k*N + i] = f[k*N + i]
                              - (f[k*N + i] - fe) / TAU;
            }
        }

        /* 3. streaming + boundary conditions */
        for (int y = 0; y < NY; ++y) {
            for (int x = 0; x < NX; ++x) {
                int i = y*NX + x;
                if (solid[i]) continue;
                for (int k = 0; k < 9; ++k) {
                    int xs = x - cx[k];
                    int ys = y - cy[k];
                    if (xs < 0) {
                        /* inlet: impose equilibrium at U0 */
                        f[k*N + i] = feq(k, 1.0, U0, 0.0);
                    } else if (xs >= NX) {
                        /* outlet: zero-gradient */
                        f[k*N + i] = fnew[k*N + i];
                    } else if (ys < 0 || ys >= NY) {
                        /* top/bottom walls: no-slip bounce-back */
                        f[k*N + i] = fnew[op[k]*N + i];
                    } else {
                        int si = ys*NX + xs;
                        if (solid[si])
                            f[k*N + i] = fnew[op[k]*N + i];  /* bounce */
                        else
                            f[k*N + i] = fnew[k*N + si];     /* advect */
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
                    r = g = b = 35;                 /* obstacle */
                } else {
                    double vm = sqrt(ux[i]*ux[i] + uy[i]*uy[i]);
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
    printf("Encode with:\n"
           "  ffmpeg -framerate 30 -i frame_%%05d.bmp "
           "-pix_fmt yuv420p out.mp4\n");

    free(f); free(fnew); free(rho); free(ux); free(uy);
    free(solid); free(img);
    return 0;
}