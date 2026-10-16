#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ------------- grid / physics ------------- */
#define NX          160          /* streamwise  (flow goes +X)   */
#define NY          64           /* vertical    (camera looks along this) */
#define NZ          64           /* spanwise                      */
#define NSTEPS      500
#define SAVE_EVERY  25
#define U0          0.08f        /* inlet speed (keep < 0.15)     */
#define TAU         0.65f        /* relaxation time -> viscosity  */
#define SPH_R       7.0f         /* sphere radius (cells)         */
#define SPH_X       (NX/3.0f)
#define SPH_Y       (NY/2.0f)
#define SPH_Z       (NZ/2.0f)

/* ------------- visualization ------------- */
#define ALPHA_SCALE 0.75f        /* wake opacity  (tune 0.4 - 1.5) */
#define VMAX_SCALE  1.5f         /* color scale: vmax = VMAX_SCALE*U0 */

/* ------------- D3Q19 lattice ------------- */
/*  k:  0 = rest
 *      1..6   = face neighbours  (w=1/18)
 *      7..18  = edge neighbours  (w=1/36)                        */
static const int cx[19] = {0, 1,-1, 0, 0, 0, 0,
                           1,-1, 1,-1, 1,-1, 1,-1,
                           0, 0, 0, 0};
static const int cy[19] = {0, 0, 0, 1,-1, 0, 0,
                           1,-1,-1, 1, 0, 0, 0, 0,
                           1,-1, 1,-1};
static const int cz[19] = {0, 0, 0, 0, 0, 1,-1,
                           0, 0, 0, 0, 1,-1,-1, 1,
                           1,-1,-1, 1};
static const int op[19] = {0, 2, 1, 4, 3, 6, 5,
                           8, 7,10, 9,12,11,14,13,
                           16,15,18,17};
static const float w[19] = {
    1.0f/3.0f,
    1.0f/18, 1.0f/18, 1.0f/18, 1.0f/18, 1.0f/18, 1.0f/18,
    1.0f/36, 1.0f/36, 1.0f/36, 1.0f/36, 1.0f/36, 1.0f/36,
    1.0f/36, 1.0f/36, 1.0f/36, 1.0f/36, 1.0f/36, 1.0f/36
};

/* ------------- state ------------- */
static float *f, *fnew;            /* distributions, AoS: f[i*19 + k] */
static float *rho, *ux, *uy, *uz;  /* macroscopic fields              */
static unsigned char *solid;

static inline int idx(int x, int y, int z) {
    return (z*NY + y)*NX + x;
}

/* Maxwell-Boltzmann equilibrium */
static inline float feq(int k, float r, float u, float v, float t) {
    float cu = cx[k]*u + cy[k]*v + cz[k]*t;
    float u2 = u*u + v*v + t*t;
    return w[k] * r * (1.0f + 3.0f*cu + 4.5f*cu*cu - 1.5f*u2);
}

/* ------------- BMP writer (24-bit, no padding drama) ------------- */
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
        int sy = H - 1 - y;                     /* BMP is bottom-up */
        for (int x = 0; x < W; ++x) {
            row[x*3+0] = img[(sy*W + x)*3 + 2]; /* B */
            row[x*3+1] = img[(sy*W + x)*3 + 1]; /* G */
            row[x*3+2] = img[(sy*W + x)*3 + 0]; /* R */
        }
        fwrite(row, 1, row_bytes, fp);
    }
    free(row);
    fclose(fp);
}

/* classic "jet" CFD colormap */
static void colormapf(float t, float *r, float *g, float *b)
{
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;
    float k = 4.0f * t;
    *r = fminf(fmaxf(fminf(k - 1.5f, -k + 4.5f), 0.0f), 1.0f);
    *g = fminf(fmaxf(fminf(k - 0.5f, -k + 3.5f), 0.0f), 1.0f);
    *b = fminf(fmaxf(fminf(k + 0.5f, -k + 2.5f), 0.0f), 1.0f);
}

int main(void)
{
    const long N = (long)NX * NY * NZ;
    printf("Grid %d x %d x %d = %ld cells (~%.0f MB)\n",
           NX, NY, NZ, N, 2.0*19.0*4.0*N/1e6);

    f     = malloc(sizeof(float) * 19 * N);
    fnew  = malloc(sizeof(float) * 19 * N);
    rho   = malloc(sizeof(float) * N);
    ux    = malloc(sizeof(float) * N);
    uy    = malloc(sizeof(float) * N);
    uz    = malloc(sizeof(float) * N);
    solid = calloc(N, 1);
    if (!f || !fnew || !rho || !ux || !uy || !uz || !solid) {
        fprintf(stderr, "Out of memory\n"); return 1;
    }

    /* --- carve out a sphere --- */
    for (int z = 0; z < NZ; ++z)
    for (int y = 0; y < NY; ++y)
    for (int x = 0; x < NX; ++x) {
        float dx = x - SPH_X, dy = y - SPH_Y, dz = z - SPH_Z;
        if (dx*dx + dy*dy + dz*dz < SPH_R*SPH_R)
            solid[idx(x,y,z)] = 1;
    }

    /* --- initialize with uniform flow + tiny asymmetric seed --- */
    for (int z = 0; z < NZ; ++z)
    for (int y = 0; y < NY; ++y)
    for (int x = 0; x < NX; ++x) {
        int i = idx(x,y,z);
        rho[i] = 1.0f;
        ux[i]  = U0;
        uy[i]  = 1e-3f * sinf(0.11f*x) * cosf(0.13f*z);
        uz[i]  = 1e-3f * cosf(0.11f*x) * sinf(0.13f*y);
        if (solid[i]) { ux[i]=uy[i]=uz[i]=0.0f; }
        for (int k = 0; k < 19; ++k)
            f[i*19 + k] = feq(k, rho[i], ux[i], uy[i], uz[i]);
    }

    unsigned char *img = malloc((size_t)NX * NZ * 3);
    char path[256];
    int frame = 0;

    printf("Simulating %d steps...\n", NSTEPS);

    for (int step = 0; step <= NSTEPS; ++step) {

        /* ---- collision + macroscopic moments (fused) ---- */
        #pragma omp parallel for schedule(static)
        for (int z = 0; z < NZ; ++z)
        for (int y = 0; y < NY; ++y)
        for (int x = 0; x < NX; ++x) {
            int i = idx(x,y,z);
            if (solid[i]) continue;

            float r = 0.0f, mx = 0.0f, my = 0.0f, mz = 0.0f;
            const float *fi = f + i*19;
            for (int k = 0; k < 19; ++k) {
                float fk = fi[k];
                r  += fk;
                mx += fk * cx[k];
                my += fk * cy[k];
                mz += fk * cz[k];
            }
            float invr = 1.0f / r;
            float u = mx * invr, v = my * invr, t = mz * invr;
            rho[i] = r; ux[i] = u; uy[i] = v; uz[i] = t;

            const float inv_tau = 1.0f / TAU;
            float *fn = fnew + i*19;
            for (int k = 0; k < 19; ++k) {
                float fe = feq(k, r, u, v, t);
                fn[k] = fi[k] - (fi[k] - fe) * inv_tau;
            }
        }

        /* ---- streaming (pull) + boundary conditions ---- */
        #pragma omp parallel for schedule(static)
        for (int z = 0; z < NZ; ++z)
        for (int y = 0; y < NY; ++y)
        for (int x = 0; x < NX; ++x) {
            int i = idx(x,y,z);
            if (solid[i]) continue;

            float *fi = f + i*19;
            const float *fn = fnew + i*19;

            for (int k = 0; k < 19; ++k) {
                int xs = x - cx[k];
                int ys = y - cy[k];
                int zs = z - cz[k];

                if (xs < 0) {
                    /* inlet: impose equilibrium at U0 */
                    fi[k] = feq(k, 1.0f, U0, 0.0f, 0.0f);
                } else if (xs >= NX) {
                    /* outlet: zero-gradient (copy own post-collision) */
                    fi[k] = fn[k];
                } else if (ys < 0 || ys >= NY || zs < 0 || zs >= NZ) {
                    /* lateral walls: no-slip bounce-back */
                    fi[k] = fn[op[k]];
                } else {
                    int si = idx(xs,ys,zs);
                    if (solid[si])
                        fi[k] = fn[op[k]];        /* bounce at sphere */
                    else
                        fi[k] = fnew[si*19 + k];  /* advect */
                }
            }
        }

        /* ---- volume render + save ---- */
        if (step % SAVE_EVERY == 0) {
            const float umax = VMAX_SCALE * U0;

            #pragma omp parallel for schedule(static)
            for (int z = 0; z < NZ; ++z) {
                for (int x = 0; x < NX; ++x) {
                    float R = 0.0f, G = 0.0f, B = 0.0f, A = 0.0f;

                    /* march the ray through Y */
                    for (int y = 0; y < NY; ++y) {
                        int i = idx(x,y,z);
                        float r, g, b, a;

                        if (solid[i]) {
                            r = 0.10f; g = 0.10f; b = 0.13f;
                            a = 1.0f;
                        } else {
                            float vm = sqrtf(ux[i]*ux[i]
                                           + uy[i]*uy[i]
                                           + uz[i]*uz[i]);
                            colormapf(vm / umax, &r, &g, &b);
                            /* opacity from deviation from freestream */
                            float dev = fabsf(vm - U0) / U0;
                            a = ALPHA_SCALE * dev;
                            if (a > 1.0f) a = 1.0f;
                        }

                        float contrib = (1.0f - A) * a;
                        R += contrib * r;
                        G += contrib * g;
                        B += contrib * b;
                        A += contrib;
                        if (A > 0.998f) break;   /* early out: opaque */
                    }

                    int p = (z*NX + x);
                    if (R > 1.0f) R = 1.0f;
                    if (G > 1.0f) G = 1.0f;
                    if (B > 1.0f) B = 1.0f;
                    img[p*3+0] = (unsigned char)(255.0f * R);
                    img[p*3+1] = (unsigned char)(255.0f * G);
                    img[p*3+2] = (unsigned char)(255.0f * B);
                }
            }

            snprintf(path, sizeof(path), "frame_%05d.bmp", frame++);
            write_bmp(path, img, NX, NZ);
            if (frame % 10 == 0)
                printf("  frame %d (step %d)\n", frame, step);
        }
    }

    printf("Done. %d frames written.\n", frame);
    printf("Encode:\n"
           "  ffmpeg -framerate 30 -i frame_%%05d.bmp "
           "-pix_fmt yuv420p out.mp4\n");

    free(f); free(fnew); free(rho);
    free(ux); free(uy); free(uz); free(solid); free(img);
    return 0;
}