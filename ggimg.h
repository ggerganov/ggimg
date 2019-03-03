/*! \file ggimg.h
 *  \brief Poor man's 2d and 3d image operations
 *  \author Georgi Gerganov
 */

#pragma once

#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>

namespace ggimg {

    //
    // Single-threaded operations
    //

    /*! \brief RGB pixel to luminance Rec. 601: Y = 0.2989 R + 0.5870 G + 0.1140 B */
    template <typename T> inline float rgb_to_luma601(T r, T g, T b) { return 0.2989f*((float)(r)) + 0.5870f*((float)(g)) + 0.1140f*((float)(b)); }
    template <typename T> inline float rgb_to_luma601(const T * p) { return rgb_to_luma601(*p, *(p + 1), *(p + 2)); }

    /*! \brief RGB pixel to luminance Rec. 709: Y = 0.2126 R + 0.7152 G + 0.0722 B */
    template <typename T> inline float rgb_to_luma709(T r, T g, T b) { return 0.2126f*((float)(r)) + 0.7152f*((float)(g)) + 0.0722f*((float)(b)); }
    template <typename T> inline float rgb_to_luma709(const T * p) { return rgb_to_luma709(*p, *(p + 1), *(p + 2)); }

    template <typename TSrc, typename TDst> bool rgb_to_luma601_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_luma601_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_luma709_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_luma709_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_gray_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_gray_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_c_2d(int nx, int ny, int c, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_c_3d(int nx, int ny, int nz, int c, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_r_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_r_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_g_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_g_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_b_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool rgb_to_b_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool gray_to_rgb_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool gray_to_rgb_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);

    template <typename TSrc, typename TDst> bool normalize_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool normalize_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool normalize_robust_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);
    template <typename TSrc, typename TDst> bool normalize_robust_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst);

    template <typename T> bool normalize_hist_2d(int nx, int ny, const T * src, T * dst, T dmax, int nw = 0, int * work = nullptr, bool wzero = true);
    template <typename T> bool normalize_hist_3d(int nx, int ny, int nz, const T * src, T * dst, T dmax, int nw = 0, int * work = nullptr, bool wzero = true);

    template <typename T> bool gradient_sobel_2d(bool doX, bool doY, int nx, int ny, const T * src, float * dst);
    template <typename TSrc, typename TDst> bool gradient_sobel_2d(int mode, int nx, int ny, const TSrc * src, TDst dmax, TDst * dst, int nw = 0, float * work = nullptr);

    template <typename T> bool gradient_sobel_3d(bool doX, bool doY, bool doZ, int nx, int ny, int nz, const T * src, float * dst);
    template <typename TSrc, typename TDst> bool gradient_sobel_3d(int mode, int nx, int ny, int nz, const TSrc * src, TDst dmax, TDst * dst, int nw = 0, float * work = nullptr);

    template <typename T> bool convolve_2d(int nx, int ny, const T * src, T * dst, int nk, const float * k, int nw = 0, T * work = nullptr, bool wzero = true);
    template <typename T> bool gaussian_filter_2d(int nx, int ny, const T * src, T * dst, float sigma, int nw = 0, T * work = nullptr, bool wzero = true);
    template <typename T> bool median_filter_2d(int nx, int ny, const T * src, T * dst, int k, int nw = 0, int * work = nullptr, bool wzero = true);

    template <typename T> bool scale_nn_2d(int snx, int sny, const T * src, float sx, float sy, int & dnx, int & dny, std::vector<T> & dst);
    template <typename T> bool scale_nn_isotropic_2d(int snx, int sny, const T * src, float s, int & dnx, int & dny, std::vector<T> & dst);

    template <typename T> bool scale_li_2d(int snx, int sny, const T * src, float sx, float sy, int & dnx, int & dny, std::vector<T> & dst);
    template <typename T> bool scale_li_isotropic_2d(int snx, int sny, const T * src, float s, int & dnx, int & dny, std::vector<T> & dst);

    //
    // Multi-threaded operations (to enable, define GGIMG_MT before including this header)
    //

    template <typename T> bool convolve_3d(int nx, int ny, int nz, const T * src, T * dst, int nk, const float * k, int nthreads = 1, int nw = 0, T * work = nullptr, bool wzero = true);
    template <typename T> bool gaussian_filter_3d(int nx, int ny, int nz, const T * src, T * dst, float sigma, int nthreads = 1, int nw = 0, T * work = nullptr, bool wzero = true);
    template <typename T> bool median_filter_3d(int nx, int ny, int nz, const T * src, T * dst, int k, int nthreads = 1, int nw = 0, int * work = nullptr, bool wzero = true);

    template <typename T> bool scale_nn_3d(int snx, int sny, int snz, const T * src, float sx, float sy, float sz, int & dnx, int & dny, int & dnz, std::vector<T> & dst, int nthreads = 1);
    template <typename T> bool scale_nn_isotropic_3d(int snx, int sny, int snz, const T * src, float s, int & dnx, int & dny, int & dnz, std::vector<T> & dst, int nthreads = 1);

    template <typename T> bool scale_li_3d(int snx, int sny, int snz, const T * src, float sx, float sy, float sz, int & dnx, int & dny, int & dnz, std::vector<T> & dst, int nthreads = 1);
    template <typename T> bool scale_li_isotropic_3d(int snx, int sny, int snz, const T * src, float s, int & dnx, int & dny, int & dnz, std::vector<T> & dst, int nthreads = 1);

}

//
// Single-threaded operations
//

namespace ggimg {

    //
    // Implementation
    //

    /*! \brief RGB image to Rec. 601 luminance image */
    template <typename TSrc, typename TDst>
        bool rgb_to_luma601_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return rgb_to_luma601_3d(nx, ny, 1, src, dmin, dmax, dst);
        }

    /*! \brief RGB image to Rec. 601 luminance image */
    template <typename TSrc, typename TDst>
        bool rgb_to_luma601_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            if (nx < 1) return false;
            if (ny < 1) return false;
            if (nz < 1) return false;
            if (src == nullptr || dst == nullptr) return false;
            if (dmax < dmin) return false;

            int n = nx*ny*nz;

            auto imin = rgb_to_luma601(src + 0);
            auto imax = imin;
            for (int i = 1; i < n; ++i) {
                auto icur = rgb_to_luma601(src + 3*i);
                if (icur > imax) imax = icur;
                if (icur < imin) imin = icur;
            }

            float scale = (imax > imin) ? ((float)(dmax - dmin))/(imax - imin) : 1;
            for (int i = 0; i < n; ++i) {
                dst[i] = dmin + ((float)(rgb_to_luma601(src + 3*i) - imin))*scale;
            }

            return true;
        }

    /*! \brief RGB image to Rec. 709 luminance image */
    template <typename TSrc, typename TDst>
        bool rgb_to_luma709_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return rgb_to_luma709_3d(nx, ny, 1, src, dmin, dmax, dst);
        }

    /*! \brief RGB image to Rec. 709 luminance image */
    template <typename TSrc, typename TDst>
        bool rgb_to_luma709_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            if (nx < 1) return false;
            if (ny < 1) return false;
            if (nz < 1) return false;
            if (src == nullptr || dst == nullptr) return false;
            if (dmax < dmin) return false;

            int n = nx*ny*nz;

            auto imin = rgb_to_luma709(src + 0);
            auto imax = imin;
            for (int i = 1; i < n; ++i) {
                auto icur = rgb_to_luma709(src + 3*i);
                if (icur > imax) imax = icur;
                if (icur < imin) imin = icur;
            }

            float scale = (imax > imin) ? ((float)(dmax - dmin))/(imax - imin) : 1;
            for (int i = 0; i < n; ++i) {
                dst[i] = dmin + ((float)(rgb_to_luma709(src + 3*i) - imin))*scale;
            }

            return true;
        }

    /*! \brief RGB image to grayscale image */
    template <typename TSrc, typename TDst>
        bool rgb_to_gray_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return rgb_to_gray_3d(nx, ny, 1, src, dmin, dmax, dst);
        }

    /*! \brief RGB image to grayscale image */
    template <typename TSrc, typename TDst>
        bool rgb_to_gray_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return rgb_to_luma709_3d(nx, ny, nz, src, dmin, dmax, dst);
        }

    template <typename TSrc, typename TDst>
        bool rgb_to_c_2d(int nx, int ny, int c, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return rgb_to_c_2d(nx, ny, 1, c, src, dmin, dmax, dst);
        }

    template <typename TSrc, typename TDst>
        bool rgb_to_c_3d(int nx, int ny, int nz, int c, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            if (nx < 1) return false;
            if (ny < 1) return false;
            if (nz < 1) return false;
            if (src == nullptr || dst == nullptr) return false;
            if (dmax < dmin) return false;

            int n = nx*ny*nz;

            auto imin = src[c];
            auto imax = imin;
            for (int i = 1; i < n; ++i) {
                auto icur = src[3*i + c];
                if (icur > imax) imax = icur;
                if (icur < imin) imin = icur;
            }

            float scale = (imax > imin) ? ((float)(dmax - dmin))/(imax - imin) : 1;
            for (int i = 0; i < n; ++i) {
                dst[i] = dmin + ((float)(src[3*i + c] - imin))*scale;
            }

            return true;
        }

    template <typename TSrc, typename TDst>
        bool rgb_to_r_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return rgb_to_c_2d(nx, ny, 0, src, dmin, dmax, dst);
        }

    template <typename TSrc, typename TDst>
        bool rgb_to_r_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return rgb_to_c_3d(nx, ny, nz, 0, src, dmin, dmax, dst);
        }

    template <typename TSrc, typename TDst>
        bool rgb_to_g_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return rgb_to_c_2d(nx, ny, 1, src, dmin, dmax, dst);
        }

    template <typename TSrc, typename TDst>
        bool rgb_to_g_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return rgb_to_c_3d(nx, ny, nz, 1, src, dmin, dmax, dst);
        }

    template <typename TSrc, typename TDst>
        bool rgb_to_b_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return rgb_to_c_2d(nx, ny, 2, src, dmin, dmax, dst);
        }

    template <typename TSrc, typename TDst>
        bool rgb_to_b_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return rgb_to_c_3d(nx, ny, nz, 2, src, dmin, dmax, dst);
        }

    template <typename TSrc, typename TDst>
        bool gray_to_rgb_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return gray_to_rgb_3d(nx, ny, 1, src, dmin, dmax, dst);
        }

    template <typename TSrc, typename TDst>
        bool gray_to_rgb_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            if (nx < 1) return false;
            if (ny < 1) return false;
            if (nz < 1) return false;
            if (src == nullptr || dst == nullptr) return false;
            if (dmax < dmin) return false;

            int n = nx*ny*nz;

            auto imin = src[0];
            auto imax = imin;
            for (int i = 1; i < n; ++i) {
                auto icur = src[i];
                if (icur > imax) imax = icur;
                if (icur < imin) imin = icur;
            }

            float scale = (imax > imin) ? ((float)(dmax - dmin))/(imax - imin) : 1;
            for (int i = 0; i < n; ++i) {
                float v = dmin + ((float)(src[i] - imin))*scale;
                dst[3*i + 0] = v;
                dst[3*i + 1] = v;
                dst[3*i + 2] = v;
            }

            return true;
        }

    /*! \brief Simple intensity scale.
     *
     * Works in-place (\a src == \a dst)
     * Destination intensities in range [dmin, dmax]
     */
    template <typename TSrc, typename TDst>
        bool normalize_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return normalize_3d(nx, ny, 1, src, dmin, dmax, dst);
        }

    /*! \brief Simple intensity scale.
     *
     * Works in-place (\a src == \a dst)
     * Destination intensities in range [dmin, dmax]
     */
    template <typename TSrc, typename TDst>
        bool normalize_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            if (nx < 1) return false;
            if (ny < 1) return false;
            if (nz < 1) return false;
            if (src == nullptr || dst == nullptr) return false;
            if (dmax < dmin) return false;

            int n = nx*ny*nz;

            auto imin = src[0];
            auto imax = src[0];
            for (int i = 1; i < n; ++i) {
                if (src[i] > imax) imax = src[i];
                if (src[i] < imin) imin = src[i];
            }

            float scale = (imax > imin) ? ((float)(dmax - dmin))/(imax - imin) : 1;
            for (int i = 0; i < n; ++i) {
                dst[i] = dmin + ((float)(src[i] - imin))*scale;
            }

            return true;
        }

    /*! \brief Intensity scale which is robust to small number of huge intensity spikes.
     *
     * Works in-place (\a src == \a dst)
     * Destination intensities in range [dmin, dmax]
     */
    template <typename TSrc, typename TDst>
        bool normalize_robust_2d(int nx, int ny, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            return normalize_robust_3d(nx, ny, 1, src, dmin, dmax, dst);
        }

    /*! \brief Intensity scale which is robust to small number of huge intensity spikes.
     *
     * Works in-place (\a src == \a dst)
     * Destination intensities in range [dmin, dmax]
     */
    template <typename TSrc, typename TDst>
        bool normalize_robust_3d(int nx, int ny, int nz, const TSrc * src, TDst dmin, TDst dmax, TDst * dst) {
            if (nx < 1) return false;
            if (ny < 1) return false;
            if (nz < 1) return false;
            if (src == nullptr || dst == nullptr) return false;
            if (dmax < dmin) return false;

            int n = nx*ny*nz;

            std::vector<TSrc> vals(n);
            for (int i = 0; i < n; ++i) vals[i] = src[i];
            std::sort(vals.begin(), vals.end());

            auto imin = vals[0.01*n];
            auto imax = vals[0.99*n];

            float scale = (imax > imin) ? ((float)(dmax - dmin))/(imax - imin) : 1;
            for (int i = 0; i < n; ++i) {
                auto v = dmin + ((float)(src[i] - imin))*scale;
                if (v < dmin) v = dmin;
                if (v >= dmax) v = dmax;
                dst[i] = v;
            }

            return true;
        }

    /*! \brief Apply Sobel operator to image
     *
     * Used to detect edges.
     * Cannot work in-place (\a src != \a dst)
     * Optional \a work buffer must be with size \a nw >= nx*ny
     *
     * Ref: https://en.wikipedia.org/wiki/Sobel_operator
     */
    template <typename T>
        bool gradient_sobel_2d(bool doX, bool doY, int nx, int ny, const T * src, float * dst) {
            if (nx < 3) return false;
            if (ny < 3) return false;
            if (src == nullptr || dst == nullptr) return false;

            int nx1 = nx - 1;
            int ny1 = ny - 1;

            int nxy = nx*ny;

            {
                for (int iy = 1; iy < ny1; ++iy) {
                    for (int ix = 1; ix < nx1; ++ix) {
                        int id = iy*nx + ix;
                        float g1 = 0;
                        if (doX) {
                            g1 += 1.0f*src[id - 1 - nx      ];
                            g1 += 2.0f*src[id - 1           ];
                            g1 += 1.0f*src[id - 1 + nx      ];
                            g1 -= 1.0f*src[id + 1 - nx      ];
                            g1 -= 2.0f*src[id + 1           ];
                            g1 -= 1.0f*src[id + 1 + nx      ];
                        }

                        float g2 = 0;
                        if (doY) {
                            g2 += 1.0f*src[id - 1 - nx      ];
                            g2 += 2.0f*src[id     - nx      ];
                            g2 += 1.0f*src[id + 1 - nx      ];
                            g2 -= 1.0f*src[id - 1 + nx      ];
                            g2 -= 2.0f*src[id     + nx      ];
                            g2 -= 1.0f*src[id + 1 + nx      ];
                        }

                        dst[id] = sqrt(g1*g1 + g2*g2);
                    }
                }
            }

            // x == 0, x == nx - 1
            for (int iy = 0; iy < ny; ++iy) {
                dst[iy*nx] = dst[iy*nx + 1];
                dst[iy*nx + nx - 1] = dst[iy*nx + nx - 2];
            }

            // y == 0, y == ny - 1
            for (int ix = 0; ix < nx; ++ix) {
                dst[ix] = dst[ix + nx];
                dst[(ny - 1)*nx + ix] = dst[(ny - 2)*nx + ix];
            }

            return true;
        }

    /*! \brief Apply Sobel operator to image
     *
     * Used to detect edges.
     * Can work in-place (\a src == \a dst)
     * Optional \a work buffer must be with size \a nw >= nx*ny
     *
     * At the end, scales image intensities in the range [0, 255]
     *
     * \param mode Do X direction if mode & 1
     *             Do Y direction if mode & 2
     *             Do both directions if mode == 0
     *
     * Ref: https://en.wikipedia.org/wiki/Sobel_operator
     */
    template <typename TSrc, typename TDst>
        bool gradient_sobel_2d(int mode, int nx, int ny, const TSrc * src, TDst dmax, TDst * dst, int nw, float * work) {
            if (nx < 3) return false;
            if (ny < 3) return false;
            if (src == nullptr || dst == nullptr) return false;

            bool wdelete = false;
            if (work) {
                if (nw < nx*ny) return false;
            } else {
                wdelete = true;
                work = new float[nx*ny];
            }

            bool res = false;
            if (mode == 0) {
                res = gradient_sobel_2d(true, true, nx, ny, src, work);
            } else {
                bool doX = mode & 1;
                bool doY = mode & 2;
                res = gradient_sobel_2d(doX, doY, nx, ny, src, work);
            }
            if (res == false) {
                if (wdelete) delete [] work;
                return false;
            }

            if (dmax != 0) {
                if (normalize_2d<float, TDst>(nx, ny, work, 0, dmax, dst) == false) {
                    if (wdelete) delete [] work;
                    return false;
                }
            } else {
                for (int i = 0; i < nx*ny; ++i) {
                    dst[i] = work[i];
                }
            }

            if (wdelete) delete [] work;
            return true;
        }

    /*! \brief Apply Sobel operator to image
     *
     * Used to detect edges.
     * Cannot work in-place (\a src != \a dst)
     * Optional \a work buffer must be with size \a nw >= nx*ny*nz
     *
     * Ref: https://en.wikipedia.org/wiki/Sobel_operator
     */
    template <typename T>
        bool gradient_sobel_3d(bool doX, bool doY, bool doZ, int nx, int ny, int nz, const T * src, float * dst) {
            if (nx < 3) return false;
            if (ny < 3) return false;
            if (nz < 3) return false;
            if (src == nullptr || dst == nullptr) return false;

            int nx1 = nx - 1;
            int ny1 = ny - 1;
            int nz1 = nz - 1;

            int nxy = nx*ny;

            {
                for (int iz = 1; iz < nz1; ++iz) {
                    for (int iy = 1; iy < ny1; ++iy) {
                        for (int ix = 1; ix < nx1; ++ix) {
                            int id = iz*ny*nx + iy*nx + ix;
                            float g1 = 0;
                            if (doX) {
                                g1 += 1.0f*src[id - 1 - nx - nxy];
                                g1 += 2.0f*src[id - 1 - nx      ];
                                g1 += 1.0f*src[id - 1 - nx + nxy];
                                g1 += 2.0f*src[id - 1      - nxy];
                                g1 += 4.0f*src[id - 1           ];
                                g1 += 2.0f*src[id - 1      + nxy];
                                g1 += 1.0f*src[id - 1 + nx - nxy];
                                g1 += 2.0f*src[id - 1 + nx      ];
                                g1 += 1.0f*src[id - 1 + nx + nxy];
                                g1 -= 1.0f*src[id + 1 - nx - nxy];
                                g1 -= 2.0f*src[id + 1 - nx      ];
                                g1 -= 1.0f*src[id + 1 - nx + nxy];
                                g1 -= 2.0f*src[id + 1      - nxy];
                                g1 -= 4.0f*src[id + 1           ];
                                g1 -= 2.0f*src[id + 1      + nxy];
                                g1 -= 1.0f*src[id + 1 + nx - nxy];
                                g1 -= 2.0f*src[id + 1 + nx      ];
                                g1 -= 1.0f*src[id + 1 + nx + nxy];
                            }

                            float g2 = 0;
                            if (doY) {
                                g2 += 1.0f*src[id - 1 - nx - nxy];
                                g2 += 2.0f*src[id - 1 - nx      ];
                                g2 += 1.0f*src[id - 1 - nx + nxy];
                                g2 += 2.0f*src[id     - nx - nxy];
                                g2 += 4.0f*src[id     - nx      ];
                                g2 += 2.0f*src[id     - nx + nxy];
                                g2 += 1.0f*src[id + 1 - nx - nxy];
                                g2 += 2.0f*src[id + 1 - nx      ];
                                g2 += 1.0f*src[id + 1 - nx + nxy];
                                g2 -= 1.0f*src[id - 1 + nx - nxy];
                                g2 -= 2.0f*src[id - 1 + nx      ];
                                g2 -= 1.0f*src[id - 1 + nx + nxy];
                                g2 -= 2.0f*src[id     + nx - nxy];
                                g2 -= 4.0f*src[id     + nx      ];
                                g2 -= 2.0f*src[id     + nx + nxy];
                                g2 -= 1.0f*src[id + 1 + nx - nxy];
                                g2 -= 2.0f*src[id + 1 + nx      ];
                                g2 -= 1.0f*src[id + 1 + nx + nxy];
                            }

                            float g3 = 0;
                            if (doZ) {
                                g3 += 1.0f*src[id - 1 - nx - nxy];
                                g3 += 2.0f*src[id - 1      - nxy];
                                g3 += 1.0f*src[id - 1 + nx - nxy];
                                g3 += 2.0f*src[id     - nx - nxy];
                                g3 += 4.0f*src[id          - nxy];
                                g3 += 2.0f*src[id     + nx - nxy];
                                g3 += 1.0f*src[id + 1 - nx - nxy];
                                g3 += 2.0f*src[id + 1      - nxy];
                                g3 += 1.0f*src[id + 1 + nx - nxy];
                                g3 -= 1.0f*src[id - 1 - nx + nxy];
                                g3 -= 2.0f*src[id - 1      + nxy];
                                g3 -= 1.0f*src[id - 1 + nx + nxy];
                                g3 -= 2.0f*src[id     - nx + nxy];
                                g3 -= 4.0f*src[id          + nxy];
                                g3 -= 2.0f*src[id     + nx + nxy];
                                g3 -= 1.0f*src[id + 1 - nx + nxy];
                                g3 -= 2.0f*src[id + 1      + nxy];
                                g3 -= 1.0f*src[id + 1 + nx + nxy];
                            }

                            dst[id] = sqrt(g1*g1 + g2*g2 + g3*g3);
                        }
                    }
                }
            }

            // x == 0, x == nx - 1
            for (int iz = 0; iz < nz; ++iz) {
                for (int iy = 0; iy < ny; ++iy) {
                    dst[iz*nxy + iy*nx] = dst[iz*nxy + iy*nx + 1];
                    dst[iz*nxy + iy*nx + nx - 1] = dst[iz*nxy + iy*nx + nx - 2];
                }
            }

            // y == 0, y == ny - 1
            for (int iz = 0; iz < nz; ++iz) {
                for (int ix = 0; ix < nx; ++ix) {
                    dst[iz*nxy + ix] = dst[iz*nxy + ix + nx];
                    dst[iz*nxy + (ny - 1)*nx + ix] = dst[iz*nxy + (ny - 2)*nx + ix];
                }
            }

            // z == 0, z == nz - 1
            for (int iy = 0; iy < ny; ++iy) {
                for (int ix = 0; ix < nx; ++ix) {
                    dst[iy*nx + ix] = dst[nxy + iy*nx + ix];
                    dst[(nz - 1)*nxy + iy*nx + ix] = dst[(nz - 2)*nxy + iy*nx + ix];
                }
            }

            return true;
        }

    /*! \brief Apply Sobel operator to image
     *
     * Used to detect edges.
     * Can work in-place (\a src == \a dst)
     * Optional \a work buffer must be with size \a nw >= nx*ny*nz
     *
     * At the end, scales image intensities in the range [0, 255]
     *
     * \param mode Do X direction if mode & 1
     *             Do Y direction if mode & 2
     *             Do Z direction if mode & 4
     *             Do all directions if mode == 0
     *
     * Ref: https://en.wikipedia.org/wiki/Sobel_operator
     */
    template <typename TSrc, typename TDst>
        bool gradient_sobel_3d(int mode, int nx, int ny, int nz, const TSrc * src, TDst dmax, TDst * dst, int nw, float * work) {
            if (nx < 3) return false;
            if (ny < 3) return false;
            if (nz < 3) return false;
            if (src == nullptr || dst == nullptr) return false;

            bool wdelete = false;
            if (work) {
                if (nw < nx*ny*nz) return false;
            } else {
                wdelete = true;
                work = new float[nx*ny*nz];
            }

            bool res = false;
            if (mode == 0) {
                res = gradient_sobel_3d(true, true, true, nx, ny, nz, src, work);
            } else {
                bool doX = mode & 1;
                bool doY = mode & 2;
                bool doZ = mode & 4;
                res = gradient_sobel_3d(doX, doY, doZ, nx, ny, nz, src, work);
            }
            if (res == false) {
                if (wdelete) delete [] work;
                return false;
            }

            if (dmax != 0) {
                if (normalize_3d<float, TDst>(nx, ny, nz, work, 0, dmax, dst) == false) {
                    if (wdelete) delete [] work;
                    return false;
                }
            } else {
                for (int i = 0; i < nx*ny*nz; ++i) {
                    dst[i] = work[i];
                }
            }

            if (wdelete) delete [] work;
            return true;
        }

    /*! \brief Histogram equalization
     *
     * Can work in-place (\a src == \a dst)
     * If \a work buffer is provided, it's size \a nw must at least (imax - imin + 1), where
     * imax and imin are the maximum and minimum intensities in the image.
     * On success, the \a dst image has intensities in the range [0, dmax]
     *
     * Ref: https://en.wikipedia.org/wiki/Histogram_equalization
     */
    template <typename T>
        bool normalize_hist_2d(int nx, int ny, const T * src, T * dst, T dmax, int nw, int * work, bool wzero) {
            return normalize_hist_3d(nx, ny, 1, src, dst, dmax, nw, work, wzero);
        }

    /*! \brief Histogram equalization
     *
     * Can work in-place (\a src == \a dst)
     * If \a work buffer is provided, it's size \a nw must at least (imax - imin + 1), where
     * imax and imin are the maximum and minimum intensities in the image.
     * On success, the \a dst image has intensities in the range [0, dmax]
     *
     * Ref: https://en.wikipedia.org/wiki/Histogram_equalization
     */
    template <typename T>
        bool normalize_hist_3d(int nx, int ny, int nz, const T * src, T * dst, T dmax, int nw, int * work, bool wzero) {
            if (nx <= 0 || ny <= 0 || nz <= 0) return false;
            if (src == nullptr || dst == nullptr) return false;

            int n = nx*ny*nz;
            if (n <= 1) return false;

            auto imin = src[0];
            auto imax = src[0];
            for (int i = 1; i < n; ++i) {
                if (src[i] < imin) imin = src[i];
                if (src[i] > imax) imax = src[i];
            }

            int nbin = imax - imin + 1;

            bool wdelete = false;
            if (work) {
                if (nw < nbin) return false;
            } else {
                wdelete = true;
                work = new int[nbin];
                wzero = true;
            }

            if (wzero) {
                for (int i = 0; i < nbin; ++i) work[i] = 0;
            }

            for (int i = 0; i < n; ++i) {
                ++work[src[i] - imin];
            }

            int sum = 0;
            for (int i = 0; i < nbin; ++i) {
                work[i] += sum;
                sum = work[i];
            }
            for (int i = 0; i < n; ++i) {
                dst[i] = (T)((((float)(work[src[i] - imin] - work[0]))/(n - 1))*dmax);
            }

            if (wdelete) delete [] work;

            return true;
        }

    template <typename T>
        bool convolve_2d(int nx, int ny, const T * src, T * dst, int nk, const float * k, int nw, T * work, bool wzero) {
            if (nx <= 0 || ny <= 0) return false;
            if (src == nullptr || dst == nullptr) return false;
            if (k == nullptr) return false;
            if (nk < 3) return false;
            if (nk % 2 == 0) return false;

            int n = nx*ny;
            if (n <= 1) return false;

            float inorm = 0.0;
            for (int i = 0; i < nk; ++i) {
                inorm += k[i];
            }
            if (std::fabs(inorm) < 1e-3) return false;
            inorm = 1.0/inorm;

            int kw = nk/2;
            int nx1 = nx + 2*kw;
            int ny1 = ny + 2*kw;
            int n1 = nx1*ny1;

            bool wdelete = false;
            if (work == nullptr) {
                wdelete = true;
                work = new T[2*n1];
                wzero = true;
            } else if (nw < 2*n1) {
                return false;
            }

            T * tmp0 = work;
            T * tmp1 = work + n1;

            if (wzero) {
                for (int iy = 0; iy < ny1; ++iy) {
                    for (int ix = 0; ix < kw; ++ix) {
                        tmp0[iy*nx1 + ix] = 0;
                        tmp0[iy*nx1 + (nx1 - ix - 1)] = 0;
                        tmp1[iy*nx1 + ix] = 0;
                        tmp1[iy*nx1 + (nx1 - ix - 1)] = 0;
                    }
                }
            }

            for (int iy = 0; iy < ny; ++iy) {
                for (int ix = 0; ix < nx; ++ix) {
                    tmp0[(iy + kw)*nx1 + (ix + kw)] = src[iy*nx + ix];
                }
            }

            // x
            {
                for (int iy = kw; iy < ny + kw; ++iy) {
                    for (int ix = kw; ix < nx + kw; ++ix) {
                        float v = 0.0;
                        for (int ik = -kw; ik <= kw; ++ik) {
                            v += tmp0[iy*nx1 + (ix + ik)]*k[ik + kw];
                        }
                        tmp1[iy*nx1 + ix] = v*inorm;
                    }
                }
            }

            // y
            {
                for (int ix = kw; ix < nx + kw; ++ix) {
                    for (int iy = kw; iy < ny + kw; ++iy) {
                        float v = 0.0;
                        for (int ik = -kw; ik <= kw; ++ik) {
                            v += tmp1[(iy + ik)*nx1 + ix]*k[ik + kw];
                        }
                        dst[(iy - kw)*nx + (ix - kw)] = v*inorm;
                    }
                }
            }

            if (wdelete) {
                delete [] work;
            }

            return true;
        }

    template <typename T>
        bool gaussian_filter_2d(int nx, int ny, const T * src, T * dst, float sigma, int nw, T * work, bool wzero) {
            if (sigma <= 0.0) {
                for (int i = 0; i < nx*ny; ++i) dst[i] = src[i];
                return true;
            }

            int nk = std::max(1.0, std::sqrt(-2.0*std::log(0.10))*sigma);

            std::vector<float> k(2*nk + 1);
            for (int i = 1; i <= nk; ++i) {
                float x = std::exp(-(i*i)/(2.0*sigma*sigma));
                k[nk - i] = k[nk + i] = x;
            }
            k[nk] = 1;

            return convolve_2d(nx, ny, src, dst, k.size(), k.data(), nw, work, wzero);
        }

    template <typename T>
        bool scale_li_2d(int snx, int sny, const T * src, float sx, float sy, int & dnx, int & dny, std::vector<T> & dst) {
            if (snx <= 0) return false;
            if (sny <= 0) return false;
            if (src == nullptr) return false;
            if (sx <= 0.0f) return false;
            if (sy <= 0.0f) return false;

            if (src == dst.data()) return false;

            dnx = std::round(snx*sx);
            dny = std::round(sny*sy);

            if (dnx <= 0) return false;
            if (dny <= 0) return false;

            dst.resize(dnx*dny);

            {
                float cx = ((float)(snx))/dnx;
                float cy = ((float)(sny))/dny;

                for (int iy = 0; iy < dny; ++iy) {
                    float fy = ((float)(iy) + 0.5f)*cy;

                    int iy1 = (int)(fy + 0.5f);
                    int iy0 = iy1 - 1;
                    if (iy0 < 0) ++iy0;
                    if (iy1 >= sny) --iy1;
                    fy -= (0.5f + (float)(iy0));
                    for (int ix = 0; ix < dnx; ++ix) {
                        float fx = ((float)(ix) + 0.5f)*cx;

                        int ix1 = (int)(fx + 0.5f);
                        int ix0 = ix1 - 1;
                        if (ix0 < 0) ++ix0;
                        if (ix1 >= snx) --ix1;
                        fx -= (0.5f + (float)(ix0));

                        auto v00 = src[iy0*snx + ix0];
                        auto v01 = src[iy0*snx + ix1];
                        auto v10 = src[iy1*snx + ix0];
                        auto v11 = src[iy1*snx + ix1];

                        auto v0 = v00 + fy*(v10 - v00);
                        auto v1 = v01 + fy*(v11 - v01);

                        dst[iy*dnx + ix] = v0 + fx*(v1 - v0);
                    }
                }
            }

            return true;
        }

    template <typename T>
        bool scale_li_isotropic_2d(int snx, int sny, const T * src, float s, int & dnx, int & dny, std::vector<T> & dst) {
            return scale_li_2d(snx, sny, src, s, s, dnx, dny, dst);
        }

    template <typename T>
        bool scale_nn_2d(int snx, int sny, const T * src, float sx, float sy, int & dnx, int & dny, std::vector<T> & dst) {
            if (snx <= 0) return false;
            if (sny <= 0) return false;
            if (src == nullptr) return false;
            if (sx <= 0.0f) return false;
            if (sy <= 0.0f) return false;

            if (src == dst.data()) return false;

            dnx = std::round(snx*sx);
            dny = std::round(sny*sy);

            if (dnx <= 0) return false;
            if (dny <= 0) return false;

            dst.resize(dnx*dny);

            {
                float cx = ((float)(snx))/dnx;
                float cy = ((float)(sny))/dny;

                for (int iy = 0; iy < dny; ++iy) {
                    float fy = ((float)(iy) + 0.5f)*cy;
                    int iy0 = std::round(fy - 0.5f);
                    fy -= (0.5f + (float)(iy0));
                    for (int ix = 0; ix < dnx; ++ix) {
                        float fx = ((float)(ix) + 0.5f)*cx;
                        int ix0 = std::round(fx - 0.5f);
                        dst[iy*dnx + ix] = src[iy0*snx + ix0];
                    }
                }
            }

            return true;
        }

    template <typename T>
        bool scale_nn_isotropic_2d(int snx, int sny, const T * src, float s, int & dnx, int & dny, std::vector<T> & dst) {
            return scale_nn_2d(snx, sny, src, s, s, dnx, dny, dst);
        }

    template <>
        bool median_filter_2d<uint8_t>(int nx, int ny, const uint8_t * src, uint8_t * dst, int k, int nw, int * work, bool wzero) {
            if (nx <= 0) return false;
            if (ny <= 0) return false;
            if (src == nullptr) return false;
            if (dst == nullptr) return false;
            if (src == dst) return false;
            if (k < 0) return false;
            if (2*k >= nx) return false;
            if (2*k >= ny) return false;

            int wsize = 256*nx;

            bool wdelete = false;
            if (work) {
                if (nw < wsize) return false;
            } else {
                wdelete = true;
                work = new int[wsize];
                wzero = true;
            }

            if (wzero) {
                for (int i = 0; i < wsize; ++i) work[i] = 0;
            }

            for (int x = 0; x < nx; ++x) {
                for (int y = 0; y <= k; ++y) {
                    ++work[256*x + src[y*nx + x]];
                }
            }

            int j = 0;
            int nker = 0;
            int hker[256];

            for (int y = 0; y < k; ++y) {
                nker = 0;
                std::fill(hker, hker + 256, 0);
                for (int i = 0; i <= k; ++i) {
                    for (j = 0; j < 256; ++j) {
                        hker[j] += work[256*i + j];
                        nker += work[256*i + j];
                    }
                }

                for (int x = 0; x < k; ++x) {
                    int ncur = 0;
                    for (j = 0; j < 255 && ncur < nker/2; ++j) {
                        ncur += hker[j];
                    }
                    dst[y*nx + x] = j;

                    for (j = 0; j < 256; ++j) {
                        hker[j] += work[256*(x + k + 1) + j];
                        nker += work[256*(x + k + 1) + j];
                    }
                }

                for (int x = k; x < nx - k - 1; ++x) {
                    int ncur = 0;
                    for (j = 0; j < 255 && ncur < nker/2; ++j) {
                        ncur += hker[j];
                    }
                    dst[y*nx + x] = j;

                    for (j = 0; j < 256; ++j) {
                        hker[j] -= work[256*(x - k) + j];
                        nker -= work[256*(x - k) + j];
                    }

                    for (j = 0; j < 256; ++j) {
                        hker[j] += work[256*(x + k + 1) + j];
                        nker += work[256*(x + k + 1) + j];
                    }
                }

                for (int x = nx - k - 1; x < nx; ++x) {
                    int ncur = 0;
                    for (j = 0; j < 255 && ncur < nker/2; ++j) {
                        ncur += hker[j];
                    }
                    dst[y*nx + x] = j;

                    for (j = 0; j < 256; ++j) {
                        hker[j] -= work[256*(x - k) + j];
                        nker -= work[256*(x - k) + j];
                    }
                }

                for (int x = 0; x < nx; ++x) {
                    ++work[256*x + src[(y + k + 1)*nx + x]];
                }
            }

            for (int y = k; y < ny - k; ++y) {
                nker = 0;
                std::fill(hker, hker + 256, 0);
                for (int i = 0; i <= k; ++i) {
                    for (j = 0; j < 256; ++j) {
                        hker[j] += work[256*i + j];
                        nker += work[256*i + j];
                    }
                }

                for (int x = 0; x < k; ++x) {
                    int ncur = 0;
                    for (j = 0; j < 255 && ncur < nker/2; ++j) {
                        ncur += hker[j];
                    }
                    dst[y*nx + x] = j;

                    for (j = 0; j < 256; ++j) {
                        hker[j] += work[256*(x + k + 1) + j];
                        nker += work[256*(x + k + 1) + j];
                    }
                }

                for (int x = k; x < nx - k - 1; ++x) {
                    int ncur = 0;
                    for (j = 0; j < 255 && ncur < nker/2; ++j) {
                        ncur += hker[j];
                    }
                    dst[y*nx + x] = j;

                    for (j = 0; j < 256; ++j) {
                        hker[j] -= work[256*(x - k) + j];
                        nker -= work[256*(x - k) + j];
                    }

                    for (j = 0; j < 256; ++j) {
                        hker[j] += work[256*(x + k + 1) + j];
                        nker += work[256*(x + k + 1) + j];
                    }
                }

                for (int x = nx - k - 1; x < nx; ++x) {
                    int ncur = 0;
                    for (j = 0; j < 255 && ncur < nker/2; ++j) {
                        ncur += hker[j];
                    }
                    dst[y*nx + x] = j;

                    for (j = 0; j < 256; ++j) {
                        hker[j] -= work[256*(x - k) + j];
                        nker -= work[256*(x - k) + j];
                    }
                }

                for (int x = 0; x < nx; ++x) {
                    --work[256*x + src[(y - k)*nx + x]];
                    ++work[256*x + src[(y + k + 1)*nx + x]];
                }
            }

            for (int y = ny - k; y < ny; ++y) {
                nker = 0;
                std::fill(hker, hker + 256, 0);
                for (int i = 0; i <= k; ++i) {
                    for (j = 0; j < 256; ++j) {
                        hker[j] += work[256*i + j];
                        nker += work[256*i + j];
                    }
                }

                for (int x = 0; x < k; ++x) {
                    int ncur = 0;
                    for (j = 0; j < 255 && ncur < nker/2; ++j) {
                        ncur += hker[j];
                    }
                    dst[y*nx + x] = j;

                    for (j = 0; j < 256; ++j) {
                        hker[j] += work[256*(x + k + 1) + j];
                        nker += work[256*(x + k + 1) + j];
                    }
                }

                for (int x = k; x < nx - k - 1; ++x) {
                    int ncur = 0;
                    for (j = 0; j < 255 && ncur < nker/2; ++j) {
                        ncur += hker[j];
                    }
                    dst[y*nx + x] = j;

                    for (j = 0; j < 256; ++j) {
                        hker[j] -= work[256*(x - k) + j];
                        nker -= work[256*(x - k) + j];
                    }

                    for (j = 0; j < 256; ++j) {
                        hker[j] += work[256*(x + k + 1) + j];
                        nker += work[256*(x + k + 1) + j];
                    }
                }

                for (int x = nx - k - 1; x < nx; ++x) {
                    int ncur = 0;
                    for (j = 0; j < 255 && ncur < nker/2; ++j) {
                        ncur += hker[j];
                    }
                    dst[y*nx + x] = j;

                    for (j = 0; j < 256; ++j) {
                        hker[j] -= work[256*(x - k) + j];
                        nker -= work[256*(x - k) + j];
                    }
                }

                for (int x = 0; x < nx; ++x) {
                    --work[256*x + src[(y - k)*nx + x]];
                }
            }

            if (wdelete) delete [] work;

            return true;
        }

}

#ifdef GGIMG_MT

//
// Multi-threaded operations
//

#include <atomic>
#include <thread>

namespace ggimg {
    template <typename T>
        bool convolve_3d(int nx, int ny, int nz, const T * src, T * dst, int nk, const float * k, int nthreads, int nw, T * work, bool wzero) {
            if (nx <= 0 || ny <= 0 || nz <= 0) return false;
            if (src == nullptr || dst == nullptr) return false;
            if (k == nullptr) return false;
            if (nk < 3) return false;
            if (nk % 2 == 0) return false;
            if (nthreads < 1) return false;

            int n = nx*ny*nz;
            if (n <= 1) return false;

            float inorm = 0.0;
            for (int i = 0; i < nk; ++i) {
                inorm += k[i];
            }
            if (std::fabs(inorm) < 1e-3) return false;
            inorm = 1.0/inorm;

            int kw = nk/2;
            int nx1 = nx + 2*kw;
            int ny1 = ny + 2*kw;
            int nz1 = nz + 2*kw;
            int n1 = nx1*ny1*nz1;

            bool wdelete = false;
            if (work == nullptr) {
                wdelete = true;
                work = new T[2*n1];
                wzero = true;
            } else if (nw < 2*n1) {
                return false;
            }

            T * tmp0 = work;
            T * tmp1 = work + n1;

            if (wzero) {
                for (int iz = 0; iz < nz1; ++iz) {
                    for (int iy = 0; iy < ny1; ++iy) {
                        for (int ix = 0; ix < kw; ++ix) {
                            tmp0[iz*ny1*nx1 + iy*nx1 + ix] = 0;
                            tmp0[iz*ny1*nx1 + iy*nx1 + (nx1 - ix - 1)] = 0;
                            tmp1[iz*ny1*nx1 + iy*nx1 + ix] = 0;
                            tmp1[iz*ny1*nx1 + iy*nx1 + (nx1 - ix - 1)] = 0;
                        }
                    }
                }

                for (int iz = 0; iz < nz1; ++iz) {
                    for (int iy = 0; iy < kw; ++iy) {
                        for (int ix = 0; ix < nx1; ++ix) {
                            tmp0[iz*ny1*nx1 + iy*nx1 + ix] = 0;
                            tmp0[iz*ny1*nx1 + (ny1 - iy - 1)*nx1 + ix] = 0;
                            tmp1[iz*ny1*nx1 + iy*nx1 + ix] = 0;
                            tmp1[iz*ny1*nx1 + (ny1 - iy - 1)*nx1 + ix] = 0;
                        }
                    }
                }

                for (int iz = 0; iz < kw; ++iz) {
                    for (int iy = 0; iy < ny1; ++iy) {
                        for (int ix = 0; ix < nx1; ++ix) {
                            tmp0[iz*ny1*nx1 + iy*nx1 + ix] = 0;
                            tmp0[(nz1 - iz - 1)*ny1*nx1 + iy*nx1 + ix] = 0;
                            tmp1[iz*ny1*nx1 + iy*nx1 + ix] = 0;
                            tmp1[(nz1 - iz - 1)*ny1*nx1 + iy*nx1 + ix] = 0;
                        }
                    }
                }
            }

            for (int iz = 0; iz < nz; ++iz) {
                for (int iy = 0; iy < ny; ++iy) {
                    for (int ix = 0; ix < nx; ++ix) {
                        tmp0[(iz + kw)*ny1*nx1 + (iy + kw)*nx1 + (ix + kw)] = src[iz*ny*nx + iy*nx + ix];
                    }
                }
            }

            // x
            {
                std::atomic<int> curz(kw);
                std::vector<std::thread> workers(nthreads);
                for (auto & worker : workers) {
                    worker = std::thread([&curz, nx, ny, nz, kw, tmp0, tmp1, nx1, ny1, nz1, inorm, k]() {
                        while (true) {
                            int iz = curz.fetch_add(1);
                            if (iz >= nz + kw) break;
                            for (int iy = kw; iy < ny + kw; ++iy) {
                                for (int ix = kw; ix < nx + kw; ++ix) {
                                    float v = 0.0;
                                    for (int ik = -kw; ik <= kw; ++ik) {
                                        v += tmp0[iz*ny1*nx1 + iy*nx1 + (ix + ik)]*k[ik + kw];
                                    }
                                    tmp1[iz*ny1*nx1 + iy*nx1 + ix] = v*inorm;
                                }
                            }
                        }
                    });
                }
                for (auto & worker : workers) worker.join();
            }

            // y
            {
                std::atomic<int> curz(kw);
                std::vector<std::thread> workers(nthreads);
                for (auto & worker : workers) {
                    worker = std::thread([&curz, nx, ny, nz, kw, tmp0, tmp1, nx1, ny1, nz1, inorm, k]() {
                        while (true) {
                            int iz = curz.fetch_add(1);
                            if (iz >= nz + kw) break;
                            for (int ix = kw; ix < nx + kw; ++ix) {
                                for (int iy = kw; iy < ny + kw; ++iy) {
                                    float v = 0.0;
                                    for (int ik = -kw; ik <= kw; ++ik) {
                                        v += tmp1[iz*ny1*nx1 + (iy + ik)*nx1 + ix]*k[ik + kw];
                                    }
                                    tmp0[iz*ny1*nx1 + iy*nx1 + ix] = v*inorm;
                                }
                            }
                        }
                    });
                }
                for (auto & worker : workers) worker.join();
            }

            // z
            {
                std::atomic<int> cury(kw);
                std::vector<std::thread> workers(nthreads);
                for (auto & worker : workers) {
                    worker = std::thread([&cury, nx, ny, nz, kw, tmp0, dst, nx1, ny1, nz1, inorm, k]() {
                        while (true) {
                            int iy = cury.fetch_add(1);
                            if (iy >= ny + kw) break;
                            for (int ix = kw; ix < nx + kw; ++ix) {
                                for (int iz = kw; iz < nz + kw; ++iz) {
                                    float v = 0.0;
                                    for (int ik = -kw; ik <= kw; ++ik) {
                                        v += tmp0[(iz + ik)*ny1*nx1 + iy*nx1 + ix]*k[ik + kw];
                                    }
                                    dst[(iz - kw)*ny*nx + (iy - kw)*nx + (ix - kw)] = v*inorm;
                                }
                            }
                        }
                    });
                }
                for (auto & worker : workers) worker.join();
            }

            if (wdelete) {
                delete [] work;
            }

            return true;
        }

    template <typename T>
        bool gaussian_filter_3d(int nx, int ny, int nz, const T * src, T * dst, float sigma, int nthreads, int nw, T * work, bool wzero) {
            if (sigma <= 0.0) {
                for (int i = 0; i < nx*ny*nz; ++i) dst[i] = src[i];
                return true;
            }

            int nk = std::max(1.0, std::sqrt(-2.0*std::log(0.10))*sigma);

            std::vector<float> k(2*nk + 1);
            for (int i = 1; i <= nk; ++i) {
                float x = std::exp(-(i*i)/(2.0*sigma*sigma));
                k[nk - i] = k[nk + i] = x;
            }
            k[nk] = 1;

            return convolve_3d(nx, ny, nz, src, dst, k.size(), k.data(), nthreads, nw, work, wzero);
        }

    template <typename T>
        bool scale_li_3d(int snx, int sny, int snz, const T * src, float sx, float sy, float sz, int & dnx, int & dny, int & dnz, std::vector<T> & dst, int nthreads) {
            if (snx <= 0) return false;
            if (sny <= 0) return false;
            if (snz <= 0) return false;
            if (src == nullptr) return false;
            if (sx <= 0.0f) return false;
            if (sy <= 0.0f) return false;
            if (sz <= 0.0f) return false;

            if (src == dst.data()) return false;

            dnx = std::round(snx*sx);
            dny = std::round(sny*sy);
            dnz = std::round(snz*sz);

            if (dnx <= 0) return false;
            if (dny <= 0) return false;
            if (dnz <= 0) return false;

            dst.resize(dnx*dny*dnz);

            {
                std::atomic<int> curz(0);
                std::vector<std::thread> workers(nthreads);
                for (auto & worker : workers) {
                    worker = std::thread([&curz, snx, sny, snz, src, dnx, dny, dnz, &dst]() {
                        float cx = ((float)(snx))/dnx;
                        float cy = ((float)(sny))/dny;
                        float cz = ((float)(snz))/dnz;

                        while (true) {
                            int iz = curz.fetch_add(1);
                            if (iz >= dnz) break;
                            float fz = ((float)(iz) + 0.5f)*cz;

                            int iz1 = (int)(fz + 0.5f);
                            int iz0 = iz1 - 1;
                            if (iz0 < 0) ++iz0;
                            if (iz1 >= snz) --iz1;
                            fz -= (0.5f + (float)(iz0));
                            for (int iy = 0; iy < dny; ++iy) {
                                float fy = ((float)(iy) + 0.5f)*cy;

                                int iy1 = (int)(fy + 0.5f);
                                int iy0 = iy1 - 1;
                                if (iy0 < 0) ++iy0;
                                if (iy1 >= sny) --iy1;
                                fy -= (0.5f + (float)(iy0));
                                for (int ix = 0; ix < dnx; ++ix) {
                                    float fx = ((float)(ix) + 0.5f)*cx;

                                    int ix1 = (int)(fx + 0.5f);
                                    int ix0 = ix1 - 1;
                                    if (ix0 < 0) ++ix0;
                                    if (ix1 >= snx) --ix1;
                                    fx -= (0.5f + (float)(ix0));

                                    auto v000 = src[iz0*sny*snx + iy0*snx + ix0];
                                    auto v001 = src[iz0*sny*snx + iy0*snx + ix1];
                                    auto v010 = src[iz0*sny*snx + iy1*snx + ix0];
                                    auto v011 = src[iz0*sny*snx + iy1*snx + ix1];
                                    auto v100 = src[iz1*sny*snx + iy0*snx + ix0];
                                    auto v101 = src[iz1*sny*snx + iy0*snx + ix1];
                                    auto v110 = src[iz1*sny*snx + iy1*snx + ix0];
                                    auto v111 = src[iz1*sny*snx + iy1*snx + ix1];

                                    auto v00 = v000 + fz*(v100 - v000);
                                    auto v01 = v001 + fz*(v101 - v001);
                                    auto v10 = v010 + fz*(v110 - v010);
                                    auto v11 = v011 + fz*(v111 - v011);

                                    auto v0 = v00 + fy*(v10 - v00);
                                    auto v1 = v01 + fy*(v11 - v01);

                                    dst[iz*dny*dnx + iy*dnx + ix] = v0 + fx*(v1 - v0);
                                }
                            }
                        }
                    });
                }
                for (auto & worker : workers) worker.join();
            }

            return true;
        }

    template <typename T>
        bool scale_li_isotropic_3d(int snx, int sny, int snz, const T * src, float s, int & dnx, int & dny, int & dnz, std::vector<T> & dst, int nthreads) {
            return scale_li_3d(snx, sny, snz, src, s, s, s, dnx, dny, dnz, dst, nthreads);
        }

    template <typename T>
        bool scale_nn_3d(int snx, int sny, int snz, const T * src, float sx, float sy, float sz, int & dnx, int & dny, int & dnz, std::vector<T> & dst, int nthreads) {
            if (snx <= 0) return false;
            if (sny <= 0) return false;
            if (snz <= 0) return false;
            if (src == nullptr) return false;
            if (sx <= 0.0f) return false;
            if (sy <= 0.0f) return false;
            if (sz <= 0.0f) return false;

            if (src == dst.data()) return false;

            dnx = std::round(snx*sx);
            dny = std::round(sny*sy);
            dnz = std::round(snz*sz);

            if (dnx <= 0) return false;
            if (dny <= 0) return false;
            if (dnz <= 0) return false;

            dst.resize(dnx*dny*dnz);

            {
                std::atomic<int> curz(0);
                std::vector<std::thread> workers(nthreads);
                for (auto & worker : workers) {
                    worker = std::thread([&curz, snx, sny, snz, src, dnx, dny, dnz, &dst]() {
                        float cx = ((float)(snx))/dnx;
                        float cy = ((float)(sny))/dny;
                        float cz = ((float)(snz))/dnz;

                        while (true) {
                            int iz = curz.fetch_add(1);
                            if (iz >= dnz) break;
                            float fz = ((float)(iz) + 0.5f)*cz;

                            int iz0 = std::round(fz - 0.5f);
                            for (int iy = 0; iy < dny; ++iy) {
                                float fy = ((float)(iy) + 0.5f)*cy;

                                int iy0 = std::round(fy - 0.5f);
                                fy -= (0.5f + (float)(iy0));
                                for (int ix = 0; ix < dnx; ++ix) {
                                    float fx = ((float)(ix) + 0.5f)*cx;

                                    int ix0 = std::round(fx - 0.5f);

                                    dst[iz*dny*dnx + iy*dnx + ix] = src[iz0*sny*snx + iy0*snx + ix0];
                                }
                            }
                        }
                    });
                }
                for (auto & worker : workers) worker.join();
            }

            return true;
        }

    template <typename T>
        bool scale_nn_isotropic_3d(int snx, int sny, int snz, const T * src, float s, int & dnx, int & dny, int & dnz, std::vector<T> & dst, int nthreads) {
            return scale_nn_3d(snx, sny, snz, src, s, s, s, dnx, dny, dnz, dst, nthreads);
        }

}

#endif
