/*! \file ggimg.cpp
 *  \brief Some tests of ggimg
 *  \author Georgi Gerganov
 */

#include "ggimg.h"

#include <cstdio>
#include <fstream>
#include <vector>
#include <string>

bool read_ppm(const char * fname, int & nx, int & ny, std::vector<uint8_t> & img) {
    std::ifstream fin(fname, std::ios::binary);
    if (fin.good() == false) return false;

    int ifmt = -1;
    int imax = -1;

    {
        std::string fmt;
        fin >> fmt;
        //if (fmt == "P1") ifmt = 1;
        //if (fmt == "P2") ifmt = 2;
        //if (fmt == "P3") ifmt = 3;
        //if (fmt == "P4") ifmt = 4;
        //if (fmt == "P5") ifmt = 5;
        if (fmt == "P6") ifmt = 6;
    }

    if (ifmt == -1) return false;

    fin >> nx >> ny >> imax;

    if (fin.eof()) return false;

    img.resize(3*nx*ny);
    fin.read((char *)(img.data()), 1);
    fin.read((char *)(img.data()), 3*nx*ny);

    if (fin.eof()) return false;

    return true;
}

bool write_ppm(const char * fname, int nx, int ny, const std::vector<uint8_t> & img, int bpp = 3) {
    std::ofstream fout(fname, std::ios::binary);

    if (bpp == 1) fout << "P5\n";
    if (bpp == 3) fout << "P6\n";

    fout << nx << " " << ny << "\n255\n";

    if (bpp == 1) fout.write((char *)(img.data()), nx*ny);
    if (bpp == 3) fout.write((char *)(img.data()), 3*nx*ny);

    fout << "\n";

    fout.close();

    return true;
}

int main(int argc, char ** argv) {
    printf("Usage: %s img.raw\n", argv[0]);
    if (argc < 2) return -1;

    int nx = 0;
    int ny = 0;

    std::string fname(argv[1]);

    using Image2D = std::vector<uint8_t>;

    Image2D img_tmp;
    Image2D img_rgb;
    Image2D img_luma601;
    Image2D img_luma709;
    Image2D img_grayscale;
    Image2D img_normalize;
    Image2D img_normalize_hist;
    Image2D img_gaussian_filter;
    Image2D img_median_filter;
    Image2D img_gradient_x;
    Image2D img_gradient_y;
    Image2D img_gradient_xy;
    Image2D img_scale_nn;
    Image2D img_scale_li;
    Image2D img_transform_homography_gray_nn;
    Image2D img_transform_homography_rgb_nn;

    {
        printf("[+] Read/write ppm\n");

        if (read_ppm(fname.c_str(), nx, ny, img_rgb) == false) {
            printf("Failed to read image '%s'\n", fname.c_str());
            return -1;
        }

        if (write_ppm("ggimg_original.ppm", nx, ny, img_rgb, 3) == false) {
            printf("Failed to write P6 ppm\n");
            return -1;
        }
    }

    {
        printf("[+] ggimg::rgb_to_luma601\n");

        img_luma601.resize(nx*ny);
        if (ggimg::rgb_to_luma601_2d(nx, ny, img_rgb.data(), img_luma601.data()) == false) {
            printf("Failed ggimg::rgb_to_luma601\n");
            return -1;
        }

        write_ppm("ggimg_luma601.ppm", nx, ny, img_luma601, 1);
    }

    {
        printf("[+] ggimg::rgb_to_luma709\n");

        img_luma709.resize(nx*ny);
        if (ggimg::rgb_to_luma709_2d(nx, ny, img_rgb.data(), img_luma709.data()) == false) {
            printf("Failed ggimg::rgb_to_luma709\n");
            return -1;
        }

        write_ppm("ggimg_luma709.ppm", nx, ny, img_luma709, 1);
    }

    {
        printf("[+] ggimg::rgb_to_gray_2d\n");

        img_grayscale.resize(nx*ny);
        if (ggimg::rgb_to_gray_2d(nx, ny, img_rgb.data(), img_grayscale.data()) == false) {
            printf("Failed ggimg::rgb_to_gray_2d\n");
            return -1;
        }

        write_ppm("ggimg_grayscale.ppm", nx, ny, img_grayscale, 1);
    }

    {
        printf("[+] ggimg::gray_to_rgb_2d\n");

        img_tmp.resize(3*nx*ny);
        if (ggimg::gray_to_rgb_2d(nx, ny, img_grayscale.data(), img_tmp.data()) == false) {
            printf("Failed ggimg::gray_to_rgb_2d\n");
            return -1;
        }

        write_ppm("ggimg_rgb.ppm", nx, ny, img_tmp, 3);
    }

    {
        printf("[+] ggimg::normalize_2d\n");

        img_normalize.resize(nx*ny);
        if (ggimg::normalize_2d(nx, ny, img_grayscale.data(), (uint8_t) 100, (uint8_t) 200, img_normalize.data()) == false) {
            printf("Failed ggimg::normalize_2d\n");
            return -1;
        }

        write_ppm("ggimg_normalize.ppm", nx, ny, img_normalize, 1);
    }

    {
        printf("[+] ggimg::normalize_hist_2d\n");

        img_normalize_hist.resize(nx*ny);
        if (ggimg::normalize_hist_2d(nx, ny, img_grayscale.data(), img_normalize_hist.data(), (uint8_t) 255) == false) {
            printf("Failed ggimg::normalize_hist_2d\n");
            return -1;
        }

        write_ppm("ggimg_normalize_hist.ppm", nx, ny, img_normalize_hist, 1);
    }

    {
        printf("[+] ggimg::gradient_sobel_2d(mode = 1)\n");

        img_gradient_x.resize(nx*ny);
        if (ggimg::gradient_sobel_2d(1, nx, ny, img_grayscale.data(), (uint8_t) 255, img_gradient_x.data()) == false) {
            printf("Failed ggimg::gradient_sobel_2d\n");
            return -1;
        }

        write_ppm("ggimg_gradient_x.ppm", nx, ny, img_gradient_x, 1);
    }

    {
        printf("[+] ggimg::gradient_sobel_2d(mode = 2)\n");

        img_gradient_y.resize(nx*ny);
        if (ggimg::gradient_sobel_2d(2, nx, ny, img_grayscale.data(), (uint8_t) 255, img_gradient_y.data()) == false) {
            printf("Failed ggimg::gradient_sobel_2d\n");
            return -1;
        }

        write_ppm("ggimg_gradient_y.ppm", nx, ny, img_gradient_y, 1);
    }

    {
        printf("[+] ggimg::gradient_sobel_2d(mode = 0)\n");

        img_gradient_xy.resize(nx*ny);
        if (ggimg::gradient_sobel_2d(0, nx, ny, img_grayscale.data(), (uint8_t) 255, img_gradient_xy.data()) == false) {
            printf("Failed ggimg::gradient_sobel_2d\n");
            return -1;
        }

        write_ppm("ggimg_gradient.ppm", nx, ny, img_gradient_xy, 1);
    }

    {
        printf("[+] ggimg::gaussian_filter_2d\n");

        img_gaussian_filter.resize(nx*ny);
        if (ggimg::gaussian_filter_2d(nx, ny, img_grayscale.data(), img_gaussian_filter.data(), 3.0f) == false) {
            printf("Failed ggimg::gaussian_filter_2d\n");
            return -1;
        }

        write_ppm("ggimg_gaussian_filter.ppm", nx, ny, img_gaussian_filter, 1);
    }

    {
        printf("[+] ggimg::median_filter_2d\n");

        img_median_filter.resize(nx*ny);
        if (ggimg::median_filter_2d(nx, ny, img_grayscale.data(), img_median_filter.data(), 5) == false) {
            printf("Failed ggimg::median_filter_2d\n");
            return -1;
        }

        write_ppm("ggimg_median_filter.ppm", nx, ny, img_median_filter, 1);
    }

    {
        printf("[+] ggimg::scale_nn_2d\n");

        int snx = -1;
        int sny = -1;
        if (ggimg::scale_nn_2d(nx, ny, img_grayscale.data(), 0.33f, 0.66f, snx, sny, img_scale_nn) == false) {
            printf("Failed ggimg::scale_nn_2d\n");
            return -1;
        }

        write_ppm("ggimg_scale_nn.ppm", snx, sny, img_scale_nn, 1);
    }

    {
        printf("[+] ggimg::scale_li_2d\n");

        int snx = -1;
        int sny = -1;
        if (ggimg::scale_li_2d(nx, ny, img_grayscale.data(), 0.33f, 0.66f, snx, sny, img_scale_li) == false) {
            printf("Failed ggimg::scale_li_2d\n");
            return -1;
        }

        write_ppm("ggimg_scale_li.ppm", snx, sny, img_scale_li, 1);
    }

    {
        printf("[+] ggimg::transform_homography_gray_nn\n");

        img_transform_homography_gray_nn.resize(nx*ny);
        if (ggimg::transform_homography_gray_nn(nx, ny, img_grayscale.data(), { 0.68, -0.48, 0.38, 0.87, 0.81, -0.43, -0.15, 0.12, 1.00 }, nx, ny, img_transform_homography_gray_nn.data()) == false) {
            printf("Failed ggimg::transform_homography_gray_nn \n");
            return -1;
        }

        write_ppm("ggimg_transform_homography_gray_nn.ppm", nx, ny, img_transform_homography_gray_nn, 1);
    }

    {
        printf("[+] ggimg::transform_homography_rgb_nn\n");

        img_transform_homography_rgb_nn.resize(3*nx*ny);
        if (ggimg::transform_homography_rgb_nn(nx, ny, img_rgb.data(), { 0.68, -0.48, 0.38, 0.87, 0.81, -0.43, -0.15, 0.12, 1.00 }, nx, ny, img_transform_homography_rgb_nn.data()) == false) {
            printf("Failed ggimg::transform_homography_rgb_nn \n");
            return -1;
        }

        write_ppm("ggimg_transform_homography_rgb_nn.ppm", nx, ny, img_transform_homography_rgb_nn, 3);
    }

    printf("All done\n");

    return 0;
}
