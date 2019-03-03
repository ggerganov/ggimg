/*! \file ggimg-mt.cpp
 *  \brief Some tests of ggimg (multi-threaded operations)
 *  \author Georgi Gerganov
 */

#define GGIMG_MT
#include "ggimg.h"

#include <cstdio>

int main(int argc, char ** argv) {
    printf("Usage: %s img.raw\n", argv[0]);
    if (argc < 2) return -1;

    printf("All done\n");

    return 0;
}
