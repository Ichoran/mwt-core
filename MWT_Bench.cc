/**********************************************************
* This file implements a worm-tracking benchmark for the  *
* MWT library.                                            *
***********************************************************
* Copyright (c) 2017 by Rex Kerr and Calico Life Sciences *
***********************************************************
* This file is part of the core Multi-Worm Tracker code.  *
* See LICENCE.txt for full copyright information.         *
**********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include <vector>

#include "MWT_Geometry.h"
#include "MWT_Image.h"
#include "MWT_Library.h"

void deltize(timeval &t) {
    timeval now; gettimeofday(&now, NULL);
    t.tv_sec = now.tv_sec - t.tv_sec;
    t.tv_usec = now.tv_usec - t.tv_usec;
    if (t.tv_usec < 0) { t.tv_usec += 1000000; t.tv_sec -= 1; }
    if (t.tv_sec < 0) { t.tv_sec = 0; t.tv_usec = 0; }
}

std::vector<Image*> names_to_data(int argn, char *argv[], int dep) {
    std::vector<Image*> images;
    images.reserve(argn-1);
    unsigned char *buffer = new unsigned char[512*512];
    for (int i = 1; i < argn; i++) {
        FILE *f = fopen(argv[i], "rb");
        Image *im = new Image(Point(512, 512), false);
        im->depth = (dep < 1) ? 1 : ((dep > 14) ? 14: dep);
        auto n = fread((void*)buffer, 1, 512*512, f);
        unsigned char *b = buffer;
        short *c = im->pixels;
        if (im->depth > 8) {
            short m = 1<<(im->depth - 8);
            for (int j = 0; j < 512*512; j++) { *c = (short)(*b) * m; b++; c++; }
        }
        else {
            for (int j = 0; j < 512*512; j++) { *c = (short)(*b); b++; c++; }
        }
        if (n != 512*512) { printf("Wrong number of pixels: %ld in %s\n", n, argv[i]); exit(1); }
        fclose(f);
        images.push_back(im);
    }
    return images;
}

std::vector<Image8*> names_to_data8(int argn, char *argv[]) {
    std::vector<Image8*> images;
    images.reserve(argn-1);
    for (int i = 1; i < argn; i++) {
        FILE *f = fopen(argv[i], "rb");
        Image8 *im = new Image8(Point(512, 512), false);
        auto n = fread((void*)im->pixels, 1, 512*512, f);
        if (n != 512*512) { printf("Wrong number of pixels: %ld in %s\n", n, argv[i]); exit(1); }
        fclose(f);
        images.push_back(im);
    }
    return images;
}

void set_the_date(TrackerLibrary &mwt, int h) {
    time_t t = time(NULL);
    tm tt; localtime_r(&t, &tt);
    mwt.setDate(h, tt.tm_year+1900, tt.tm_mon+1, tt.tm_mday, tt.tm_hour, tt.tm_min, tt.tm_sec);
}

void yell(int i, int j, int e) {
    if (i != j) printf("%d!%d ", e, j);
}

#define EIGHT

int main(int argn, char *argv[]) {
    timeval t0; gettimeofday(&t0, NULL);
#ifdef EIGHT
    auto im = new Image8(Point(512, 512), false);
    auto images = names_to_data8(argn, argv);
#else
    int nbits = 8;
    auto im = new Image(Point(512, 512), false);
    auto images = names_to_data(argn, argv, nbits);
#endif
    deltize(t0);
    printf("Read %d files in %ld.%03lds\n", argn-1, t0.tv_sec, t0.tv_usec/1000);
    TrackerLibrary mwt;
    auto h = mwt.getNewHandle();
    int hh = h;
    hh = mwt.setCombineBlobs(h, true); yell(h, hh, -1);
#ifdef EIGHT
    hh = mwt.setImageInfo(h, 8, 512, 512); yell(h, hh, -2);
#else
    hh = mwt.setImageInfo(h, nbits, 512, 512); yell(h, hh, -2);
#endif
    set_the_date(mwt, h);
    hh = mwt.setOutput(h, "fake-data", "bench", true, false, false); yell(h, hh, -3);
    hh = mwt.setRectangle(h, 1, 511, 1, 511); yell(h, hh, -4);
    hh = mwt.setRefIntensityThreshold(h, 1, 2);
    hh = mwt.setDancerBorderSize(h, 8); yell(h, hh, -5);
#ifdef EIGHT
    hh = mwt.setObjectIntensityThresholds(h, 240, 230); yell(h, hh, -6);
#else
    hh = mwt.setObjectIntensityThresholds(h, 960/4, 920/4); yell(h, hh, -6);
#endif
    hh = mwt.setObjectSizeThresholds(h, 5, 10, 2000, 3000); yell(h, hh, -7);
    hh = mwt.setObjectPersistenceThreshold(h, 8); yell(h, hh, -8);
    hh = mwt.setAdaptationRate(h, 4); yell(h, hh, -9);
    hh = mwt.enableOutlining(h, true); yell(h, hh, -10);
    hh = mwt.enableSkeletonization(h, true); yell(h, hh, -11);

    timeval t1; gettimeofday(&t1, NULL);
    hh = mwt.beginOutput(h); yell(h, hh, -12);
    hh = mwt.setUpdateBandNumber(h, 12); yell(h, hh, -13);
    hh = mwt.setVelocityIntegrationTime(h, 0.5);  yell(h, hh, -14);
    for (unsigned int i = 0; i < images.size(); i++) {
        if (i < 1) { 
#ifdef EIGHT
            int n = mwt.scanObjects8(h, *images[i]); if (n < 0) yell(h, n, i);
#else
            int n = mwt.scanObjects(h, *images[i]); if (n < 0) yell(h, n, i);
#endif
            mwt.all_trackers[h]->performance.background->writeTiff("bg0.tiff");
            mwt.all_trackers[h]->performance.foreground->writeTiff("fg0.tiff");
        }
        else {
#ifdef EIGHT
            hh = mwt.loadImage8(h, *images[i], 0.03f*i); yell(h, hh, 2*i);
#else
            hh = mwt.loadImage(h, *images[i], 0.03f*i); yell(h, hh, 2*i);
#endif
            if (i == 1) {
                mwt.all_trackers[h]->performance.background->writeTiff("bg1a.tiff");
                mwt.all_trackers[h]->performance.foreground->writeTiff("fg1a.tiff");
            }
            int n = mwt.processImage(h); yell(h, hh, 2*i+1); if (n < 0) yell(h, n, 2*i+1);
            if (i == 1) {
                mwt.all_trackers[h]->performance.background->writeTiff("bg1b.tiff");
                mwt.all_trackers[h]->performance.foreground->writeTiff("fg1b.tiff");
            }
        }
    }
    mwt.complete(h);
    deltize(t1);
    printf("Processed %ld images in %ld.%03lds\n", images.size(), t1.tv_sec, t1.tv_usec/1000);
}
