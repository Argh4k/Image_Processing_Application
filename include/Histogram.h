//
// Created by Maciej on 01.03.2019.
//

#ifndef IMAGEPROCESSING_HISTOGRAM_H
#define IMAGEPROCESSING_HISTOGRAM_H

#include "CImg.h"
#include <array>

namespace cimg = cimg_library;

static constexpr int HISTOGRAM_SIZE = 256; //max possible level of intensity

template <typename T>
std::array<int, HISTOGRAM_SIZE> createHistogram(cimg::CImg<T> & image, int channel)
{
    std::array<int, HISTOGRAM_SIZE> histogram {0};
    for (int x = 0; x < image.width(); x++)
    {
        for (int y = 0; y < image.height(); y++)
        {
            histogram[(int)image(x, y,0, channel)]++;
        }
    }
    return histogram;
}

template <typename T>
cimg::CImg<uint8_t> visualizeHistogram(const cimg::CImg<T> & image, int imageChannel)
{
    auto histogram = createHistogram(image, imageChannel);
    int maxHistogramHeight = 0;
    for (int i = 0; i < HISTOGRAM_SIZE; i++)
    {
        if (histogram[i] > maxHistogramHeight)
            maxHistogramHeight = histogram[i];
    }
    cimg::CImg<uint8_t> histogramPicture(HISTOGRAM_SIZE, maxHistogramHeight/10 + 1, 1, 3);
    for (int x = 0; x < HISTOGRAM_SIZE; x+=2)
    {
        for (float y = (maxHistogramHeight)/10 ; y > maxHistogramHeight/10-histogram[x]/10 ; y--)
        {
            histogramPicture(x, y, 0, imageChannel) = 255;
            histogramPicture(x+1, y, 0, imageChannel) = 255;
        }
    }
    return histogramPicture;
}

#endif //IMAGEPROCESSING_HISTOGRAM_H
