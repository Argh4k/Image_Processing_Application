//
// Created by Maciej on 01.03.2019.
//

#ifndef IMAGEPROCESSING_IMAGEUTILS_H
#define IMAGEPROCESSING_IMAGEUTILS_H

#include "CImg.h"
#include <string>
#include <fstream>

namespace cimg = cimg_library;

template<typename T>
cimg::CImg<T> loadImageFromFile(const std::string& fileName) {
    cimg::CImg<T> image;
    std::ifstream file(fileName);
    if(file.good()) {
        image.load(fileName);
        return image;
    } else {
        throw std::invalid_argument("Could not open file");
    }
}

template <typename T>
void SaveImage(const cimg::CImg<T> & image, const std::string& fileName) {
    image.save(fileName.c_str());
}

#endif //IMAGEPROCESSING_IMAGEUTILS_H
