
#ifndef IMAGE_H_
#define IMAGE_H_

#include "IMAGE.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include "CImg.h"


using namespace cimg_library;
using namespace std;


constexpr int PIXEL_RANGE = 256;
constexpr int MAX_PIXEL_VALUE = 255;


static void showErrorData (double clean_picture_data, double median_filter_data, double gmean_filter_data) {

	cout << "\nBefore processing " << clean_picture_data << endl;
	cout << "After median filter " << median_filter_data << endl;
	cout << "After geometric filter " << gmean_filter_data << endl;
}


template<typename T>
T getMaximumPixelValue(const CImg<T> & image) {

	T maximumValue = 0;
	for (int c = 0; c < image.spectrum(); c++) {
		for (int x = 0; x < image.width(); x++) {
			for (int y = 0; y < image.height(); y++) {
				if (maximumValue < image(x, y, 0, c)) {
					if (image(x, y, 0, c) == MAX_PIXEL_VALUE)
						return MAX_PIXEL_VALUE;
					else
						maximumValue = image(x, y, 0, c);
				}
			}
		}
	}
	return maximumValue;
}

static uint8_t clampToPixelValue(int newValue) {
	if(newValue < 0)
		return 0;
	else if (newValue > 255)
		return 255;
	return static_cast<uint8_t>(newValue);
}

// FUNCTION CREATING LUT DEPENDIND ON THE OPERATION
std::array<uint8_t, PIXEL_RANGE> createLookUpTable(uint8_t(*operation)(uint8_t))
{
	std::array<uint8_t, PIXEL_RANGE> lookUpTable{};
	for(size_t i = 0; i < lookUpTable.size(); ++i) {
		lookUpTable[i] = operation(i);
	}
	return lookUpTable;;
}

// POSSIBLE OPERATIONS FOR LUT

uint8_t changeBrightness(uint8_t oldValue, uint8_t increment) {
	return clampToPixelValue(oldValue + increment);
}

uint8_t changeContrast(uint8_t oldValue, uint8_t ratio) {
	return clampToPixelValue((oldValue - 128)*ratio + 128);
}

uint8_t makeNegative(uint8_t oldValue) {
	return clampToPixelValue(255 - oldValue);
}

void applyOperationToImage(CImg<uint8_t> & image, uint8_t(*operation)(uint8_t)) {
	auto lookUpTable = createLookUpTable(operation);
	for (int c = 0; c < image.spectrum(); c++) {
		for (int x = 0; x < image.width(); x++) {
			for (int y = 0; y < image.height(); y++) {
				image(x, y, 0, c) = lookUpTable[image(x, y, 0, c)];
			}
		}
	}
}

//FLIPS

template <typename T>
void verticalFlip(CImg<T> &image) {
	for (int c = 0; c < image.spectrum(); c++) {
		for (int x = 0; x < image.width(); x++) {
			for (int y = 0; y < image.height() / 2; y++) {
				swap(image(x, y, 0, c), image(x, image.height() - y - 1, 0, c));
			}
		}
	}
}

template <typename T>
void horizontalFlip(CImg<T> & image) {

	for (int y = 0; y < image.height(); y++)
	{
		for (int x = 0; x < image.width() / 2; x++)
		{
			for (int c = 0; c < image.spectrum(); c++) {
				swap(image(x, y, 0, c), image(image.width() - x - 1, y, 0, c));
			}
		}
	}

}

template<typename T>
void diagonalFlip(CImg<T>& image) {

	for (int x = 0; x < image.width(); x++)
	{
		for (int y = 0; y < image.height() / 2; y++)
		{
			for (int c = 0; c < image.spectrum(); c++)
			{
				swap(image(x, y, 0, c), image(image.width() - 1 - x, image.height() - 1 - y, 0, c));
			}
		}
	}
}

CImg<float> shrink(const CImg<float> & image) {

	CImg<float> shrinked_image((image.width() / 2), (image.height() / 2), 1, 3);

	for (int x = 0; x < shrinked_image.width(); x++)  //to not run out of table
	{
		for (int y = 0; y < shrinked_image.height(); y++)
		{
			for (int c = 0; c < image.spectrum(); c++)
			{

				shrinked_image(x, y, 0, c) = image(2 * x, 2 * y, 0, c); //every second pixel (0,0)(2,2)(4,4) etc

			}
		}
	}

	return shrinked_image;
}

CImg<float> enlarge(const CImg<float> & image) {

	CImg<float> enlarged_image((image.width() * 2), (image.height() * 2), 1, 3);

	for (int x = 0; x < enlarged_image.width() - 1; x += 2)
	{
		for (int y = 0; y < enlarged_image.height() - 1; y += 2)
		{
			for (int c = 0; c < image.spectrum(); c++)
			{

				enlarged_image(x, y, 0, c) = image(x / 2, y / 2, 0, c);
				enlarged_image(x + 1, y + 1, 0, c) = image(x / 2, y / 2, 0, c);

			}
		}
	}
	return enlarged_image;
}

template<typename T>
T median(const CImg<T> & image, int x, int y, int c, int maskSize)
{
	std::vector<T> values;
	values.reserve(maskSize * maskSize);
	int maskRange = maskSize / 2;
	for (int i = -maskRange; i <= maskRange; i++)
	{
		for (int j = -maskRange; j <= maskRange; j++)
		{
			values.push_back(image(x + i, y + j, 0, c));
		}
	}
	std::sort(values.begin(), values.end());
	return values.at(values.size()/2);
}

template<typename T>
T geometricMean(const CImg<T>& image, int x, int y, int c, int maskSize)
{
	double sum = 1;
	int maskRange = maskSize / 2;
	for (int i = -maskRange; i <= maskRange; i++)
	{
		for (int j = -maskRange; j <= maskRange; j++)
		{
			sum*=image(x + i, y + j, 0, c);
		}
	}
	double root = 1.0f / (maskSize*maskSize);
	return pow(sum, root);
}

template<typename T>
CImg<T> medianfilter(const CImg<T> & image, int maskSize) {
	CImg<T> filteredImage(image);
	int maskRange = maskSize / 2;
	if(!(maskSize % 2)) {
		throw std::invalid_argument("Mask must be odd");
	}
	for (int c = 0; c < image.spectrum(); c++) {
		for (int x = maskRange; x < image.width() - maskRange; x++) {
			for (int y = maskRange; y < image.height() - maskRange; y++) {
				filteredImage(x, y, 0, c) = median(image, x, y, c);
			}
		}
	}
	return filteredImage;
}

template<typename T>
CImg<T> * geometricfilter(const CImg<T> & image, int maskSize) {
	CImg<T> filteredImage(image);
	int maskRange = maskSize / 2;
	if(!(maskSize % 2)) {
		throw std::invalid_argument("Mask must be odd");
	}
	for (int c = 0; c < image.spectrum(); c++) {
		for (int x = maskRange; x < image.width() - maskRange; x++) {
			for (int y = maskRange; y < image.height() - maskRange; y++) {
				filteredImage(x, y, 0, c) = geometricMean(image, x, y, c);
			}
		}
	}
	return filteredImage;
}

template <typename T>
void SaveImage(const CImg<T> & image) {
	string name("output.bmp");
	image.save(name.c_str());
}

template<typename T>
void computeMeanSquareError(const CImg<T> & imageWithoutNoise, const CImg<T> & imageWithNoise)
{
	CImg<T> imageAfterMedianFilter = medianfilter(imageWithNoise);
	CImg<T> imageAfterGeometricFilter = geometricfilter(imageWithNoise);
	double meanError = 0;
	double meanErrorAfterMedianFilter = 0;
	double meanErrorAfterGeometric = 0;

	for (int c = 0; c < imageWithoutNoise.spectrum(); c++)
	{
		for (int x = 0; x < imageWithoutNoise.width(); x++)
		{
			for (int y = 0; y < imageWithoutNoise.height(); y++)
			{
				meanError += pow((imageWithoutNoise(x, y, 0, c) - imageWithNoise(x, y, 0, c)), 2);
				meanErrorAfterGeometric += pow((imageWithoutNoise(x, y, 0, c) - imageAfterGeometricFilter(x, y, 0, c)), 2);
				meanErrorAfterMedianFilter += pow((imageWithoutNoise(x, y, 0, c) - imageAfterMedianFilter(x, y, 0, c)), 2);
			}
		}

	}


	double coefficient = (1.0 / (imageWithoutNoise.width()*imageWithoutNoise.height() * imageWithoutNoise.spectrum()));
	meanError = coefficient*meanError;
	meanErrorAfterGeometric = coefficient*meanErrorAfterGeometric;
	meanErrorAfterMedianFilter = coefficient*meanErrorAfterMedianFilter;

	cout << "\nMean square error" << endl;
	showErrorData(meanError, meanErrorAfterMedianFilter, meanErrorAfterGeometric);

}

template<typename T>
void computePeakMeanSquareError(const CImg<T> &imageWithoutNoise, const CImg<T> &imageWithNoise)
{
	CImg<T> imageAfterMedianFilter = medianfilter(imageWithNoise);
	CImg<T> imageAfterGeometricFilter = geometricfilter(imageWithNoise);
	double meanError = 0;
	double meanErrorAfterMedianFilter = 0;
	double meanErrorAfterGeometric = 0;

	for (int c = 0; c < imageWithoutNoise.spectrum(); c++)
	{
		for (int x = 0; x < imageWithoutNoise.width(); x++)
		{
			for (int y = 0; y < imageWithoutNoise.height(); y++)
			{
				meanError += pow((imageWithoutNoise(x, y, 0, c) - imageWithNoise(x, y, 0, c)), 2);
				meanErrorAfterGeometric += pow((imageWithoutNoise(x, y, 0, c) - imageAfterGeometricFilter(x, y, 0, c)), 2);
				meanErrorAfterMedianFilter += pow((imageWithoutNoise(x, y, 0, c) - imageAfterMedianFilter(x, y, 0, c)), 2);
			}
		}

	}
	int pow_max_value = pow(getMaximumPixelValue(imageWithoutNoise), 2);
	double coefficient = (1.0 / (imageWithoutNoise.width()*imageWithoutNoise.height() * imageWithoutNoise.spectrum()));
	meanError = (coefficient*meanError) / pow_max_value;
	meanErrorAfterMedianFilter = (coefficient*meanErrorAfterMedianFilter) / pow_max_value;
	meanErrorAfterGeometric = (coefficient*meanErrorAfterGeometric) / pow_max_value;

	cout << "\nPeak mean square error" << endl;
	showErrorData(meanError,meanErrorAfterMedianFilter,meanErrorAfterGeometric);

}

template <typename T>
void calculateSignalToNoiseError(const CImg<T> & imageWithoutNoise, const CImg<T> & imageWithNoise)
{
	CImg<T> imageAfterMedianFilter = medianfilter(imageWithNoise);
	CImg<T> imageAfterGeometricFilter = geometricfilter(imageWithNoise);
	double error = 0;
	double errorAfterMedianFilter = 0;
	double errorAfterGeometricFilter = 0;
	double sum = 0;

	for (size_t c = 0; c < imageWithoutNoise.spectrum(); c++)
	{
		for (size_t x = 0; x < imageWithoutNoise.width(); x++)
		{
			for (size_t y = 0; y < imageWithoutNoise.height(); y++)
			{
				sum += pow(imageWithoutNoise(x, y, 0, c), 2);
				error += pow((imageWithoutNoise(x, y, 0, c) - imageWithNoise(x, y, 0, c)), 2);
				errorAfterGeometricFilter += pow((imageWithoutNoise(x, y, 0, c) - imageAfterGeometricFilter(x, y, 0, c)), 2);
				errorAfterMedianFilter += pow((imageWithoutNoise(x, y, 0, c) - imageAfterMedianFilter(x, y, 0, c)), 2);
			}
		}

	}

	error = 10*log10(sum/error);
	errorAfterGeometricFilter = 10*log10(sum/ errorAfterGeometricFilter);
	errorAfterMedianFilter = 10*log10(sum / errorAfterMedianFilter);

	cout << "Signal to noise error" << endl;
	showErrorData(error,errorAfterMedianFilter,errorAfterGeometricFilter);

}

template<typename T>
void computePeakSignalToNoiseError(const CImg<T> & imageWithoutNoise, const CImg<T> & imageWithNoise)
{
	CImg<T> imageAfterMedianFilter = medianfilter(imageWithNoise);
	CImg<T> imageAfterGeometricFilter = geometricfilter(imageWithNoise);
	double error = 0;
	double errorAfterMedianFilter = 0;
	double errorAfterGeometricFilter = 0;
	double sum = 0;

	for (size_t c = 0; c < imageWithoutNoise.spectrum(); c++)
	{
		for (int x = 0; x < imageWithoutNoise.width(); x++)
		{
			for (int y = 0; y < imageWithoutNoise.height(); y++)
			{
				error += pow((imageWithoutNoise(x, y, 0, c) - imageWithNoise(x, y, 0, c)), 2);
				errorAfterGeometricFilter += pow((imageWithoutNoise(x, y, 0, c) - imageAfterGeometricFilter(x, y, 0, c)), 2);
				errorAfterMedianFilter += pow((imageWithoutNoise(x, y, 0, c) - imageAfterMedianFilter(x, y, 0, c)), 2);
			}
		}
	}
	double coefficient = (1.0 / (imageWithoutNoise.width()*imageWithoutNoise.height() * imageWithoutNoise.spectrum()));
	int pow_max_value = pow(getMaximumPixelValue(imageWithoutNoise), 2);
	error = 10 * log10(pow_max_value / (error*coefficient));
	errorAfterGeometricFilter = 10 * log10(pow_max_value / (errorAfterGeometricFilter*coefficient));
	errorAfterMedianFilter = 10 * log10(pow_max_value / (errorAfterMedianFilter*coefficient));

	cout << "Peak signal to noise error" << endl;
	showErrorData(error, errorAfterMedianFilter, errorAfterGeometricFilter);
}

template<typename T>
void calculateMaximumDifference(const CImg<T> & imageWithoutNoise, const CImg<T> & imageWithNoise)
{
	if(imageWithNoise.width() != imageWithoutNoise.width() ||
	   imageWithNoise.height() != imageWithoutNoise.height() ||
	   imageWithNoise.spectrum() != imageWithoutNoise.spectrum()) {
		throw std::invalid_argument("Images must have the same dimensions");
	}

	CImg<T> imageAfterMedianFilter = medianfilter(imageWithNoise);
	CImg<T> imageAfterGeometricFilter = geometricfilter(imageWithNoise);
	double maxDiff = 0;
	double maxDiffAfterMedianFilter = 0;
	double maxDiffAfterGeometric = 0;


	for (size_t c = 0; c < imageWithoutNoise.spectrum(); c++)
	{
		for (size_t x = 0; x < imageWithoutNoise.width(); x++)
		{
			for (size_t y = 0; y < imageWithoutNoise.height(); y++)
			{
				if (maxDiff < abs ( (imageWithoutNoise(x, y, 0, c) - imageWithNoise(x, y, 0, c))))
					maxDiff = abs((imageWithoutNoise(x, y, 0, c) - imageWithNoise(x, y, 0, c)));

				if (maxDiffAfterMedianFilter < abs ( (imageWithoutNoise(x, y, 0, c) - imageAfterMedianFilter(x, y, 0, c))))
					maxDiffAfterMedianFilter = abs((imageWithoutNoise(x, y, 0, c) - imageAfterMedianFilter(x, y, 0, c)));

				if (maxDiffAfterGeometric < abs ((imageWithoutNoise(x, y, 0, c) - imageAfterGeometricFilter(x, y, 0, c))))
					maxDiffAfterGeometric = abs((imageWithoutNoise(x, y, 0, c) - imageAfterGeometricFilter(x, y, 0, c)));
			}
		}

	}

	cout << "Maximum difference" << endl;
	showErrorData(maxDiff, maxDiffAfterMedianFilter, maxDiffAfterGeometric);
}

#endif