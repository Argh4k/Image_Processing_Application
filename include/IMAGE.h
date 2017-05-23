

#ifndef IMAGE_H_
#define IMAGE_H_

#include "CImg.h"
#include <string>


using namespace cimg_library;

const int num_of_layers = 3;
const int MAX_PIXEL_VALUE = 255;

// FUNCTION CREATING LOOKUPTABLE
int * lookuptable(float level, float(*operation)(float, float));
// FUNCTIONS FOR LOOKUPTABLE

float brightlut(float value, float level = 0);
float contrastlut(float value, float level = 0);
float negativelut(float value, float level = 0);

// FUNCTION USING LOOKUPTABLE
void basicoperations(float level, CImg<float> & image, float(*operation)(float, float));
void Comparisions();

// FLIPS
void verticalFlip(CImg<float> & image);
void horizontalFlip(CImg<float> & image);
void diagonalFlip(CImg<float>& image);

// SHRINK AND ENLARGE
CImg<float> shrink(CImg<float> & image);
CImg<float> enlarge(CImg<float> & image);

// MEDIANFILTER AND MEDIANA
CImg<float> * medianfilter(CImg<float> & image);
float mediana(CImg<float> & image, int x, int y, int c);

//GEOMETRICMEAN FILTER AND GEOMETRICMEAN
float geometricmean(CImg<float> & image, int x, int y, int c);
CImg<float> * geometricfilter(CImg<float> & image);

// SAVING IMAGE FUNCTION
void SaveImage(CImg<float> & image);
void SaveImage2(CImg<float> & image);

CImg<float> Load_Image(const char* name);

// ERRORS CALCULATIONS
void Mean_square_error(CImg<float> & image_without_noise, CImg<float> & image_with_noise);
void Peak_mean_square_error(CImg<float> & image_without_noise, CImg<float> & image_with_noise);
void Signal_to_noise_error(CImg<float> & image_without_noise, CImg<float> & image_with_noise);
void Peak_signal_to_noise_error(CImg<float> & image_without_noise, CImg<float> & image_with_noise);
void Maximum_difference(CImg<float> & image_without_noise,CImg<float> & image_with_noise);
double Find_maximum_value(CImg<float> & image_without_noise);
void Show_error_data(double clean_picture_data, double median_filter_data, double gmean_filter_data);



 ////////TASK 2 FUNCTIONS
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//HISTOGRAM
int * createhistogramtable(CImg<float> & image, int channel = 0 );
void Createhistogramimage(CImg<float> & image, int channel = 0);

//H1
void UniformFinalProbabilityDensityFunction(CImg<float> & image);

//S1
CImg<float> * Low_pass_filter(CImg<float> & image, int mask_chosen = 1);
float lpFilter(CImg<float> &image, int x, int y, int c, int * values_mask, double dividor);

// AMEAN AS OPTIMALIZATION FOR MASK 1
CImg<float> * ameanfilter(CImg<float> & image);
int amean(CImg<float> & image, int x, int y, int c);



//C1 - C6
float Cmean(CImg<float> & image, int channel = 0 );
float Cvariance(CImg<float> & image, int channel = 0);
float Cstdev(CImg<float> & image, int channel = 0);
float Cvarcoi(CImg<float> & image, int channel = 0);
float Casyco(CImg<float> & image, int channel = 0);
float Cfasyco(CImg<float> & image, int channel = 0);
float Centropy(CImg<float> & image, int channel = 0);
float Cvarcoii(CImg<float> & image, int channel = 0);

// 01
CImg<float> * Roberts_operator(CImg<float> & image);
float roFilter(CImg<float> & image, int x, int y, int c);
float roFilter2(CImg<float> & image, int x, int y, int c);


// M

class StructuralElement;
CImg<float> * Dilation(CImg<float> & image, StructuralElement & Elementk);
void dilate(CImg<float> * image, int x, int y, StructuralElement & Element);
CImg<float> * Erosion(CImg<float> & image, StructuralElement & Element);
bool erosecheck(CImg<float> & image, int x, int y, StructuralElement & Element);
CImg<float> * HMT(CImg<float> & image, StructuralElement & Element);
bool HMTcheck(CImg<float> & image, int x, int y, StructuralElement & StEl);

void Difference(CImg<float> & image1, CImg<float> & image2);



//region growing,  _x _y are coordinates of seed point
CImg<float>* Region_growing(CImg<float>* source_image, int tolerancy, int _x, int _y); 
void Grow_region(int** seed_array, CImg<float>* image,CImg<float>* output,  int tolerancy, int x, int y);

class StructuralElement {
private:
	unsigned char tab[9];
public:
	StructuralElement(int n);
	int operator()(int x, int y) { return tab[y * 3 + x]; }
};

#endif // !IMAGE_H_

