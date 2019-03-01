/*
#include "IMAGE.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include "CImg.h"
#include <sstream>


using namespace cimg_library;
using namespace std;

const int HISTOGRAM_SIZE = 256; //max possible level of intensity

class ImageMask {
private:
	float dividor;
	int * mask_values;
	string mask_1 = "1 1 1 1 1 1 1 1 1";
	string mask_2 = "1 1 1 1 2 1 1 1 1";
	string mask_3 = "1 2 1 2 4 2 1 2 1";
public:
	float Get_dividor() { return dividor; }
	int *Get_mask() { return mask_values; }

	ImageMask(int i) {
		stringstream ss;
		mask_values = new int[9];
		if (i == 1)
		{
			 dividor = 9.0;
			 ss << mask_1;
		}
		else if (i == 2)
		{
			dividor = 10.0;
			ss << mask_2;
		}
		else if (i == 3)
		{
			dividor = 16.0;
			ss << mask_3;
		}

		for (int i = 0; i < 9; i++)
		{
			ss >> mask_values[i];
		}
	}
	~ImageMask() {
		delete[] mask_values;
}
};


CImg<float>  Load_Image(const char* name) {

	ifstream file(name);

	if(file.good())
	{
		CImg<float> image;
		image.load(name);
		return image;
	}
	else
	{
		cout << " Couldn't find file: " << name << endl;
		exit(0);
	}
}


int * createhistogramtable(CImg<float> & image, int channel)
{
	int * histogram = new int[HISTOGRAM_SIZE]();
	for (int x = 0; x < image.width(); x++)
	{
		for (int y = 0; y < image.height(); y++)
		{
			histogram[(int)image(x, y,0, channel)]++;
		}
	}
	cout << endl << endl;
	return histogram;
}


void Createhistogramimage(CImg<float> & image,int channel)
{
	int * histogram = createhistogramtable(image, channel);
	int histogramheight = 0;
	for (int i = 0; i < HISTOGRAM_SIZE; i++)
	{
		if (histogram[i] > histogramheight)
			histogramheight = histogram[i];
	}
	CImg<float> histogrampicture(HISTOGRAM_SIZE, histogramheight/10 + 1, 1, 3);
	for (int x = 0; x < HISTOGRAM_SIZE; x+=2)
	{
		for (float y = (histogramheight)/10 ; y > histogramheight/10-histogram[x]/10 ; y--)
		{
			histogrampicture(x, y, 0, channel) = 255;
			histogrampicture(x+1, y, 0, channel) = 255;
		}
	}
	SaveImage(histogrampicture);
	delete[] histogram;
}


void UniformFinalProbabilityDensityFunction(CImg<float> & image)
{
	double numofpixels;
	int gmax, gmin, sum;
	gmax = 255;
	gmin = 0;
	numofpixels = image.width()*image.height();
	int * histogram;
	float * lut = new float[HISTOGRAM_SIZE];
	for (int c = 0; c < image.spectrum(); c++)
	{
		histogram = createhistogramtable(image, c);
		for (int i = 0; i < HISTOGRAM_SIZE; i++)
		{
			sum = 0;
			for (int m = 0; m <= i; m++)
			{
				sum += histogram[m];
			}
			lut[i] = gmin + (gmax - gmin) * (1.0 / numofpixels) * sum;
		}

		for (int x = 0; x < image.width(); x++)
		{
			for (int y = 0; y < image.height(); y++)
			{
					image(x, y, 0 , c) = lut[(int)image(x, y, 0 , c)];
			}
		}
		delete[] histogram;
	}
	SaveImage(image); 
	delete[] lut;

}




CImg<float> * Low_pass_filter(CImg<float> & image, int mask_chosen) {
	CImg <float> * filtredimage = new CImg<float>(image);

	ImageMask mask(mask_chosen);
	int * values_mask = mask.Get_mask();

	double dividor = mask.Get_dividor();

	for (int x = 1; x < image.width() - 1; x++)			//in such case we avoid the border pixels
	{
		for (int y = 1; y < image.height() - 1; y++)
		{
			for (int c = 0; c < image.spectrum(); c++) 
			{
				
				(*filtredimage)(x, y, 0, c) = lpFilter(image, x, y, c, values_mask, dividor);
			}
		}
	}
	return filtredimage;
}


//mask_chose is argv[3]!  btw fajny bajer, ten komentarz pojawia siê przy funkcji :)
float lpFilter(CImg<float> &image, int x, int y, int c, int values_mask[], double dividor ) {

	int iterator = 0;
	double sum = 1.0;
	
	for (int _y = y - 1; _y < y + 2; _y++)
	{

		for (int _x = x - 1; _x < x + 2; _x++)
		{
			sum += (values_mask[iterator] * image(_x, _y, 0, c));
			iterator++;
		}
	}


	return (sum/dividor);  
	
}


CImg<float> * ameanfilter(CImg<float> & image) {
	CImg <float> * filtredimage = new CImg<float>(image);

	for (int x = 1; x < image.width() - 1; x++)			//in such case we avoid the border pixels
	{
		for (int y = 1; y < image.height() - 1; y++)
		{
			for (int c = 0; c < image.spectrum(); c++)
			{

				(*filtredimage)(x, y, 0, c) = amean(image, x, y, c);
			}
		}
	}
	return filtredimage;
}

int amean(CImg<float> & image, int x, int y, int c)
{
	int sum = 0;
	for (int _y = y - 1; _y < y + 2; _y++)
	{

		for (int _x = x - 1; _x < x + 2; _x++)
		{

			sum += image(_x, _y, 0, c);
		}
	}
	return sum / 9;
}

float Cmean(CImg<float> & image,int  channel ) {
	float number_of_pixels = image.height()*image.width();
	int* histogram;
	float sum = 0;

		histogram = createhistogramtable(image, channel);
		for (int i = 0; i < HISTOGRAM_SIZE - 1; i++)
		{
			sum += i*histogram[i];
		}
		delete[] histogram;
	
		//cout << "\n\n channel:  " << channel << endl;
	
	return (1.0 /number_of_pixels)*sum;
}

float Cvariance(CImg<float> & image, int channel)
{
	float number_of_pixels = image.height()*image.width();
	float mean = Cmean(image, channel);

	int* histogram;
	float sum = 0;

	histogram = createhistogramtable(image,channel);

	for (int i = 0; i < HISTOGRAM_SIZE - 1; i++)
	{
		sum += pow((i - mean),2) * histogram[i];
	}
	delete[] histogram;

	return (1.0/number_of_pixels)*sum;
}

float Cstdev(CImg<float> & image, int channel) {
	return pow(Cvariance(image, channel), 0.5);
}

float Cvarcoi(CImg<float> &image, int channel)
{
	return Cstdev(image, channel) / Cmean(image, channel);
}

float Casyco(CImg<float>& image, int channel)
{
	float number_of_pixels = image.height()*image.width();
	float sigma = Cstdev(image, channel);
	int* histogram;
	float sum = 0;
	float mean = Cmean(image,channel);

	histogram = createhistogramtable(image, channel);

	for (int i = 0; i < HISTOGRAM_SIZE - 1; i++)
	{
		sum += pow((i - mean), 3) * histogram[i];
	}
	delete[] histogram;


	return ((1 / pow(sigma, 3))*(1/number_of_pixels)*sum);
}

float Cfasyco(CImg<float> & image, int channel) {

	float number_of_pixels = image.height()*image.width();
	float sigma = Cstdev(image,channel);
	int* histogram;
	float sum = 0;
	float mean = Cmean(image, channel);

	histogram = createhistogramtable(image, channel);

	for (int i = 0; i < HISTOGRAM_SIZE - 1; i++)
	{
		sum += pow((i - mean), 4) * histogram[i];
	}
	delete[] histogram;


	return ((1 / pow(sigma, 4))*(1 / number_of_pixels)*sum)-3;
}

float Cvarcoii(CImg<float> & image, int channel) {

	float number_of_pixels = image.height()*image.width();
	int* histogram;
	float sum = 0;
	float mean = Cmean(image, channel);

	histogram = createhistogramtable(image, channel);

	for (int i = 0; i < HISTOGRAM_SIZE - 1; i++)
	{
		sum += pow( histogram[i],2);
	}
	delete[] histogram;


	return pow((1 / number_of_pixels),2)*sum;
}

float Centropy(CImg<float> & image, int channel) {

	float number_of_pixels = image.height()*image.width();
	int* histogram;
	float sum = 0;
	float mean = Cmean(image, channel);
	double log;

	histogram = createhistogramtable(image, channel);

	for (int i = 0; i < HISTOGRAM_SIZE - 1; i++)
	{
		if (histogram[i] == 0)
			log = 0;
		else log = log2(histogram[i] / number_of_pixels);
		sum += histogram[i]*log;
	}
	delete[] histogram;
	return (-1/number_of_pixels)*sum;

}

CImg<float> * Roberts_operator(CImg<float> & image)
{
	CImg<float> * filtredimage = new CImg<float>(image);
	for (int x = 1; x < image.width() - 1; x++)			//in such case we avoid the border pixels
	{
		for (int y = 1; y < image.height() - 1; y++)
		{
			for (int c = 0; c < image.spectrum(); c++)
			{
				
				(*filtredimage)(x, y, 0, c) = roFilter(image, x, y, c);
				if ((*filtredimage)(x, y, 0, c)>255)
					(*filtredimage)(x, y, 0, c) = 255;
			}
		}
	}
	return filtredimage;
}

float roFilter(CImg<float> & image, int x, int y, int c)
{
	float pow1 = pow((image(x, y, 0, c) - image(x + 1, y + 1, 0 , c)), 2);
	float pow2 = pow((image(x, y + 1, 0, c) - image(x + 1, y, 0, c)), 2);
	return pow((pow1 + pow2), 0.5);
}


float roFilter2(CImg<float> & image, int x, int y, int c)
{
	float abs1 = abs(image(x, y, 0, c) - image(x + 1, y + 1, 0, c));
	float abs2 = abs(image(x, y + 1, 0, c) - image(x + 1, y, 0, c));
	return abs1 + abs2;
}

*/