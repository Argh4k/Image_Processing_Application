#include "IMAGE.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include "CImg.h"


using namespace cimg_library;
using namespace std;


const int ArrSize = 255;


// FUNCTION CREATING LUT DEPENDIND ON THE OPERATION
int * lookuptable(float level, float(*operation)(float, float))
{
	int * lut = new int[ArrSize+1];
	{
		for (int i = 0; i < ArrSize +1 ; i++)
		{
			lut[i] = operation(i, level);
		}

	}
	return lut;
}

// POSSIBLE OPERATIONS FOR LUT

float brightlut(float value, float level)
{
	if (value + level < 0) { return 0; }
	else if (value + level > 255) { return 255; }
	else return value + level;
}

float contrastlut(float value, float level)
{
	

	if (((value - 128)*level) + 128 > 255)
	{
		return 255;
	}
	else if (((value - 128)*level) + 128 < 0) {
		return 0;
	}
	else {
		return ((value - 128)*level) + 128;
	}

}

float negativelut(float value, float level )
{
	float newvalue = 255 - value;
	if (newvalue < 0) return 0;
	else
		return newvalue;
	//if (value == 255)
	//	return 0;
	//else return 255;
}

// ITERATES THROUGH EVERY PIXEL AND CREATES SUTIABLE LUT
void basicoperations(float level, CImg<float> & image, float(*operation)(float, float)) {

	int * lut = lookuptable(level, operation);


	for (int x = 0; x < image.width(); x++)
	{
		for (int y = 0; y < image.height(); y++)
		{
			for (int c = 0; c < num_of_layers; c++)
			{
				image(x, y, 0, c) = lut[(int)image(x, y, 0, c)];
			}
		}
	}
	delete[] lut;
}

//FLIPS

void verticalFlip(CImg<float> & image) {

	for (int x = 0; x < image.width(); x++)
	{
		for (int y = 0; y < image.height() / 2; y++)
		{
			for (int c = 0; c < num_of_layers; c++) {
				swap(image(x, y, 0, c), image(x, image.height() - y - 1, 0, c));
			}
		}
	}

}

void horizontalFlip(CImg<float> & image) {

	for (int y = 0; y < image.height(); y++)
	{
		for (int x = 0; x < image.width() / 2; x++)
		{
			for (int c = 0; c < num_of_layers; c++) {
				swap(image(x, y, 0, c), image(image.width() - x - 1, y, 0, c));
			}
		}
	}

}

void diagonalFlip(CImg<float>& image) {

	for (int x = 0; x < image.width(); x++)
	{
		for (int y = 0; y < image.height() / 2; y++)
		{
			for (int c = 0; c < num_of_layers; c++)
			{
				swap(image(x, y, 0, c), image(image.width() - 1 - x, image.height() - 1 - y, 0, c));
			}
		}
	}
}

CImg<float> shrink(CImg<float> & image) {

	CImg<float> shrinked_image((image.width() / 2), (image.height() / 2), 1, 3);

	for (int x = 0; x < shrinked_image.width(); x++)  //to not run out of table
	{
		for (int y = 0; y < shrinked_image.height(); y++)
		{
			for (int c = 0; c < num_of_layers; c++)
			{

				shrinked_image(x, y, 0, c) = image(2 * x, 2 * y, 0, c); //every second pixel (0,0)(2,2)(4,4) etc

			}
		}
	}

	return shrinked_image;
}

CImg<float> enlarge(CImg<float> & image) {

	CImg<float> enlarged_image((image.width() * 2), (image.height() * 2), 1, 3);

	for (int x = 0; x < enlarged_image.width() - 1; x += 2)
	{
		for (int y = 0; y < enlarged_image.height() - 1; y += 2)
		{
			for (int c = 0; c < num_of_layers; c++)
			{

				enlarged_image(x, y, 0, c) = image(x / 2, y / 2, 0, c);
				enlarged_image(x + 1, y + 1, 0, c) = image(x / 2, y / 2, 0, c);

			}
		}
	}
	return enlarged_image;
}

float mediana(CImg<float> & image, int x, int y, int c)
{
	int * arr = new int[9];
	int temp = 0;
	for (int i = -1; i < 2; i++)
	{
		for (int j = -1; j < 2; j++)
		{
			arr[temp] = image(x + i, y + j, 0, c);
			temp++;
		}
	}
	std::sort(arr, arr + 9);
	int mediana = arr[4];
	delete[] arr;
	return mediana;
}

float geometricmean(CImg<float> & image, int x, int y, int c)
{
	double sum = 1;

	for (int i = -1; i < 2; i++)
	{
		for (int j = -1; j < 2; j++)
		{
			sum*=image(x + i, y + j, 0, c);
		}
	}
	double root = 1.0f / 9.0f;
	return pow(sum, root);
}

CImg<float> * medianfilter(CImg<float> & image) {
	CImg<float> * filterimage = new CImg<float>;
	*filterimage = image;
	for (int x = 1; x < image.width() - 1; x++)
	{
		for (int y = 1; y < image.height() - 1; y++)
		{
			for (int c = 0; c < num_of_layers; c++)
			{
				(*filterimage)(x, y, 0, c) = mediana(image, x, y, c);
			}
		}
	}
	return filterimage;
}

CImg<float> * geometricfilter(CImg<float> & image) {
	CImg<float> * filterimage = new CImg<float>;
	*filterimage = image;
	for (int x = 1; x < image.width() - 1; x++)
	{
		for (int y = 1; y < image.height() - 1; y++)
		{
			for (int c = 0; c < num_of_layers; c++)
			{
				(*filterimage)(x, y, 0, c) = geometricmean(image, x, y, c);
			}
		}
	}
	return filterimage;
}

void SaveImage(CImg<float> & image) {
	string name;
	name += "output.bmp";
	image.save(name.c_str());
}

void Mean_square_error(CImg<float> & image_without_noise, CImg<float> & image_with_noise)
{
	CImg<float> image_median = *medianfilter(image_with_noise);
	CImg<float> image_geometric = *geometricfilter(image_with_noise);
	double mean_error = 0;
	double mean_error_median = 0;
	double mean_error_geometric = 0;
	
	for (int c = 0; c < num_of_layers; c++)
	{
		for (int x = 0; x < image_without_noise.width(); x++)
		{
			for (int y = 0; y < image_without_noise.height(); y++)
			{
				mean_error += pow((image_without_noise(x, y, 0, c) - image_with_noise(x, y, 0, c)), 2);
				mean_error_geometric += pow((image_without_noise(x, y, 0, c) - image_geometric(x, y, 0, c)), 2);
				mean_error_median += pow((image_without_noise(x, y, 0, c) - image_median(x, y, 0, c)), 2);
			}
		}
		
	}


	double coefficient = (1.0 / (image_without_noise.width()*image_without_noise.height() * num_of_layers));
	mean_error = coefficient*mean_error;
	mean_error_geometric = coefficient*mean_error_geometric;
	mean_error_median = coefficient*mean_error_median;
	
	cout << "\nMean square error" << endl;
	Show_error_data(mean_error, mean_error_median, mean_error_geometric);

}

void Peak_mean_square_error(CImg<float> & image_without_noise, CImg<float> & image_with_noise)
{
	CImg<float> image_median = *medianfilter(image_with_noise);
	CImg<float> image_geometric = *geometricfilter(image_with_noise);
	double mean_error = 0;
	double mean_error_median = 0;
	double mean_error_geometric = 0;
	
	for (int c = 0; c < num_of_layers; c++)
	{
		for (int x = 0; x < image_without_noise.width(); x++)
		{
			for (int y = 0; y < image_without_noise.height(); y++)
			{
				mean_error += pow((image_without_noise(x, y, 0, c) - image_with_noise(x, y, 0, c)), 2);
				mean_error_geometric += pow((image_without_noise(x, y, 0, c) - image_geometric(x, y, 0, c)), 2);
				mean_error_median += pow((image_without_noise(x, y, 0, c) - image_median(x, y, 0, c)), 2);
			}
		}

	}
	int pow_max_value = pow(Find_maximum_value(image_without_noise), 2);
	double coefficient = (1.0 / (image_without_noise.width()*image_without_noise.height() * num_of_layers));
	mean_error = (coefficient*mean_error) / pow_max_value;
	mean_error_geometric = (coefficient*mean_error_geometric) / pow_max_value;
	mean_error_median = (coefficient*mean_error_median) / pow_max_value;

	cout << "\nPeak mean square error" << endl;
	Show_error_data(mean_error,mean_error_median,mean_error_geometric);

}

void Signal_to_noise_error(CImg<float> & image_without_noise, CImg<float> & image_with_noise)
{
	CImg<float> image_median = *medianfilter(image_with_noise);
	CImg<float> image_geometric = *geometricfilter(image_with_noise);
	double error = 0;
	double error_median = 0;
	double error_geometric = 0;
	double sum = 0;

	for (int c = 0; c < num_of_layers; c++)
	{
		for (int x = 0; x < image_without_noise.width(); x++)
		{
			for (int y = 0; y < image_without_noise.height(); y++)
			{
				sum += pow(image_without_noise(x, y, 0, c), 2);
				error += pow((image_without_noise(x, y, 0, c) - image_with_noise(x, y, 0, c)), 2);
				error_geometric += pow((image_without_noise(x, y, 0, c) - image_geometric(x, y, 0, c)), 2);
				error_median += pow((image_without_noise(x, y, 0, c) - image_median(x, y, 0, c)), 2);
			}
		}

	}
	
	error = 10*log10(sum/error);
	error_geometric = 10*log10(sum/ error_geometric);
	error_median = 10*log10(sum / error_median);

	cout << "Signal to noise error" << endl;
	Show_error_data(error,error_median,error_geometric);

}


void Peak_signal_to_noise_error(CImg<float> & image_without_noise, CImg<float> & image_with_noise)
{
	CImg<float> image_median = *medianfilter(image_with_noise);
	CImg<float> image_geometric = *geometricfilter(image_with_noise);
	double error = 0;
	double error_median = 0;
	double error_geometric = 0;
	double sum = 0;

	for (int c = 0; c < num_of_layers; c++)
	{
		for (int x = 0; x < image_without_noise.width(); x++)
		{
			for (int y = 0; y < image_without_noise.height(); y++)
			{
				
				error += pow((image_without_noise(x, y, 0, c) - image_with_noise(x, y, 0, c)), 2);
				error_geometric += pow((image_without_noise(x, y, 0, c) - image_geometric(x, y, 0, c)), 2);
				error_median += pow((image_without_noise(x, y, 0, c) - image_median(x, y, 0, c)), 2);
			}
		}
	}
	double coefficient = (1.0 / (image_without_noise.width()*image_without_noise.height() * num_of_layers));
	int pow_max_value = pow(Find_maximum_value(image_without_noise), 2);
	error = 10 * log10(pow_max_value / (error*coefficient));
	error_geometric = 10 * log10(pow_max_value / (error_geometric*coefficient));
	error_median = 10 * log10(pow_max_value / (error_median*coefficient));

	cout << "Peak signal to noise error" << endl;
	Show_error_data(error, error_median, error_geometric);
}

void Maximum_difference(CImg<float> & image_without_noise, CImg<float> & image_with_noise)
{
	CImg<float> image_median = *medianfilter(image_with_noise);
	CImg<float> image_geometric = *geometricfilter(image_with_noise);
	double max_diff_clean = 0;
	double max_diff_median = 0;
	double max_diff_geometric = 0;
	

	for (int c = 0; c < num_of_layers; c++)
	{
		for (int x = 0; x < image_without_noise.width(); x++)
		{
			for (int y = 0; y < image_without_noise.height(); y++)
			{
				if (max_diff_clean < abs ( (image_without_noise(x, y, 0, c) - image_with_noise(x, y, 0, c))))
					max_diff_clean = abs((image_without_noise(x, y, 0, c) - image_with_noise(x, y, 0, c)));

				if (max_diff_median < abs ( (image_without_noise(x, y, 0, c) - image_median(x, y, 0, c))))
					max_diff_median = abs((image_without_noise(x, y, 0, c) - image_median(x, y, 0, c)));

				if (max_diff_geometric < abs ((image_without_noise(x, y, 0, c) - image_geometric(x, y, 0, c))))
					max_diff_geometric = abs((image_without_noise(x, y, 0, c) - image_geometric(x, y, 0, c)));
			}
		}

	}

	

	cout << "Maximum difference" << endl;
	Show_error_data(max_diff_clean, max_diff_median, max_diff_geometric);
}
double Find_maximum_value(CImg<float> & image) {

	double maximum_value = 0;
	
	for (int x = 0; x  < image.width(); x++)
	{
		for (int y = 0; y < image.height(); y++)
		{
			for (int c = 0; c < num_of_layers; c++)
			{
				if (maximum_value < image(x, y, 0, c))
				{
					if (image(x, y, 0, c) == MAX_PIXEL_VALUE)
						return	MAX_PIXEL_VALUE;
					else
					maximum_value = image(x, y, 0, c);
				}
			}
		}
	}
	return maximum_value;
}

void Show_error_data(double clean_picture_data, double median_filter_data, double gmean_filter_data) {

	cout << "\nBefore processing " << clean_picture_data << endl;
	cout << "After median filter " << median_filter_data << endl;
	cout << "After geometric filter " << gmean_filter_data << endl;
}
