#define _USE_MATH_DEFINES
#include "Fourier.h"
#include <cmath>
#include <iostream>

using namespace std;

const complex<double> i(0.0, 1.0);
double pi = M_PI;


complex<double> ** DFT(CImg<float> & image)
{
	complex<double> ** Arr = new complex<double>*[image.height()];
	for (int i = 0; i < image.height(); ++i)
	{
		Arr[i] = new complex<double>[image.width()];
	}

	complex<double> ** Arr2 = new complex<double>*[image.height()];
	for (int i = 0; i < image.height(); ++i)
	{
		Arr2[i] = new complex<double>[image.width()];
	}

	for (int i = 0; i < image.height(); i++)
	{

		for (int x = 0; x < image.width(); x++)
		{
			Arr2[i][x] = image(x, i);
		}
	}


	for (int i = 0; i < image.height(); i++)
	{

		for (int x = 0; x < image.width(); x++)
		{
			Arr[i][x] = First_Transform(Arr2, x, i, image.width(), false);
			
		}
	}


	for (int i = 0; i < image.width(); i++)
	{
		for (int x = 0; x < image.height(); x++)
		{
			Arr2[x][i] = Second_Transform(Arr, i, x, image.height(), false);
		}
	}



	for (int i = 0; i < image.height(); ++i) {
		delete[] Arr[i];
	}
	delete[] Arr;
	
	return Arr2;


}

complex<double> ** IDFT(complex<double> ** Arr, int M, int N)
{

	complex<double> ** Arr1 = new complex<double>*[N];
	for (int i = 0; i < N; i++)
	{
		Arr1[i] = new complex<double>[M];
	}

	complex<double> ** Arr2 = new complex<double>*[N];
	for (int i = 0; i < N; i++)
	{
		Arr2[i] = new complex<double>[M];
	}

	for (int i = 0; i < M; i++)
	{
		for (int x = 0; x < N; x++)
		{
	
			Arr1[i][x] = First_Transform(Arr, x, i, M, true);
		}
	}


	for (int i = 0; i < N; i++)
	{

		for (int x = 0; x < M; x++)
		{
			Arr2[x][i] = Second_Transform(Arr1,i, x, N, true);;
		}
	}

	for (int i = 0; i < N; ++i) {
		delete[] Arr1[i];
	}
	delete[] Arr1;

	return Arr2;
}

complex<double> First_Transform(complex<double> ** Arr, int column, int row, double N, bool inverse)
{
	complex<double> sum(0.0, 0.0);
	double k;
	if (inverse)
		k = -1;
	else
		k = 1;
	double alfa;
	for (int x = 0; x < N; ++x)
	{
		alfa = ((2 * pi* column * x) / N);
		sum += Arr[row][x] * (cos(alfa) - k*i*sin(alfa));
	}
	if(inverse)
		sum = sum / N;
	return sum;
}

complex<double> Second_Transform(complex<double> ** Arr, int column, int row, double N, bool inverse)
{
	complex<double> sum(0.0, 0.0);
	double alfa;
	double k;
	if (inverse)
		k = -1;
	else
		k = 1;
	for (int x = 0; x < N; ++x)
	{
		alfa = ((2 * pi* row * x) / N);
		sum += Arr[x][column] * (cos(alfa) - k*i*sin(alfa));
	}
	if (inverse)
		sum = sum / N;
	return sum;

}

void imageswap(CImg<float> & image)
{
	for (int i = 0; i < image.width() / 2; i++)
	{
		for (int j = 0; j < image.height() / 2; j++)
		{
			std::swap(image(i, j), image(i + image.width() / 2, j + image.height() / 2));
		}
	}
	for (int i = image.width() / 2; i < image.width(); i++)
	{
		for (int j =0; j < image.height()/2; j++) 
		{
			std::swap(image(i, j), image(i - image.width() / 2, j + image.height() / 2));
		}
	}
}

void savefourier(complex<double> ** Arr, int width, int height)
{
	CImg<float> * tosave = new CImg<float>(width, height);
	

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			(*tosave)(x, y) = abs(Arr[y][x]);
		}
	}
	imageswap((*tosave));
	tosave->save("test.bmp");
}

complex<double> ** FFT(CImg<float> & image)
{
	complex<double> ** Arr = new complex<double>*[image.height()];
	for (int i = 0; i < image.height(); ++i)
	{
		Arr[i] = new complex<double>[image.width()];
	}


	for (int i = 0; i < image.height(); i++)
	{

		for (int x = 0; x < image.width(); x++)
		{
			Arr[i][x] = image(x, i);
		}
	}

	//Transform each row

	for (int i = 0; i < image.height(); i++)
	{
		FFT_ROW(Arr[i], image.width());
	}


	//Transform each colum
	
	for (int i = 0; i < image.width(); i++)
	{
		FFT_COLUMN(Arr, image.height(), i);
	}

	for (int i = 0; i < (image.width() / 2); ++i)
	{
		for (int j = 0; j < (image.height() / 2); ++j)
		{
			std::swap(Arr[j][i], Arr[j + image.height() / 2][i + image.width() / 2]);
		}
	}
	for (int i = image.width() / 2; i < image.width(); ++i)
	{
		for (int j = 0; j < (image.height() / 2); ++j)
		{
			std::swap(Arr[j][i], Arr[j + image.height() / 2][i - image.width() / 2]);
		}
	}

	return Arr;
}

complex<double> ** IFFT(complex<double> ** Arr, int M, int N) {
	

	
	conjugate(Arr, M, N);


	for (int i = 0; i < M; i++)
	{
		FFT_COLUMN(Arr, N, i);
	}
	


	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			Arr[j][i] /= M;

	//conjugate(Arr, M, N);
	//Transform each colum

	for (int i = 0; i < N; i++)
	{
		FFT_ROW(Arr[i], M);
	}

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			Arr[j][i] /=  N;
	
	return Arr;
}

void FFT_ROW(complex<double> Arr[], int length)
{
	if (length <= 1)
		return;

	complex<double> * ArrOdd = new complex<double>[length / 2];
	complex<double> * ArrEven = new complex<double>[length / 2];
	for (int i = 0, x = 0; i < length / 2; x = x + 2, i++)
	{
		ArrEven[i] = Arr[x];
		ArrOdd[i] = Arr[x + 1];
	}


	FFT_ROW(ArrEven, length / 2);
	FFT_ROW(ArrOdd, length / 2);

	
	complex<double> alfa;
	complex<double> t;
	for (int k = 0; k < length / 2; ++k)
	{
		alfa = ((2 * pi * k) / length);
		t = (cos(alfa) - i*sin(alfa)) * ArrOdd[k];
		Arr[k] = ArrEven[k] + t;
		Arr[k + length / 2] = ArrEven[k] - t;
	}
}

void FFT_COLUMN(complex<double> ** Arr, int length, int column)
{

	complex<double> * Array = new complex<double>[length];
	for (int i = 0; i < length; i++)
	{
		Array[i] = Arr[i][column];
	}


	FFT_ROW(Array, length);

	for (int i = 0; i < length; i++)
	{
		Arr[i][column] = Array[i];
	}
	delete[] Array;
}

complex<double>** LowPassFilter(CImg<float> &image,int radius) {

	//little idiotproof solution :P
	if (radius > image.width() / 2 || radius > image.height() / 2)
	{
		cout << "You have chosen too big radius" << endl;
		return 0;
	}

	complex<double>  **mask = FFT(image);
	

	//just 4 debbuging
	CImg<float> mask_image(image.width(), image.height());
		mask_image = *PrintMask(mask,image.width(),image.height());
		mask_image.save("mask.bmp");

		for (int x = 0; x < image.width(); ++x)
		{
			for (int y = 0; y < image.height(); ++y)
			{
				if (!checkRadius(x, y, image.width() / 2, image.height() / 2, radius))
					mask[y][x] = 0;
			}
		}

		mask_image = *PrintMask(mask, image.width(), image.height());
		mask_image.save("mask2.bmp");

		return mask;
}

CImg<float>* PrintMask(complex<double>**Arr, int N, int M)
{
	CImg<float> * toreturn = new CImg<float>(M, N);
	//double coefficient = 1 / sqrt(M*N);
	for (int x = 0; x < toreturn->width(); x++)
	{
		for (int y = 0; y < toreturn->height(); y++)
		{
			//(*toreturn)(x, y) = 10*log(1+abs(Arr[y][x]));
			(*toreturn)(x, y) = log(1+abs(Arr[y][x]));

		}
	}
	return toreturn;
}

CImg<float>* Print(complex<double>**Arr, int N, int M)
{
	CImg<float> * toreturn = new CImg<float>(M, N);
	for (int x = 0; x < toreturn->width(); x++)
	{
		for (int y = 0; y < toreturn->height(); y++)
		{
			(*toreturn)(x, y) = abs(Arr[y][x]);
		}
	}
	return toreturn;
}

bool checkRadius(int x, int y, int x_0, int y_0, int radius)
{
	if (pow(x - x_0, 2) + pow(y - y_0, 2) <= radius*radius)
		return true;
	else
		return false;
}

complex<double>** HighPassFilter(CImg<float> &image, int radius) {

	//little idiotproof solution :P
	if (radius > image.width() / 2 || radius > image.height() / 2)
	{
		cout << "You have chosen too big radius" << endl;
		return 0;
	}

	complex<double>  **mask = FFT(image);


	//just 4 debbuging
	CImg<float> mask_image(image.width(), image.height());
	mask_image = *PrintMask(mask, image.width(), image.height());
	mask_image.save("mask.bmp");

	for (int x = 0; x < image.width(); ++x)
	{
		for (int y = 0; y < image.height(); ++y)
		{
			if (checkRadius(x, y, image.width() / 2, image.height() / 2, radius))
				mask[y][x] = 0;
		}
	}

	mask_image = *PrintMask(mask, image.width(), image.height());
	mask_image.save("mask2.bmp");

	return mask;
}

complex<double>** BandPassFilter(CImg<float> &image, int radius_start, int radius_end) {

	//little idiotproof solution :P
	if (radius_end > image.width() / 2 || radius_end > image.height() / 2)
	{
		cout << "You have chosen too big radius" << endl;
		return 0;
	}

	complex<double>  **mask = FFT(image);


	//just 4 debbuging
	CImg<float> mask_image(image.width(), image.height());
	mask_image = *PrintMask(mask, image.width(), image.height());
	mask_image.save("mask.bmp");

	for (int x = 0; x < image.width(); ++x)
	{
		for (int y = 0; y < image.height(); ++y)
		{
			if (!checkRadiusRegion(x, y, image.width() / 2, image.height() / 2, radius_start,radius_end))
				mask[y][x] = 0;
		}
	}

	mask_image = *PrintMask(mask, image.width(), image.height());
	mask_image.save("mask2.bmp");

	return mask;
}

complex<double>** BandCutFilter(CImg<float> &image, int radius_start, int radius_end) {

	//little idiotproof solution :P
	if (radius_end > image.width() / 2 || radius_end > image.height() / 2)
	{
		cout << "You have chosen too big radius" << endl;
		return 0;
	}

	complex<double>  **mask = DFT(image);


	//just 4 debbuging
	CImg<float> mask_image(image.width(), image.height());
	mask_image = *PrintMask(mask, image.width(), image.height());
	mask_image.save("mask.bmp");

	for (int x = 0; x < image.width(); ++x)
	{
		for (int y = 0; y < image.height(); ++y)
		{
			if (checkRadiusRegion(x, y, image.width() / 2, image.height() / 2, radius_start, radius_end))
				mask[y][x] = 0;
		}
	}

	mask_image = *PrintMask(mask, image.width(), image.height());
	mask_image.save("mask2.bmp");

	return mask;
}

bool checkRadiusRegion(int x, int y, int x_0, int y_0, int radius_start, int radius_end)
{
	if (pow(radius_start,2) <= pow(x - x_0, 2) + pow(y - y_0, 2) && pow(x - x_0, 2) + pow(y - y_0, 2) <= pow(radius_end,2))
		return true;
	else
		return false;
}

complex<double>** MaskFilter(int variant, CImg<float>& image)
{
	//TO DO:
	// flip mask for proper variant
	CImg<float> mask_img;
	switch (variant)
	{
	case 1:
		mask_img.load("F51.bmp");
		break;	
	case 2:
		mask_img.load("F52.bmp");
		break;
	case 3:
		mask_img.load("F53.bmp");
		break;
	case 4:
		mask_img.load("F54.bmp");
		break;
	case 5:
		mask_img.load("F55.bmp");
		break;
	}

	
	complex<double>  **mask = FFT(image);

	//just 4 debbuging
	CImg<float> mask_image(image.width(), image.height());
	mask_image = *PrintMask(mask, image.width(), image.height());
	mask_image.save("mask.bmp");

	for (int x = 0; x < image.width(); ++x)
	{
		for (int y = 0; y < image.height(); ++y)
		{
			if (mask_img(x,y) == 0)
				mask[y][x] = 0;
		}
	}

	mask_image = *PrintMask(mask, image.width(), image.height());
	mask_image.save("mask2.bmp");

	return mask;

}

complex<double>** PhaseMod(CImg<float>& image, int k, int l)
{
	int M = image.width();
	int N = image.height();

	complex<double> ** Arr = new complex<double>*[N];
	for (int i = 0; i < N; i++)
	{
		Arr[i] = new complex<double>[M];
	}

	complex<double> ** Arr2;
	
	
	for (int y = 0; y < N; ++y)
	{
		for (int x = 0; x < M; ++x)
		{
			Arr[y][x] = exp(i*((-x*k*2*pi/N)+(-y*l*2*pi/M)+(k+l)*pi));
		}
	}


	/*for (int i = 0; i < (image.width() / 2); ++i)
	{
		for (int j = 0; j < (image.height() / 2); ++j)
		{
			std::swap(Arr[j][i], Arr[j + image.height() / 2][i + image.width() / 2]);
		}
	}
	for (int i = image.width() / 2; i < image.width(); ++i)
	{
		for (int j = 0; j < (image.height() / 2); ++j)
		{
			std::swap(Arr[j][i], Arr[j + image.height() / 2][i - image.width() / 2]);
		}
	}*/


	Arr2 = FFT(image);

	CImg<float> mask_image(image.width(), image.height());
	mask_image = *PrintMask(Arr2, image.width(), image.height());
	mask_image.save("mask.bmp");

	for (int y = 0; y < N; ++y)
	{
		for (int x = 0; x < M; ++x)
		{
			Arr2[y][x] *= Arr[y][x];
		}
	}


	mask_image = *PrintMask(Arr2, image.width(), image.height());
	mask_image.save("mask2.bmp");

	for (int i = 0; i < N; ++i) {
		delete[] Arr[i];
	}
	delete[] Arr;

	return Arr2;

}

void conjugate(complex<double> ** Arr, int M, int N)
{
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			Arr[j][i] = std::conj(Arr[j][i]);
}

void rescale(CImg<float> & image)
{
	double coefficient = sqrt(image.width()*image.height());
	for (int x = 0; x < image.width(); x++)
	{
		for (int y = 0; y < image.height(); y++)
		{
			image(x, y) = coefficient*image(x,y);
		}
	}
}