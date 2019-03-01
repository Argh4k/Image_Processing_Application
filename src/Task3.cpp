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

const int FOREGROUND = 0;
const int BACKGROUND = 255;
const int INACTIVE = 50;

const int VISITED = 2;
const int SEED = 1;



StructuralElement::StructuralElement(int n)
{
	if (n == 3)
	{
		for (int i = 0; i < 9; i++)
		{
			tab[i] = FOREGROUND;
		}
	}
	if (n == 1)
	{
		for (int i = 0; i < 9; i++)
		{
			tab[i] = BACKGROUND;
		}
		tab[4] = FOREGROUND;
		tab[5] = FOREGROUND;
	}
	if( n == 2)
	{
		for (int i = 0; i < 9; i++)
		{
			tab[i] = BACKGROUND;
		}
		tab[4] = FOREGROUND;
		tab[7] = FOREGROUND;
	}
	if (n == 4)
	{
		for (int i = 0; i < 9; i++)
		{
			tab[i] = BACKGROUND;
		}
		tab[1] = FOREGROUND;
		tab[3] = FOREGROUND;
		tab[4] = FOREGROUND;
		tab[5] = FOREGROUND;
		tab[7] = FOREGROUND;
	}
	if (n == 111)
	{
		tab[0] = FOREGROUND;
		tab[1] = INACTIVE;
		tab[2] = INACTIVE;
		tab[3] = FOREGROUND;
		tab[4] = BACKGROUND;
		tab[5] = INACTIVE;
		tab[6] = FOREGROUND;
		tab[7] = INACTIVE;
		tab[8] = INACTIVE;
	}
	if (n == 121)
	{
		tab[0] = BACKGROUND;
		tab[1] = BACKGROUND;
		tab[2] = BACKGROUND;
		tab[3] = INACTIVE;
		tab[4] = FOREGROUND;
		tab[5] = INACTIVE;
		tab[6] = FOREGROUND;
		tab[7] = FOREGROUND;
		tab[8] = FOREGROUND;
	}
	if (n == 111)
	{
		tab[0] = FOREGROUND;
		tab[1] = INACTIVE;
		tab[2] = INACTIVE;
		tab[3] = FOREGROUND;
		tab[4] = BACKGROUND;
		tab[5] = INACTIVE;
		tab[6] = FOREGROUND;
		tab[7] = INACTIVE;
		tab[8] = INACTIVE;
	}
	if (n == 7)
	{

		for (int i = 0; i < 9; i++)
		{
			tab[i] = BACKGROUND;
		}
		tab[3] = FOREGROUND;
		tab[4] = FOREGROUND;
		tab[5] = FOREGROUND;
	}
	if (n == 8)
	{
		tab[0] = FOREGROUND;
		tab[1] = FOREGROUND;
		tab[2] = BACKGROUND;
		tab[3] = BACKGROUND;
		tab[4] = FOREGROUND;
		tab[5] = BACKGROUND;
		tab[6] = BACKGROUND;
		tab[7] = FOREGROUND;
		tab[8] = BACKGROUND;
	}
}



CImg<float> * Dilation(CImg<float> & image, StructuralElement & Element) {
	CImg <float> * filtredimage = new CImg<float>(image);

	
	for (int x = 1; x < image.width()-1; x++)			//in such case we avoid the border pixels
	{
		for (int y = 1; y < image.height()-1; y++)
		{
			if (image(x, y, 0, 0) == FOREGROUND)
			{
				dilate(filtredimage, x, y, Element);
			}
		}
	}
	return filtredimage;
}


void dilate(CImg<float> * image, int x, int y, StructuralElement & Element)
{
	for (int _x = 0; _x < 3; _x++)
	{
		for (int _y = 0; _y < 3; _y++)
		{
			if (Element(_x, _y) == FOREGROUND)
			{
				(*image)(x - 1 + _x, y - 1 + _y) = FOREGROUND;
			}
		}
	}
	// If the element in the middle is BACKGROUND (VI, VIII)
	if (Element(1, 1) == BACKGROUND)
	{
		(*image)(x, y) = BACKGROUND;
	}
}

CImg<float> * Erosion(CImg<float> & image, StructuralElement & Element) {
	CImg <float> * filtredimage = new CImg<float>(image);


	for (int x = 1; x < image.width() - 1; x++)			//in such case we avoid the border pixels
	{
		for (int y = 1; y < image.height() - 1; y++)
		{
			if (image(x, y, 0, 0) == FOREGROUND)
			{
				if (!erosecheck(image, x, y, Element))
					(*filtredimage)(x, y, 0, 0) = BACKGROUND;
			}
		}
	}
	return filtredimage;
}


bool erosecheck(CImg<float> & image, int x, int y, StructuralElement & Element)
{
	for (int _x = 0; _x < 3; _x++)
	{
		for (int _y = 0; _y < 3; _y++)
		{
			if (Element(_x, _y) == FOREGROUND)
			{
				if (!(image(x - 1 + _x, y - 1 + _y) == FOREGROUND))
				{
					return false;
				}
			}
		}
	}
	return true;
}


CImg<float> * HMT(CImg<float> & image, StructuralElement & Element)
{
	CImg<float> * filtredimage = new CImg<float>(image);
	for (int x = 1; x < image.width()-1; x++)
	{
		for (int y = 1; y < image.height()-1; y++)
		{
			HMTcheck(image, x, y, Element) ? (*filtredimage)(x, y) = FOREGROUND : (*filtredimage)(x, y) = BACKGROUND;
		}
	}
	return filtredimage;
}


bool HMTcheck(CImg<float> & image, int x, int y, StructuralElement & StEl)
{
	for (int _x = 0; _x < 3; _x++)
	{
		for (int _y = 0; _y < 3; _y++)
		{
			if (!(StEl(_x, _y) == INACTIVE))
			{
				if (!(StEl(_x, _y) == image(x - 1 + _x, y - 1 + _y)))
				{
					return false;
				}
			}
		}
	}
	return true;
}


void Difference(CImg<float> & image1, CImg<float> & image2)
{
	//CImg<float> * newimg = new CImg<float>(image1);
	for (int x = 0; x < image1.width(); x++)
	{
		for (int y = 0; y < image1.width(); y++)
		{
			if (image1(x, y) == image2(x, y))
				image1(x, y) = BACKGROUND;
		}
	}
	//return image1;
}


unsigned int counter = 0;

CImg<float>* Region_growing(CImg<float>* source_image, int tolerancy, int _x, int _y) 
{

	CImg <float> *output = new CImg<float>((*source_image).width(), (*source_image).height());


	int** seed_array = new int*[(*source_image).width()];

	for (int i = 0; i < (*source_image).width(); i++)
		seed_array[i] = new int[(*source_image).height()];


	float treshold = (*source_image)(_x,_y);

	//cout << "he treshold is: "<< treshold << endl;

	for (int x = 0; x < (*source_image).width(); x++)
	{
		for (int y = 0; y < (*source_image).height(); y++)
		{
			if (abs((*source_image)(x, y) - treshold)<=tolerancy)
				seed_array[x][y] = SEED;
			else
				seed_array[x][y] = BACKGROUND;
		}

	}

		Grow_region(seed_array, source_image, output, tolerancy, _x, _y);

	for (int i = 0; i < (*source_image).width(); i++)
			delete seed_array[i];

	delete[] seed_array;
	return output;

}

void Grow_region(int** seed_array, CImg<float> *image, CImg<float>* output, int tolerancy, int x, int y)
{
	//cout << counter++ << endl;
		//cout << "test2" << endl;
			(*output)(x, y) = FOREGROUND;
			seed_array[x][y] = VISITED;
			

			if ( seed_array[x+1][y-1] == SEED && x > 1 && x < (*image).width() - 2 && y > 1 && y < (*image).height() - 1)	Grow_region(seed_array, image, output,  tolerancy, x + 1, y - 1);
			if (seed_array[x-1][y-1] == SEED && x > 1 && x < (*image).width() - 2 && y > 1 && y < (*image).height() - 1) 	Grow_region(seed_array, image, output,  tolerancy, x - 1, y - 1);
			if (seed_array[x+1][y+1] == SEED && x > 1 && x < (*image).width() - 2 && y > 1 && y < (*image).height() - 1) 	Grow_region(seed_array, image, output,  tolerancy, x + 1, y + 1);
			if (seed_array[x-1][y+1] == SEED && x > 1 && x < (*image).width() - 2 && y > 1 && y < (*image).height() - 1) 	Grow_region(seed_array, image, output,  tolerancy, x - 1, y + 1);
			if (seed_array[x][y+1] == SEED && x > 1 && x < (*image).width() - 2 && y > 1 && y < (*image).height() - 1)		Grow_region(seed_array, image, output,  tolerancy, x, y + 1);
			if (seed_array[x][y-1] == SEED && x > 1 && x < (*image).width() - 2 && y > 1 && y < (*image).height() - 1)		Grow_region(seed_array, image, output, tolerancy, x, y - 1);
			if (seed_array[x+1][y] == SEED && x > 1 && x < (*image).width() - 2 && y > 1 && y < (*image).height() - 1)		Grow_region(seed_array, image, output,  tolerancy, x + 1, y);
			if (seed_array[x-1][y] == SEED && x > 1 && x < (*image).width() - 2 && y > 1 && y < (*image).height() - 1)		Grow_region(seed_array, image, output,  tolerancy, x - 1, y);

	
	
}
*/