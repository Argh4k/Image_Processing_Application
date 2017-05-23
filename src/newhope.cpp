
#include <iostream>
#include <string>
#include <cstdlib>
#include <regex>
#include <fstream>
#include "CImg.h"
#include "IMAGE.h"
#include "Fourier.h"


using namespace cimg_library;
using namespace std;




int main(int argc, char * argv[]) {

	
	//CImg<float> lena = Load_Image("lena.bmp");

	//complex<double>  ** After = FFT(lena);
	//complex<double>  ** After = FFT(lena);
	//cout << endl << "done";
	//savefourier(After, lena.width(), lena.height());

	//CImg<float> mask = Load_Image("F5mask1.bmp");
	//cout << "Lena loadaded" << endl;
	//complex<double>  ** After = PhaseMod(lena,5,1);
	/*
	double newvalue;
	for (int x = 0; x < lena.width(); x++)
	{
		for (int y = 0; y < lena.height(); y++)
		{
			cout << "At coordinates " << x << y << " Real Value - " << After[y][x].real() << " Imaginary - " << After[y][x].imag() << endl;
		}
	}
	*/
	//lena = *PrintMask(IFFT(After, lena.width(), lena.height()), lena.width(), lena.height());
	//lena.save("output.bmp");
	//cin.get();

	//return 0;

	
	CImg<float> lena(400, 400, 1, 3, 0);
	lena.save("output.bmp");
	return 0;

	float(*operations[3])(float, float) = { brightlut, contrastlut, negativelut };

	string operation_to_do;
	string path_to_file;
	float modificator;
	regex pattern_comparisionts("(?:--)(mse|pmse|snr|psnr|md)");
	regex pattern_basic_operations("(?:--)(brightness|contrast|negative|vflip|hflip|dflip|shrink|enlarge|median|gmean|orobertsi)");
	regex pattern_histograms("(?:--)(histogram|slowpass|huniform|cmean|cvariance|cstdev|cvarcoi|casyco|cfasyco|cvarcoii|centropy)");
	regex pattern_morphology("(?:--)(dilation|erosion|opening|closing|HMT|M11|M12|M13|rgrowing)");
	regex pattern_fourier("(?:--)(dft|idft|F1|F2|F3|F4|F5|F6|)");

	if (argc < 2 || argc > 6)
	{
		cout << "Invalid number of arguments, if you need help use --help argument";
		return 0;
	}

	//CHECK THE ARGUMENTS

#pragma region HELP
	if (argc == 2)
		{
			if ((string)argv[1] == "--help")
			{
				ifstream help_file;
				help_file.open("help.txt");
				if (help_file.is_open())
				{

					string line;
					while (!help_file.eof())
					{
						getline(help_file, line);
						cout << line << endl;
					}
					help_file.close();
				}
				else
					cout << "Couldnt find help file";
				return 0;
			}

		}

#pragma endregion


	//comparisions 

	if (regex_search(argv[1], pattern_comparisionts)) //inf regex match go into coparisions part of main 
	{

		CImg<float> image_with_noise = Load_Image(argv[2]);
		CImg<float> image_without_noise = Load_Image(argv[3]);


		if ((string)argv[1] == "--mse")
		{

			Mean_square_error(image_without_noise, image_with_noise);

			return 0;
		}
		else if ((string)argv[1] == "--pmse") {

			Peak_mean_square_error(image_without_noise, image_with_noise);
			return 0;

		}
		else if ((string)argv[1] == "--snr")
		{

			Signal_to_noise_error(image_without_noise, image_with_noise);
			return 0;
		}
		else if ((string)argv[1] == "--psnr")
		{

			Peak_signal_to_noise_error(image_without_noise, image_with_noise);
			return 0;
		}
		else if ((string)argv[1] == "--md")
		{

			Maximum_difference(image_without_noise, image_with_noise);
			return 0;
		}
		else
		{
			cout << "(comparisions)Invalid operation\n If you need help, run program with --help parameter" << endl;
			return 0;
		}
	}

	//basic operations

	else if (regex_search(argv[1], pattern_basic_operations))
	{
		CImg <float> image = Load_Image(argv[2]);
		//operation
		operation_to_do = argv[1];

		// CHOISE OF THE OPERATION

		//OPERATIONS PERFORMED ON THE ORIGINAL IMAGE
		if ((string)operation_to_do == "--brightness")
		{
			modificator = stod(argv[3]);
			basicoperations(modificator, image, operations[0]);
			SaveImage(image);

		}
		else if ((string)operation_to_do == "--contrast")
		{
			modificator = stof(argv[3]);
			basicoperations(modificator, image, operations[1]);
			SaveImage(image);
		}
		else if ((string)operation_to_do == "--negative")
		{
			//modificator = stod(argv[3]);
			basicoperations(0, image, operations[2]);
			SaveImage(image);
		}
		else if ((string)operation_to_do == "--vflip") {

			verticalFlip(image);
			SaveImage(image);
		}
		else if ((string)operation_to_do == "--hflip")
		{
			horizontalFlip(image);
			SaveImage(image);
		}
		else if ((string)operation_to_do == "--dflip") {
			diagonalFlip(image);
			SaveImage(image);
		}

		//OPERATIONS PERFORMED ON THE COPY 
		else if ((string)operation_to_do == "--shrink")
		{
			SaveImage(shrink(image));
		}
		else if ((string)operation_to_do == "--enlarge")
		{
			SaveImage(enlarge(image));
		}
		else if ((string)operation_to_do == "--median")
		{
			CImg<float> * filterimage = medianfilter(image);
			SaveImage(*filterimage);
			delete filterimage;
		}
		else if ((string)operation_to_do == "--gmean")
		{
			CImg<float> * filterimage = geometricfilter(image);
			SaveImage(*filterimage);
			delete filterimage;
		}
		else if ((string)operation_to_do == "--orobertsi")
		{
			CImg<float> *filtredimage = Roberts_operator(image);
			SaveImage(*filtredimage);
			delete filtredimage;
		}
	}

	else if (regex_search(argv[1], pattern_histograms)) {

		CImg<float> image = Load_Image(argv[2]);
		int third_argument = 0;

		if (argc == 4)
			third_argument = stod(argv[3]); //channel or mask

		if ((string)argv[1] == "--histogram")
		{
			if (third_argument < 0 || third_argument > 2)
				cout << "You can create histogram only for 0, 1, 2 channel";
			else
				Createhistogramimage(image, third_argument);
			return 0;
		}
		else if ((string)argv[1] == "--huniform")
		{
			UniformFinalProbabilityDensityFunction(image);
			return 0;
		}
		else if ((string)argv[1] == "--slowpass")
		{
			if (third_argument < 1 || third_argument > 3)
			{
				cout << endl << "Masks can be choosen form 1 to 3";
				return 0;
			}
			else if (third_argument == 1)
			{
				CImg<float> *filtredimage = ameanfilter(image);
				SaveImage(*filtredimage);
				return 0;
			}
			else
			{
				CImg<float> *filtredimage = Low_pass_filter(image, third_argument);
				SaveImage(*filtredimage);
				delete filtredimage;
				return 0;
			}
		}
		else if ((string)argv[1] == "--cmean")
		{
			cout << Cmean(image, third_argument) << endl;
			return 0;
		}
		else if ((string)argv[1] == "--cvariance")
		{
			cout << Cvariance(image, third_argument) << endl;
			return 0;
		}
		else if ((string)argv[1] == "--cstdev")
		{
			cout << Cstdev(image, third_argument) << endl;
			return 0;
		}
		else if ((string)argv[1] == "--cvarcoi")
		{
			cout << Cvarcoi(image, third_argument) << endl;
			return 0;
		}
		else if ((string)argv[1] == "--casyco")
		{
			cout << Casyco(image, third_argument) << endl;
			return 0;
		}
		else if ((string)argv[1] == "--cfasyco")
		{
			cout << Cfasyco(image, third_argument) << endl;
			return 0;
		}
		else if ((string)argv[1] == "--cvarcoii")
		{
			cout << Cvarcoii(image, third_argument) << endl;
			return 0;
		}
		else if ((string)argv[1] == "--centropy")
		{
			cout << Centropy(image, third_argument) << endl;
			return 0;
		}
		return 0;
	}

	else if (regex_search(argv[1], pattern_morphology))
	{
		if (argc < 4)
		{
			cout << "Lack of some arguments \n";
			return 0;
		}

		CImg<float> image;
		image.load(argv[2]);
		image.channel(0);
		int third_arg = 0;
		third_arg = stoi(argv[3]);
		StructuralElement Element(third_arg);

		if ((string)(argv[1]) == "--erosion")
		{
			CImg<float> * aimage = Erosion(image, Element);
			SaveImage(*aimage);
			delete aimage;
			return 0;
		}

		if ((string)(argv[1]) == "--dilation")
		{
			CImg<float> * aimage = Dilation(image, Element);
			SaveImage(*aimage);
			delete aimage;
			return 0;
		}

		if ((string)(argv[1]) == "--opening")
		{
			CImg<float> * aimage = Erosion(image, Element);
			CImg<float> * aimage2 = Dilation(*aimage, Element);
			SaveImage(*aimage2);
			delete aimage;
			delete aimage2;
			return 0;
		}

		if ((string)(argv[1]) == "--closing")
		{
			CImg<float> * aimage = Dilation(image, Element);
			CImg<float> * aimage2 = Erosion(*aimage, Element);
			SaveImage(*aimage2);
			delete aimage;
			delete aimage2;
			return 0;
		}

		if ((string)(argv[1]) == "--HMT")
		{
			CImg<float> * aimage = HMT(image, Element);
			SaveImage(*aimage);
			delete aimage;
			return 0;
		}

		if ((string)(argv[1]) == "--M11")
		{
			CImg<float> * aimg = Dilation(image, Element);
			Difference((*aimg), image);
			SaveImage(*aimg);
			delete aimg;
			return 0;
		}

		if ((string)(argv[1]) == "--M13")
		{
			CImg<float> * aimg = Dilation(image, Element);
			CImg<float> * aimg2 = Erosion(image, Element);
			Difference((*aimg), (*aimg2));
			SaveImage(*aimg);
			delete aimg;
			delete aimg2;
			return 0;
		}

		if ((string)(argv[1]) == "--M12")
		{
			CImg<float> * aimg = Erosion(image, Element);
			Difference(image, (*aimg));
			SaveImage(image);
			delete aimg;
			return 0;
		}
		if ((string)(argv[1]) == "--rgrowing")
		{
			float fourth_arg = atof(argv[3]);
			int x = atoi(argv[4]);
			int y = atoi(argv[5]);
			CImg<float>* output_image = Region_growing(&image, third_arg,x,y);
			SaveImage(*output_image);
			delete output_image;
			return 0;
		}
	}
	else if (regex_search(argv[1], pattern_fourier))
	{
		if ((string)(argv[1]) == "--dft")
		{
			CImg<float> image = Load_Image(argv[2]);
			cout << "picure loadaded" << endl;
			complex<double>  ** After = DFT(image);
			image = *PrintMask(After, image.width(),image.height());
			imageswap(image);
			SaveImage(image);

			for (int i = 0; i < image.height(); ++i) {
				delete[] After[i];
			}
			delete[] After;
		}


		else if ((string)(argv[1]) == "--fft")
		{
			CImg<float> image = Load_Image(argv[2]);
			cout << "picure loadaded" << endl;
			complex<double>  ** After = FFT(image);
			image = *PrintMask(After, image.width(), image.height());
			SaveImage(image);

			for (int i = 0; i < image.height(); ++i) {
				delete[] After[i];
			}
			delete[] After;
		}
		else if ((string)(argv[1]) == "--f1")
		{
			CImg<float> image = Load_Image(argv[2]);
			cout << "picure loadaded" << endl;
			int radius = atoi(argv[3]);
			complex<double>  ** After = LowPassFilter(image, radius);
			image = *Print(IFFT(After, image.width(), image.height()), image.width(), image.height());
			SaveImage(image);

			for (int i = 0; i < image.height(); ++i) {
				delete[] After[i];
			}
			delete[] After;
		}
		else if ((string)(argv[1]) == "--f2")
		{
			CImg<float> image = Load_Image(argv[2]);
			cout << "picure loadaded" << endl;
			int radius = atoi(argv[3]);
			complex<double>  ** After = HighPassFilter(image, radius);
			image = *Print(IFFT(After, image.width(), image.height()), image.width(), image.height());
			SaveImage(image);

			for (int i = 0; i < image.height(); ++i) {
				delete[] After[i];
			}
			delete[] After;
		}
		else if ((string)(argv[1]) == "--f3")
		{
			CImg<float> image = Load_Image(argv[2]);
			cout << "picure loadaded" << endl;
			int arg3 = atoi(argv[3]);
			int arg4 = atoi(argv[4]);
			complex<double>  ** After = BandPassFilter(image, arg3, arg4);
			image = *Print(IFFT(After, image.width(), image.height()), image.width(), image.height());
			SaveImage(image);

			for (int i = 0; i < image.height(); ++i) {
				delete[] After[i];
			}
			delete[] After;
		}
		else if ((string)(argv[1]) == "--f4")
		{
			CImg<float> image = Load_Image(argv[2]);
			cout << "picure loadaded" << endl;
			int arg3 = atoi(argv[3]);
			int arg4 = atoi(argv[4]);
			complex<double>  ** After = BandCutFilter(image, arg3, arg4);
			image = *Print(IFFT(After, image.width(), image.height()), image.width(), image.height());
			SaveImage(image);

			for (int i = 0; i < image.height(); ++i) {
				delete[] After[i];
			}
			delete[] After;
		}
		else if ((string)(argv[1]) == "--f5")
		{
			CImg<float> image = Load_Image(argv[2]);
			cout << "picure loadaded" << endl;
			int arg3 = atoi(argv[3]);
			complex<double>  ** After = MaskFilter( arg3, image);
			image = *Print(IFFT(After, image.width(), image.height()), image.width(), image.height());
			SaveImage(image);

			for (int i = 0; i < image.height(); ++i) {
				delete[] After[i];
			}
			delete[] After;
		}
		else if ((string)(argv[1]) == "--f6")
		{
			CImg<float> image = Load_Image(argv[2]);
			cout << "picure loadaded" << endl;
			int arg3 = atoi(argv[3]);
			int arg4 = atoi(argv[4]);
			complex<double>  ** After = PhaseMod(image, arg3, arg4);
			image = *Print(IFFT(After, image.width(), image.height()), image.width(), image.height());
			SaveImage(image);
			for (int i = 0; i < image.height(); ++i) {
				delete[] After[i];
			}
			delete[] After;
		}



	}
	else { cout << "Error couldnt match to any known function :C " << endl;  return 0; }

	return 0;
}