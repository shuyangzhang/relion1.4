/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
#include "src/projector.h"
//#define DEBUG

void Projector::initialiseData(int current_size)
{
	// By default r_max is half ori_size
	if (current_size < 0)
		r_max = ori_size / 2;
	else
		r_max = current_size / 2;

	// Never allow r_max beyond Nyquist...
	r_max = XMIPP_MIN(r_max, ori_size / 2);

	// Set pad_size
	pad_size = 2 * (padding_factor * r_max + 1) + 1;

	// Short side of data array
	switch (ref_dim)
	{
	case 2:
	   data.resize(pad_size, pad_size / 2 + 1);
	   break;
	case 3:
	   data.resize(pad_size, pad_size, pad_size / 2 + 1);
	   break;
	default:
	   REPORT_ERROR("Projector::resizeData%%ERROR: Dimension of the data array should be 2 or 3");
	}

	// Set origin in the y.z-center, but on the left side for x.
	data.setXmippOrigin();
	data.xinit=0;

}

void Projector::initZeros(int current_size)
{
	initialiseData(current_size);
	data.initZeros();
}

long int Projector::getSize()
{
	// Short side of data array
	switch (ref_dim)
	{
		case 2:
		   return pad_size * (pad_size / 2 + 1);
		   break;
		case 3:
		   return pad_size * pad_size * (pad_size / 2 + 1);
		   break;
		default:
		   REPORT_ERROR("Projector::resizeData%%ERROR: Dimension of the data array should be 2 or 3");
	}

}

// Fill data array with oversampled Fourier transform, and calculate its power spectrum
void Projector::computeFourierTransformMap(MultidimArray<DOUBLE> &vol_in, MultidimArray<DOUBLE> &power_spectrum, int current_size, int nr_threads, bool do_gridding)
{

	MultidimArray<DOUBLE> Mpad;
	MultidimArray<Complex > Faux;
    FourierTransformer transformer;
    // DEBUGGING: multi-threaded FFTWs are giving me a headache?
	// For a long while: switch them off!
	//transformer.setThreadsNumber(nr_threads);
    DOUBLE normfft;

	// Size of padded real-space volume
	int padoridim = padding_factor * ori_size;

	// Initialize data array of the oversampled transform
	ref_dim = vol_in.getDim();

	// Make Mpad
	switch (ref_dim)
	{
	case 2:
	   Mpad.initZeros(padoridim, padoridim);
	   normfft = (DOUBLE)(padding_factor * padding_factor);
	   break;
	case 3:
	   Mpad.initZeros(padoridim, padoridim, padoridim);
	   if (data_dim ==3)
		   normfft = (DOUBLE)(padding_factor * padding_factor * padding_factor);
	   else
		   normfft = (DOUBLE)(padding_factor * padding_factor * padding_factor * ori_size);
	   break;
	default:
	   REPORT_ERROR("Projector::computeFourierTransformMap%%ERROR: Dimension of the data array should be 2 or 3");
	}

	// First do a gridding pre-correction on the real-space map:
	// Divide by the inverse Fourier transform of the interpolator in Fourier-space
	// 10feb11: at least in 2D case, this seems to be the wrong thing to do!!!
	// TODO: check what is best for subtomo!
	if (do_gridding)// && data_dim != 3)
		griddingCorrect(vol_in);

	// Pad translated map with zeros
	vol_in.setXmippOrigin();
	Mpad.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_in) // This will also work for 2D
		A3D_ELEM(Mpad, k, i, j) = A3D_ELEM(vol_in, k, i, j);

	// Translate padded map to put origin of FT in the center
	CenterFFT(Mpad, true);

	// Calculate the oversampled Fourier transform
	transformer.FourierTransform(Mpad, Faux, false);

	// Free memory: Mpad no longer needed
	Mpad.clear();

	// Resize data array to the right size and initialise to zero
	initZeros(current_size);

	// Fill data only for those points with distance to origin less than max_r
	// (other points will be zero because of initZeros() call above
	// Also calculate radial power spectrum
	power_spectrum.initZeros(ori_size / 2 + 1);
	MultidimArray<DOUBLE> counter(power_spectrum);
	counter.initZeros();

	int max_r2 = r_max * r_max * padding_factor * padding_factor;
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux) // This will also work for 2D
	{
		int r2 = kp*kp + ip*ip + jp*jp;
		// The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
		if (r2 <= max_r2)
		{
			// Set data array
			A3D_ELEM(data, kp, ip, jp) = DIRECT_A3D_ELEM(Faux, k, i, j) * normfft;

			// Calculate power spectrum
			int ires = ROUND( sqrt((DOUBLE)r2) / padding_factor );
			// Factor two because of two-dimensionality of the complex plane
			DIRECT_A1D_ELEM(power_spectrum, ires) += norm(A3D_ELEM(data, kp, ip, jp)) / 2.;
			DIRECT_A1D_ELEM(counter, ires) += 1.;
		}
	}

	// Calculate radial average of power spectrum
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(power_spectrum)
	{
		if (DIRECT_A1D_ELEM(counter, i) < 1.)
			DIRECT_A1D_ELEM(power_spectrum, i) = 0.;
		else
			DIRECT_A1D_ELEM(power_spectrum, i) /= DIRECT_A1D_ELEM(counter, i);
	}

	transformer.cleanup();

}

void Projector::griddingCorrect(MultidimArray<DOUBLE> &vol_in)
{

// #define DEBUG_GRID_TYPE

	// Correct real-space map by dividing it by the Fourier transform of the interpolator(s)
	vol_in.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_in)
	{
		DOUBLE r = sqrt((DOUBLE)(k*k+i*i+j*j));
		// if r==0: do nothing (i.e. divide by 1)
		if (r > 0.)
		{
			DOUBLE rval = r / (ori_size * padding_factor);
			DOUBLE sinc = sin(PI * rval) / ( PI * rval);
			//DOUBLE ftblob = blob_Fourier_val(rval, blob) / blob_Fourier_val(0., blob);
			// Interpolation (goes with "interpolator") to go from arbitrary to fine grid
			if (interpolator==NEAREST_NEIGHBOUR && r_min_nn == 0)
			{
				// NN interpolation is convolution with a rectangular pulse, which FT is a sinc function
            	A3D_ELEM(vol_in, k, i, j) /= sinc;

#ifdef DEBUG_GRID_TYPE
    std::cerr << "Now using NN interpolator to do gridding correct " << std::endl;
#endif

			}
			else if (interpolator==TRILINEAR || (interpolator==NEAREST_NEIGHBOUR && r_min_nn > 0) )
			{
				// trilinear interpolation is convolution with a triangular pulse, which FT is a sinc^2 function
				A3D_ELEM(vol_in, k, i, j) /= sinc * sinc;

#ifdef DEBUG_GRID_TYPE
    std::cerr << "Now using linear interpolator to do gridding correct " << std::endl;
#endif

			}
			else if (interpolator==CUBIC)
			{

#ifdef DEBUG_GRID_TYPE
    std::cerr << "Now using cubic interpolator to do gridding correct " << std::endl;
#endif

				// tricubic interpolation is convolution with a kernel as follow:
				// W(x) = (a + 2) * |x|^3 - (a - 3) * |x|^2 + 1, for |x| <= 1;
				//      = a * |x|^3 - 5a * |x|^2 + 8a * |x| - 4a, for 1 < |x| < 2; 
				//      = 0 , otherwise

				// which FT is a function as follow:
				// F[W(x)] = - sinc(pi * x) * 3 * (1 - 1/2 * pi^2 * x^2 * sinc((pi * x/2)^2) ) / (pi^2 * x^2) +
				//           - 5a * sinc(pi * x) * (1 - 1/2 * pi^2 * x^2 * sinc((pi * x/2)^2) ) / (pi^2 * x^2) +
				//             3 * sinc(pi^2 * x^2) / (pi^2 * x^2) +
				//             3a * sinc(pi^2 * x^2) / (pi^2 * x^2) +
				//             3a * sinc(pi^2 * x^2) * (1 - 2 * pi^2 * x^2 * sinc(pi^2 * x^2) ) / (pi^2 * x^2) +
				//           - a * sinc(pi * x) * (1 - 9/2 * pi^2 * x^2 * sinc((3/2 * pi * x)^2) ) / (pi^2 * x^2)
				
				// cubic_factor set -0.5 for test
				DOUBLE cubic_factor = -0.5;
				DOUBLE sincd2 = sin(PI * rval * 0.5) / (PI * rval * 0.5);
				DOUBLE sinc3d2 = sin(PI * rval * 1.5) / (PI * rval * 1.5);
				DOUBLE ft_of_kernel = 0.25 * sincd2 * sincd2 * sinc + 3 * sinc * sinc * sinc * sinc - 2.25 * sinc * sinc3d2 * sinc3d2;

#ifdef DEBUG_OTHER_CUBIC			//  this region of code might be wrong, so annotate first 	
				DOUBLE sinc_square = sin(PI * PI * rval * rval) / (PI * PI * rval * rval);
				DOUBLE sinc_1square2 = sin(PI * PI * rval * rval * 1/2 * 1/2) / (PI * PI * rval * rval * 1/2 * 1/2);
				DOUBLE sinc_3square2 = sin(PI * PI * rval * rval * 3/2 * 3/2) / (PI * PI * rval * rval * 3/2 * 3/2);
				DOUBLE ft_of_kernel = - 3 * (1 - 1/2 * PI * PI * rval * rval * sinc_3square2) * sinc / (PI * PI * rval * rval)
									  - 5 * cubic_factor * (1 - 1/2 * PI * PI * rval * rval * sinc_1square2) * sinc / (PI * PI * rval *rval)
									  + 3 * sinc_square / (PI * PI * rval * rval)
									  + 3 * cubic_factor * sinc_square / (PI * PI * rval * rval)
									  + 3 * cubic_factor * sinc_square * (1 - 2 * PI * PI * rval * rval * sinc_square) / (PI * PI * rval * rval)
									  - cubic_factor * sinc * (1 - 9/2 * PI * PI * rval * rval * sinc_3square2) / (PI * PI * rval * rval);
#endif

				A3D_ELEM(vol_in, k, i, j) /= ft_of_kernel;
			}
			else
				REPORT_ERROR("BUG Projector::griddingCorrect: unrecognised interpolator scheme.");
//#define DEBUG_GRIDDING_CORRECT
#ifdef DEBUG_GRIDDING_CORRECT
			if (k==0 && i==0 && j > 0)
				std::cerr << " j= " << j << " sinc= " << sinc << std::endl;
#endif
		}
	}
}

void Projector::project(MultidimArray<Complex > &f2d, Matrix2D<DOUBLE> &A, bool inv)
{
	DOUBLE fx, fy, fz, xp, yp, zp;
	int x0, x1, y0, y1, z0, z1, y, y2, r2;
	bool is_neg_x;
	Complex d000, d001, d010, d011, d100, d101, d110, d111;
	Complex dx00, dx01, dx10, dx11, dxy0, dxy1;
	Matrix2D<DOUBLE> Ainv;

	DOUBLE u0, u1, u2, u3, v0, v1, v2, v3, w0, w1, w2, w3;
	int intXX, intYY, intZZ;
	// DOUBLE dblXX, dblYY, dblZZ;
	Complex f000, f001, f002, f003, f010, f011, f012, f013, f020, f021, f022, f023, f030, f031, f032, f033,
			f100, f101, f102, f103, f110, f111, f112, f113, f120, f121, f122, f123, f130, f131, f132, f133,
			f200, f201, f202, f203, f210, f211, f212, f213, f220, f221, f222, f223, f230, f231, f232, f233,
			f300, f301, f302, f303, f310, f311, f312, f313, f320, f321, f322, f323, f330, f331, f332, f333;
	Complex fx00, fx01, fx02, fx03, fx10, fx11, fx12, fx13, fx20, fx21, fx22, fx23, fx30, fx31, fx32, fx33,
			fxy0, fxy1, fxy2, fxy3;
	// cubic_factor set -0.5 for test
	DOUBLE cubic_factor = -0.5;
	

    // f2d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside r_max should already be zero...
    // f2d.initZeros();

	// Use the inverse matrix
    if (inv)
    	Ainv = A;
    else
    	Ainv = A.transpose();

    // The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
    int my_r_max = XMIPP_MIN(r_max, XSIZE(f2d) - 1);

    // Go from the 2D slice coordinates to the 3D coordinates
    Ainv *= (DOUBLE)padding_factor;  // take scaling into account directly
    int max_r2 = my_r_max * my_r_max;
    int min_r2_nn = r_min_nn * r_min_nn;

//define DEBUG
#ifdef DEBUG
    std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
    std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
    std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
    std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
	std::cerr << " ZSIZE(data)= "<< ZSIZE(data) << std::endl;
    std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
    std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
    std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
#endif

	for (int i=0; i < YSIZE(f2d); i++)
	{
		// Dont search beyond square with side max_r
		if (i <= my_r_max)
		{
			y = i;
		}
		else if (i >= YSIZE(f2d) - my_r_max)
		{
			y = i - YSIZE(f2d);
		}
		else
			continue;

		y2 = y * y;
		for (int x=0; x <= my_r_max; x++)
		{
	    	// Only include points with radius < max_r (exclude points outside circle in square)
			r2 = x * x + y2;
			if (r2 > max_r2)
				continue;

			// Get logical coordinates in the 3D map
			xp = Ainv(0,0) * x + Ainv(0,1) * y;
			yp = Ainv(1,0) * x + Ainv(1,1) * y;
			zp = Ainv(2,0) * x + Ainv(2,1) * y;

			if (interpolator == TRILINEAR || r2 < min_r2_nn)
			{

				// Only asymmetric half is stored
				if (xp < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
					zp = -zp;
					is_neg_x = true;
				}
				else
				{
					is_neg_x = false;
				}

				// Trilinear interpolation (with physical coords)
				// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
				// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
    			x0 = FLOOR(xp);
				fx = xp - x0;
				x1 = x0 + 1;

				y0 = FLOOR(yp);
				fy = yp - y0;
				y0 -=  STARTINGY(data);
				y1 = y0 + 1;

				z0 = FLOOR(zp);
				fz = zp - z0;
				z0 -= STARTINGZ(data);
				z1 = z0 + 1;

				// Matrix access can be accelerated through pre-calculation of z0*xydim etc.
				d000 = DIRECT_A3D_ELEM(data, z0, y0, x0);
				d001 = DIRECT_A3D_ELEM(data, z0, y0, x1);
				d010 = DIRECT_A3D_ELEM(data, z0, y1, x0);
				d011 = DIRECT_A3D_ELEM(data, z0, y1, x1);
				d100 = DIRECT_A3D_ELEM(data, z1, y0, x0);
				d101 = DIRECT_A3D_ELEM(data, z1, y0, x1);
				d110 = DIRECT_A3D_ELEM(data, z1, y1, x0);
				d111 = DIRECT_A3D_ELEM(data, z1, y1, x1);

				// Set the interpolated value in the 2D output array
				dx00 = LIN_INTERP(fx, d000, d001);
				dx01 = LIN_INTERP(fx, d100, d101);
				dx10 = LIN_INTERP(fx, d010, d011);
				dx11 = LIN_INTERP(fx, d110, d111);
				dxy0 = LIN_INTERP(fy, dx00, dx10);
				dxy1 = LIN_INTERP(fy, dx01, dx11);
				DIRECT_A2D_ELEM(f2d, i, x) = LIN_INTERP(fz, dxy0, dxy1);

				// Take complex conjugated for half with negative x
				if (is_neg_x)
					DIRECT_A2D_ELEM(f2d, i, x) = conj(DIRECT_A2D_ELEM(f2d, i, x));

			} // endif TRILINEAR
			else if (interpolator == NEAREST_NEIGHBOUR )
			{
				x0 = ROUND(xp);
				y0 = ROUND(yp);
				z0 = ROUND(zp);
				if (x0 < 0)
					DIRECT_A2D_ELEM(f2d, i, x) = conj(A3D_ELEM(data, -z0, -y0, -x0));
				else
					DIRECT_A2D_ELEM(f2d, i, x) = A3D_ELEM(data, z0, y0, x0);

			} // endif NEAREST_NEIGHBOUR
			else if (interpolator == CUBIC )
			{
				if (xp < 0)
				{
					xp = -xp;
					yp = -yp;
					zp = -zp;
					is_neg_x = true;
				}
				else
				{
					is_neg_x = false;
				}

				intXX = FLOOR(xp);
				intYY = FLOOR(yp);
				intZZ = FLOOR(zp);

				u1 = xp - intXX;
				v1 = yp - intYY;
				w1 = zp - intZZ;

				// doing this because the index of y and z axis is not correct!
				intYY -=  STARTINGY(data);
				intZZ -=  STARTINGZ(data);

				u0 = 1 + u1;
				u2 = 1 - u1;
				u3 = 2 - u1;

				v0 = 1 + v1;
				v2 = 1 - v1;
				v3 = 2 - v1;

				w0 = 1 + w1;
				w2 = 1 - w1;
				w3 = 2 - w1;

				u0 = - 4 * cubic_factor + 8 * cubic_factor * u0 - 5 * cubic_factor * u0 * u0 + cubic_factor * u0 * u0 * u0;
				u1 = 1 - (cubic_factor + 3) * u1 * u1 + (cubic_factor + 2) * u1 * u1 * u1;
				u2 = 1 - (cubic_factor + 3) * u2 * u2 + (cubic_factor + 2) * u2 * u2 * u2;
				u3 = - 4 * cubic_factor + 8 * cubic_factor * u3 - 5 * cubic_factor * u3 * u3 + cubic_factor * u3 * u3 * u3;
				
				v0 = - 4 * cubic_factor + 8 * cubic_factor * v0 - 5 * cubic_factor * v0 * v0 + cubic_factor * v0 * v0 * v0;
				v1 = 1 - (cubic_factor + 3) * v1 * v1 + (cubic_factor + 2) * v1 * v1 * v1;
				v2 = 1 - (cubic_factor + 3) * v2 * v2 + (cubic_factor + 2) * v2 * v2 * v2;
				v3 = - 4 * cubic_factor + 8 * cubic_factor * v3 - 5 * cubic_factor * v3 * v3 + cubic_factor * v3 * v3 * v3;

				w0 = - 4 * cubic_factor + 8 * cubic_factor * w0 - 5 * cubic_factor * w0 * w0 + cubic_factor * w0 * w0 * w0;
				w1 = 1 - (cubic_factor + 3) * w1 * w1 + (cubic_factor + 2) * w1 * w1 * w1;
				w2 = 1 - (cubic_factor + 3) * w2 * w2 + (cubic_factor + 2) * w2 * w2 * w2;
				w3 = - 4 * cubic_factor + 8 * cubic_factor * w3 - 5 * cubic_factor * w3 * w3 + cubic_factor * w3 * w3 * w3;

#ifdef DEBUG
	std::cerr << " entering a new turn "<< std::endl;
    std::cerr << " i= "<< i << std::endl;
    std::cerr << " x= "<< x << std::endl;
    std::cerr << " xp= "<< xp << std::endl;
    std::cerr << " yp= "<< yp << std::endl;
    std::cerr << " zp= "<< zp << std::endl;
    std::cerr << " intXX= "<< intXX << std::endl;
    std::cerr << " intYY= "<< intYY << std::endl;
    std::cerr << " intZZ= "<< intZZ << std::endl;
//    std::cerr << " Ainv= " << Ainv << std::endl;
#endif


				f000 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX - 1);
				f001 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX );
				f002 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX + 1);
				f003 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX + 2);
				f010 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX - 1);
				f011 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX );
				f012 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX + 1);
				f013 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX + 2);
				f020 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX - 1);
				f021 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX );
				f022 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX + 1);
				f023 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX + 2);
				f030 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX - 1);
				f031 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX );
				f032 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX + 1);
				f033 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX + 2);
				f100 = DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX - 1);
				f101 = DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX );
				f102 = DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX + 1);
				f103 = DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX + 2);
				f110 = DIRECT_A3D_ELEM(data, intZZ , intYY , intXX - 1);
				f111 = DIRECT_A3D_ELEM(data, intZZ , intYY , intXX );
				f112 = DIRECT_A3D_ELEM(data, intZZ , intYY , intXX + 1);
				f113 = DIRECT_A3D_ELEM(data, intZZ , intYY , intXX + 2);
				f120 = DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX - 1);
				f121 = DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX );
				f122 = DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX + 1);
				f123 = DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX + 2);
				f130 = DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX - 1);
				f131 = DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX );
				f132 = DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX + 1);
				f133 = DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX + 2);
				f200 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX - 1);
				f201 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX );
				f202 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX + 1);
				f203 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX + 2);
				f210 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX - 1);
				f211 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX );
				f212 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX + 1);
				f213 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX + 2);
				f220 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX - 1);
				f221 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX );
				f222 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX + 1);
				f223 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX + 2);
				f230 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX - 1);
				f231 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX );
				f232 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX + 1);
				f233 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX + 2);
				f300 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX - 1);
				f301 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX );
				f302 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX + 1);
				f303 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX + 2);
				f310 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX - 1);
				f311 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX );
				f312 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX + 1);
				f313 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX + 2);
				f320 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX - 1);
				f321 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX );
				f322 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX + 1);
				f323 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX + 2);
				f330 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX - 1);
				f331 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX );
				f332 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX + 1);
				f333 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX + 2);

#ifdef DEBUG
	std::cerr << " 64 points for tricubic interpolation is done! "<< std::endl;
    std::cerr << " f000= "<< f000.real << " + " << f000.imag << " i " << std::endl;
	std::cerr << " f111= "<< f000.real << " + " << f000.imag << " i " << std::endl;
	std::cerr << " f222= "<< f000.real << " + " << f000.imag << " i " << std::endl;
	std::cerr << " f333= "<< f000.real << " + " << f000.imag << " i " << std::endl;
#endif

				fx00 = u0 * f000 + u1 * f001 + u2 * f002 + u3 * f003;
				fx01 = u0 * f010 + u1 * f011 + u2 * f012 + u3 * f013;
				fx02 = u0 * f020 + u1 * f021 + u2 * f022 + u3 * f023;
				fx03 = u0 * f030 + u1 * f031 + u2 * f032 + u3 * f033;
				fx10 = u0 * f100 + u1 * f101 + u2 * f102 + u3 * f103;
				fx11 = u0 * f110 + u1 * f111 + u2 * f112 + u3 * f113;
				fx12 = u0 * f120 + u1 * f121 + u2 * f122 + u3 * f123;
				fx13 = u0 * f130 + u1 * f131 + u2 * f132 + u3 * f133;
				fx20 = u0 * f200 + u1 * f201 + u2 * f202 + u3 * f203;
				fx21 = u0 * f210 + u1 * f211 + u2 * f212 + u3 * f213;
				fx22 = u0 * f220 + u1 * f221 + u2 * f222 + u3 * f223;
				fx23 = u0 * f230 + u1 * f231 + u2 * f232 + u3 * f233;
				fx30 = u0 * f300 + u1 * f301 + u2 * f302 + u3 * f303;
				fx31 = u0 * f310 + u1 * f311 + u2 * f312 + u3 * f313;
				fx32 = u0 * f320 + u1 * f321 + u2 * f322 + u3 * f323;
				fx33 = u0 * f330 + u1 * f331 + u2 * f332 + u3 * f333;

#ifdef DEBUG
	std::cerr << " x axis cubic interpolation is done! "<< std::endl;
    std::cerr << " fx00= "<< f000.real << " + " << f000.imag << " i " << std::endl;
	std::cerr << " fx01= "<< f000.real << " + " << f000.imag << " i " << std::endl;
	std::cerr << " fx02= "<< f000.real << " + " << f000.imag << " i " << std::endl;
	std::cerr << " fx03= "<< f000.real << " + " << f000.imag << " i " << std::endl;
	std::cerr << " u0= "<< u0 << std::endl;
    std::cerr << " u1= "<< u1 << std::endl;
    std::cerr << " u2= "<< u2 << std::endl;
    std::cerr << " u3= "<< u3 << std::endl;
#endif

				fxy0 = v0 * fx00 + v1 * fx01 + v2 * fx02 + v3 * fx03;
				fxy1 = v0 * fx10 + v1 * fx11 + v2 * fx12 + v3 * fx13;
				fxy2 = v0 * fx20 + v1 * fx21 + v2 * fx22 + v3 * fx23;
				fxy3 = v0 * fx30 + v1 * fx31 + v2 * fx32 + v3 * fx33;

#ifdef DEBUG
	std::cerr << " y axis cubic interpolation is done! "<< std::endl;
    std::cerr << " fxy0= "<< f000.real << " + " << f000.imag << " i " << std::endl;
	std::cerr << " fxy1= "<< f000.real << " + " << f000.imag << " i " << std::endl;
	std::cerr << " fxy2= "<< f000.real << " + " << f000.imag << " i " << std::endl;
	std::cerr << " fxy3= "<< f000.real << " + " << f000.imag << " i " << std::endl;
	std::cerr << " v0= "<< v0 << std::endl;
    std::cerr << " v1= "<< v1 << std::endl;
    std::cerr << " v2= "<< v2 << std::endl;
    std::cerr << " v3= "<< v3 << std::endl;
#endif

				DIRECT_A2D_ELEM(f2d, i, x) = w0 * fxy0 + w1 * fxy1 + w2 * fxy2 + w3 * fxy3;

				if (is_neg_x)
					DIRECT_A2D_ELEM(f2d, i, x) = conj(DIRECT_A2D_ELEM(f2d, i, x));

			} // endif CUBIC
			else
				REPORT_ERROR("Unrecognized interpolator in Projector::project");

		} // endif x-loop
	} // endif y-loop


#ifdef DEBUG
    std::cerr << "done with project..." << std::endl;
#endif
}

void Projector::rotate2D(MultidimArray<Complex > &f2d, Matrix2D<DOUBLE> &A, bool inv)
{
	DOUBLE fx, fy, xp, yp;
	int x0, x1, y0, y1, y, y2, r2;
	bool is_neg_x;
	Complex d00, d01, d10, d11, dx0, dx1;
	Matrix2D<DOUBLE> Ainv;

	DOUBLE u0, u1, u2, u3, v0, v1, v2, v3;
	int intXX, intYY, intZZ;
	Complex f00, f01, f02, f03, f10, f11, f12, f13, f20, f21, f22, f23, f30, f31, f32, f33;
	Complex fx0, fx1, fx2, fx3;

	// cubic_factor set -0.5 for test
	DOUBLE cubic_factor = -0.5;

    // f2d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero...
    // f2d.initZeros();
	// Use the inverse matrix
    if (inv)
    	Ainv = A;
    else
    	Ainv = A.transpose();

    // The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
    int my_r_max = XMIPP_MIN(r_max, XSIZE(f2d) - 1);

    // Go from the 2D slice coordinates to the map coordinates
    Ainv *= (DOUBLE)padding_factor;  // take scaling into account directly
    int max_r2 = my_r_max * my_r_max;
    int min_r2_nn = r_min_nn * r_min_nn;
#ifdef DEBUG
    std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
    std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
    std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
    std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
    std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
    std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
    std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
#endif
	for (int i=0; i < YSIZE(f2d); i++)
	{
		// Don't search beyond square with side max_r
		if (i <= my_r_max)
		{
			y = i;
		}
		else if (i >= YSIZE(f2d) - my_r_max)
		{
			y = i - YSIZE(f2d);
		}
		else
			continue;
		y2 = y * y;
		for (int x=0; x <= my_r_max; x++)
		{
	    	// Only include points with radius < max_r (exclude points outside circle in square)
			r2 = x * x + y2;
			if (r2 > max_r2)
				continue;

			// Get logical coordinates in the 3D map
			xp = Ainv(0,0) * x + Ainv(0,1) * y;
			yp = Ainv(1,0) * x + Ainv(1,1) * y;
			if (interpolator == TRILINEAR || r2 < min_r2_nn)
			{
				// Only asymmetric half is stored
				if (xp < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
					is_neg_x = true;
				}
				else
				{
					is_neg_x = false;
				}

				// Trilinear interpolation (with physical coords)
				// Subtract STARTINGY to accelerate access to data (STARTINGX=0)
				// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
    			x0 = FLOOR(xp);
				fx = xp - x0;
				x1 = x0 + 1;

				y0 = FLOOR(yp);
				fy = yp - y0;
				y0 -=  STARTINGY(data);
				y1 = y0 + 1;

				// Matrix access can be accelerated through pre-calculation of z0*xydim etc.
				d00 = DIRECT_A2D_ELEM(data, y0, x0);
				d01 = DIRECT_A2D_ELEM(data, y0, x1);
				d10 = DIRECT_A2D_ELEM(data, y1, x0);
				d11 = DIRECT_A2D_ELEM(data, y1, x1);

				// Set the interpolated value in the 2D output array
				dx0 = LIN_INTERP(fx, d00, d01);
				dx1 = LIN_INTERP(fx, d10, d11);
				DIRECT_A2D_ELEM(f2d, i, x) = LIN_INTERP(fy, dx0, dx1);
				// Take complex conjugated for half with negative x
				if (is_neg_x)
					DIRECT_A2D_ELEM(f2d, i, x) = conj(DIRECT_A2D_ELEM(f2d, i, x));
			} // endif TRILINEAR
			else if (interpolator == NEAREST_NEIGHBOUR )
			{
				x0 = ROUND(xp);
				y0 = ROUND(yp);
				if (x0 < 0)
					DIRECT_A2D_ELEM(f2d, i, x) = conj(A2D_ELEM(data, -y0, -x0));
				else
					DIRECT_A2D_ELEM(f2d, i, x) = A2D_ELEM(data, y0, x0);
			} // endif NEAREST_NEIGHBOUR
			else if (interpolator == CUBIC )
			{
				if (xp < 0)
				{
					xp = -xp;
					yp = -yp;
					is_neg_x = true;
				}
				else 
				{
					is_neg_x = false;
				}

				intXX = FLOOR(xp);
				intYY = FLOOR(yp);

				u1 = xp - intXX;
				v1 = yp - intYY;

				intYY -=  STARTINGY(data);

				u0 = 1 + u1;
				u2 = 1 - u1;
				u3 = 2 - u1;

				v0 = 1 + v1;
				v2 = 1 - v1;
				v3 = 2 - v1;

				u0 = - 4 * cubic_factor + 8 * cubic_factor * u0 - 5 * cubic_factor * u0 * u0 + cubic_factor * u0 * u0 * u0;
				u1 = 1 - (cubic_factor + 3) * u1 * u1 + (cubic_factor + 2) * u1 * u1 * u1;
				u2 = 1 - (cubic_factor + 3) * u2 * u2 + (cubic_factor + 2) * u2 * u2 * u2;
				u3 = - 4 * cubic_factor + 8 * cubic_factor * u3 - 5 * cubic_factor * u3 * u3 + cubic_factor * u3 * u3 * u3;
				
				v0 = - 4 * cubic_factor + 8 * cubic_factor * v0 - 5 * cubic_factor * v0 * v0 + cubic_factor * v0 * v0 * v0;
				v1 = 1 - (cubic_factor + 3) * v1 * v1 + (cubic_factor + 2) * v1 * v1 * v1;
				v2 = 1 - (cubic_factor + 3) * v2 * v2 + (cubic_factor + 2) * v2 * v2 * v2;
				v3 = - 4 * cubic_factor + 8 * cubic_factor * v3 - 5 * cubic_factor * v3 * v3 + cubic_factor * v3 * v3 * v3;

				f00 = DIRECT_A2D_ELEM(data, intYY - 1, intXX - 1);
				f01 = DIRECT_A2D_ELEM(data, intYY - 1, intXX );
				f02 = DIRECT_A2D_ELEM(data, intYY - 1, intXX + 1);
				f03 = DIRECT_A2D_ELEM(data, intYY - 1, intXX + 2);
				f10 = DIRECT_A2D_ELEM(data, intYY , intXX - 1);
				f11 = DIRECT_A2D_ELEM(data, intYY , intXX );
				f12 = DIRECT_A2D_ELEM(data, intYY , intXX + 1);
				f13 = DIRECT_A2D_ELEM(data, intYY , intXX + 2);
				f20 = DIRECT_A2D_ELEM(data, intYY + 1, intXX - 1);
				f21 = DIRECT_A2D_ELEM(data, intYY + 1, intXX );
				f22 = DIRECT_A2D_ELEM(data, intYY + 1, intXX + 1);
				f23 = DIRECT_A2D_ELEM(data, intYY + 1, intXX + 2);
				f30 = DIRECT_A2D_ELEM(data, intYY + 2, intXX - 1);
				f31 = DIRECT_A2D_ELEM(data, intYY + 2, intXX );
				f32 = DIRECT_A2D_ELEM(data, intYY + 2, intXX + 1);
				f33 = DIRECT_A2D_ELEM(data, intYY + 2, intXX + 2);

				fx0 = u0 * f00 + u1 * f01 + u2 * f02 + u3 * f03;
				fx1 = u0 * f10 + u1 * f11 + u2 * f12 + u3 * f13;
				fx2 = u0 * f20 + u1 * f21 + u2 * f22 + u3 * f23;
				fx3 = u0 * f30 + u1 * f31 + u2 * f32 + u3 * f33;

				DIRECT_A2D_ELEM(f2d, i, x) = v0 * fx0 + v1 * fx1 + v2 * fx2 + v3 * fx3;

				if (is_neg_x)
					DIRECT_A2D_ELEM(f2d, i, x) = conj(DIRECT_A2D_ELEM(f2d, i, x));
			} // endif CUBIC
			else
				REPORT_ERROR("Unrecognized interpolator in Projector::project");
		} // endif x-loop
	} // endif y-loop
}


void Projector::rotate3D(MultidimArray<Complex > &f3d, Matrix2D<DOUBLE> &A, bool inv)
{
	DOUBLE fx, fy, fz, xp, yp, zp;
	int x0, x1, y0, y1, z0, z1, y, z, y2, z2, r2;
	bool is_neg_x;
	Complex d000, d010, d100, d110, d001, d011, d101, d111, dx00, dx10, dxy0, dx01, dx11, dxy1;
	Matrix2D<DOUBLE> Ainv;

	DOUBLE u0, u1, u2, u3, v0, v1, v2, v3, w0, w1, w2, w3;
	int intXX, intYY, intZZ;

	Complex f000, f001, f002, f003, f010, f011, f012, f013, f020, f021, f022, f023, f030, f031, f032, f033,
			f100, f101, f102, f103, f110, f111, f112, f113, f120, f121, f122, f123, f130, f131, f132, f133,
			f200, f201, f202, f203, f210, f211, f212, f213, f220, f221, f222, f223, f230, f231, f232, f233,
			f300, f301, f302, f303, f310, f311, f312, f313, f320, f321, f322, f323, f330, f331, f332, f333;
	Complex fx00, fx01, fx02, fx03, fx10, fx11, fx12, fx13, fx20, fx21, fx22, fx23, fx30, fx31, fx32, fx33,
			fxy0, fxy1, fxy2, fxy3;

	// cubic_factor set -0.5 for test
	DOUBLE cubic_factor = -0.5;

    // f3d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero...
    // f3d.initZeros();
	// Use the inverse matrix
    if (inv)
    	Ainv = A;
    else
    	Ainv = A.transpose();

    // The f3d image may be smaller than r_max, in that case also make sure not to fill the corners!
    int my_r_max = XMIPP_MIN(r_max, XSIZE(f3d) - 1);

    // Go from the 3D rotated coordinates to the original map coordinates
    Ainv *= (DOUBLE)padding_factor;  // take scaling into account directly
    int max_r2 = my_r_max * my_r_max;
    int min_r2_nn = r_min_nn * r_min_nn;
#ifdef DEBUG
    std::cerr << " XSIZE(f3d)= "<< XSIZE(f3d) << std::endl;
    std::cerr << " YSIZE(f3d)= "<< YSIZE(f3d) << std::endl;
    std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
    std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
    std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
    std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
    std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
#endif

// #define DEBUG_INTER_TYPE

	for (int k=0; k < ZSIZE(f3d); k++)
	{
		// Don't search beyond square with side max_r
		if (k <= my_r_max)
		{
			z = k;
		}
		else if (k >= ZSIZE(f3d) - my_r_max)
		{
			z = k - ZSIZE(f3d);
		}
		else
			continue;
		z2 = z * z;

		for (int i=0; i < YSIZE(f3d); i++)
		{
			// Don't search beyond square with side max_r
			if (i <= my_r_max)
			{
				y = i;
			}
			else if (i >= YSIZE(f3d) - my_r_max)
			{
				y = i - YSIZE(f3d);
			}
			else
				continue;
			y2 = y * y;

			for (int x=0; x <= my_r_max; x++)
			{
				// Only include points with radius < max_r (exclude points outside circle in square)
				r2 = x * x + y2 + z2;
				if (r2 > max_r2)
					continue;

				// Get logical coordinates in the 3D map
				xp = Ainv(0,0) * x + Ainv(0,1) * y + Ainv(0,2) * z;
				yp = Ainv(1,0) * x + Ainv(1,1) * y + Ainv(1,2) * z;
				zp = Ainv(2,0) * x + Ainv(2,1) * y + Ainv(2,2) * z;

				if (interpolator == TRILINEAR || r2 < min_r2_nn)
				{

#ifdef DEBUG_INTER_TYPE
    std::cerr << "Now using linear interpolator to rotate3D " << std::endl;
#endif

					// Only asymmetric half is stored
					if (xp < 0)
					{
						// Get complex conjugated hermitian symmetry pair
						xp = -xp;
						yp = -yp;
						zp = -zp;
						is_neg_x = true;
					}
					else
					{
						is_neg_x = false;
					}

					// Trilinear interpolation (with physical coords)
					// Subtract STARTINGY to accelerate access to data (STARTINGX=0)
					// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
					x0 = FLOOR(xp);
					fx = xp - x0;
					x1 = x0 + 1;

					y0 = FLOOR(yp);
					fy = yp - y0;
					y0 -=  STARTINGY(data);
					y1 = y0 + 1;

					z0 = FLOOR(zp);
					fz = zp - z0;
					z0 -=  STARTINGZ(data);
					z1 = z0 + 1;

					// Matrix access can be accelerated through pre-calculation of z0*xydim etc.
					d000 = DIRECT_A3D_ELEM(data, z0, y0, x0);
					d001 = DIRECT_A3D_ELEM(data, z0, y0, x1);
					d010 = DIRECT_A3D_ELEM(data, z0, y1, x0);
					d011 = DIRECT_A3D_ELEM(data, z0, y1, x1);
					d100 = DIRECT_A3D_ELEM(data, z1, y0, x0);
					d101 = DIRECT_A3D_ELEM(data, z1, y0, x1);
					d110 = DIRECT_A3D_ELEM(data, z1, y1, x0);
					d111 = DIRECT_A3D_ELEM(data, z1, y1, x1);

					// Set the interpolated value in the 2D output array
					// interpolate in x
					dx00 = LIN_INTERP(fx, d000, d001);
					dx01 = LIN_INTERP(fx, d100, d101);
					dx10 = LIN_INTERP(fx, d010, d011);
					dx11 = LIN_INTERP(fx, d110, d111);
					// interpolate in y
					dxy0 = LIN_INTERP(fy, dx00, dx10);
					dxy1 = LIN_INTERP(fy, dx01, dx11);
					//interpolate in z
					DIRECT_A3D_ELEM(f3d, k, i, x) = LIN_INTERP(fz, dxy0, dxy1);

					// Take complex conjugated for half with negative x
					if (is_neg_x)
						DIRECT_A3D_ELEM(f3d, k, i, x) = conj(DIRECT_A3D_ELEM(f3d, k, i, x));

				} // endif TRILINEAR
				else if (interpolator == NEAREST_NEIGHBOUR )
				{

#ifdef DEBUG_INTER_TYPE
    std::cerr << "Now using NN interpolator to rotate3D " << std::endl;
#endif

					x0 = ROUND(xp);
					y0 = ROUND(yp);
					z0 = ROUND(zp);

					if (x0 < 0)
						DIRECT_A3D_ELEM(f3d, k, i, x) = conj(A3D_ELEM(data, -z0, -y0, -x0));
					else
						DIRECT_A3D_ELEM(f3d, k, i, x) = A3D_ELEM(data, z0, y0, x0);

				} // endif NEAREST_NEIGHBOUR
				else if (interpolator == CUBIC )
				{

#ifdef DEBUG_INTER_TYPE
    std::cerr << "Now using cubic interpolator to rotate3D " << std::endl;
#endif

					if (xp < 0)
					{
						xp = -xp;
						yp = -yp;
						zp = -zp;
						is_neg_x = true;
					}
					else
					{
						is_neg_x = false;
					}

					intXX = FLOOR(xp);
					intYY = FLOOR(yp);
					intZZ = FLOOR(zp);

					u1 = xp - intXX;
					v1 = yp - intYY;
					w1 = zp - intZZ;

					intYY -=  STARTINGY(data);
					intZZ -=  STARTINGZ(data);

					u0 = 1 + u1;
					u2 = 1 - u1;
					u3 = 2 - u1;

					v0 = 1 + v1;
					v2 = 1 - v1;
					v3 = 2 - v1;

					w0 = 1 + w1;
					w2 = 1 - w1;
					w3 = 2 - w1;

					u0 = - 4 * cubic_factor + 8 * cubic_factor * u0 - 5 * cubic_factor * u0 * u0 + cubic_factor * u0 * u0 * u0;
					u1 = 1 - (cubic_factor + 3) * u1 * u1 + (cubic_factor + 2) * u1 * u1 * u1;
					u2 = 1 - (cubic_factor + 3) * u2 * u2 + (cubic_factor + 2) * u2 * u2 * u2;
					u3 = - 4 * cubic_factor + 8 * cubic_factor * u3 - 5 * cubic_factor * u3 * u3 + cubic_factor * u3 * u3 * u3;
					
					v0 = - 4 * cubic_factor + 8 * cubic_factor * v0 - 5 * cubic_factor * v0 * v0 + cubic_factor * v0 * v0 * v0;
					v1 = 1 - (cubic_factor + 3) * v1 * v1 + (cubic_factor + 2) * v1 * v1 * v1;
					v2 = 1 - (cubic_factor + 3) * v2 * v2 + (cubic_factor + 2) * v2 * v2 * v2;
					v3 = - 4 * cubic_factor + 8 * cubic_factor * v3 - 5 * cubic_factor * v3 * v3 + cubic_factor * v3 * v3 * v3;

					w0 = - 4 * cubic_factor + 8 * cubic_factor * w0 - 5 * cubic_factor * w0 * w0 + cubic_factor * w0 * w0 * w0;
					w1 = 1 - (cubic_factor + 3) * w1 * w1 + (cubic_factor + 2) * w1 * w1 * w1;
					w2 = 1 - (cubic_factor + 3) * w2 * w2 + (cubic_factor + 2) * w2 * w2 * w2;
					w3 = - 4 * cubic_factor + 8 * cubic_factor * w3 - 5 * cubic_factor * w3 * w3 + cubic_factor * w3 * w3 * w3;

					f000 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX - 1);
					f001 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX );
					f002 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX + 1);
					f003 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX + 2);
					f010 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX - 1);
					f011 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX );
					f012 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX + 1);
					f013 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX + 2);
					f020 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX - 1);
					f021 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX );
					f022 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX + 1);
					f023 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX + 2);
					f030 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX - 1);
					f031 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX );
					f032 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX + 1);
					f033 = DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX + 2);
					f100 = DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX - 1);
					f101 = DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX );
					f102 = DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX + 1);
					f103 = DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX + 2);
					f110 = DIRECT_A3D_ELEM(data, intZZ , intYY , intXX - 1);
					f111 = DIRECT_A3D_ELEM(data, intZZ , intYY , intXX );
					f112 = DIRECT_A3D_ELEM(data, intZZ , intYY , intXX + 1);
					f113 = DIRECT_A3D_ELEM(data, intZZ , intYY , intXX + 2);
					f120 = DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX - 1);
					f121 = DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX );
					f122 = DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX + 1);
					f123 = DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX + 2);
					f130 = DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX - 1);
					f131 = DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX );
					f132 = DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX + 1);
					f133 = DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX + 2);
					f200 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX - 1);
					f201 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX );
					f202 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX + 1);
					f203 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX + 2);
					f210 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX - 1);
					f211 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX );
					f212 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX + 1);
					f213 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX + 2);
					f220 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX - 1);
					f221 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX );
					f222 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX + 1);
					f223 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX + 2);
					f230 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX - 1);
					f231 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX );
					f232 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX + 1);
					f233 = DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX + 2);
					f300 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX - 1);
					f301 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX );
					f302 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX + 1);
					f303 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX + 2);
					f310 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX - 1);
					f311 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX );
					f312 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX + 1);
					f313 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX + 2);
					f320 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX - 1);
					f321 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX );
					f322 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX + 1);
					f323 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX + 2);
					f330 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX - 1);
					f331 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX );
					f332 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX + 1);
					f333 = DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX + 2);
					
					fx00 = u0 * f000 + u1 * f001 + u2 * f002 + u3 * f003;
					fx01 = u0 * f010 + u1 * f011 + u2 * f012 + u3 * f013;
					fx02 = u0 * f020 + u1 * f021 + u2 * f022 + u3 * f023;
					fx03 = u0 * f030 + u1 * f031 + u2 * f032 + u3 * f033;
					fx10 = u0 * f100 + u1 * f101 + u2 * f102 + u3 * f103;
					fx11 = u0 * f110 + u1 * f111 + u2 * f112 + u3 * f113;
					fx12 = u0 * f120 + u1 * f121 + u2 * f122 + u3 * f123;
					fx13 = u0 * f130 + u1 * f131 + u2 * f132 + u3 * f133;
					fx20 = u0 * f200 + u1 * f201 + u2 * f202 + u3 * f203;
					fx21 = u0 * f210 + u1 * f211 + u2 * f212 + u3 * f213;
					fx22 = u0 * f220 + u1 * f221 + u2 * f222 + u3 * f223;
					fx23 = u0 * f230 + u1 * f231 + u2 * f232 + u3 * f233;
					fx30 = u0 * f300 + u1 * f301 + u2 * f302 + u3 * f303;
					fx31 = u0 * f310 + u1 * f311 + u2 * f312 + u3 * f313;
					fx32 = u0 * f320 + u1 * f321 + u2 * f322 + u3 * f323;
					fx33 = u0 * f330 + u1 * f331 + u2 * f332 + u3 * f333;

					fxy0 = v0 * fx00 + v1 * fx01 + v2 * fx02 + v3 * fx03;
					fxy1 = v0 * fx10 + v1 * fx11 + v2 * fx12 + v3 * fx13;
					fxy2 = v0 * fx20 + v1 * fx21 + v2 * fx22 + v3 * fx23;
					fxy3 = v0 * fx30 + v1 * fx31 + v2 * fx32 + v3 * fx33;

					DIRECT_A3D_ELEM(f3d, k, i, x) = w0 * fxy0 + w1 * fxy1 + w2 * fxy2 + w3 * fxy3;

					if (is_neg_x)
						DIRECT_A3D_ELEM(f3d, k, i, x) = conj(DIRECT_A3D_ELEM(f3d, k, i, x));
				} // endif CUBIC
				else
					REPORT_ERROR("Unrecognized interpolator in Projector::project");
			} // endif x-loop
		} // endif y-loop
	} // endif z-loop
}







