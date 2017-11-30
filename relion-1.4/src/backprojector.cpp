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
/*
 * backprojector.cpp
 *
 *  Created on: 24 Aug 2010
 *      Author: scheres
 */

#include "src/backprojector.h"

void BackProjector::initialiseDataAndWeight(int current_size)
{

	initialiseData(current_size);
	weight.resize(data);

}

void BackProjector::initZeros(int current_size)
{

	initialiseDataAndWeight(current_size);
	data.initZeros();
	weight.initZeros();
}

void BackProjector::backproject(const MultidimArray<Complex > &f2d,
		                        const Matrix2D<DOUBLE> &A, bool inv,
		                        const MultidimArray<DOUBLE> *Mweight)
{
	DOUBLE fx, fy, fz, mfx, mfy, mfz, xp, yp, zp;
	int first_x, x0, x1, y0, y1, z0, z1, y, y2, r2;
	bool is_neg_x;
	DOUBLE dd000, dd001, dd010, dd011, dd100, dd101, dd110, dd111;
	Complex my_val;
	Matrix2D<DOUBLE> Ainv;
	DOUBLE my_weight = 1.;

	int intXX, intYY, intZZ;
	DOUBLE u0, u1, u2, u3, v0, v1, v2, v3, w0, w1, w2, w3;
	DOUBLE df000, df001, df002, df003, df010, df011, df012, df013, df020, df021, df022, df023, df030, df031, df032, df033,
		   df100, df101, df102, df103, df110, df111, df112, df113, df120, df121, df122, df123, df130, df131, df132, df133,
		   df200, df201, df202, df203, df210, df211, df212, df213, df220, df221, df222, df223, df230, df231, df232, df233,
		   df300, df301, df302, df303, df310, df311, df312, df313, df320, df321, df322, df323, df330, df331, df332, df333;

	// cubic_factor set -0.5 for test
	DOUBLE cubic_factor = -0.5;

	// f2d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero...

	// Use the inverse matrix
    if (inv)
    	Ainv = A;
    else
    	Ainv = A.transpose();

    // Go from the 2D slice coordinates to the 3D coordinates
    Ainv *= (DOUBLE)padding_factor;  // take scaling into account directly
    int max_r2 = r_max * r_max;
    int min_r2_nn = r_min_nn * r_min_nn;

#define DEBUG_BACKP
#ifdef DEBUG_BACKP
    std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
    std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
    std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
    std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
    std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
    std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
    std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
	std::cerr << " Interpolator= " << interpolator << std::endl;
#endif

    for (int i=0; i < YSIZE(f2d); i++)
	{
		// Dont search beyond square with side max_r
		if (i <= r_max)
		{
			y = i;
			first_x = 0;
		}
		else if (i >= YSIZE(f2d) - r_max)
		{
			y = i - YSIZE(f2d);
			// x==0 plane is stored twice in the FFTW format. Dont set it twice in BACKPROJECTION!
			first_x = 1;
		}
		else
			continue;

		y2 = y * y;
		for (int x=first_x; x <= r_max; x++)
		{
	    	// Only include points with radius < max_r (exclude points outside circle in square)
			r2 = x * x + y2;
			if (r2 > max_r2)
				continue;

			// Get the relevant value in the input image
			my_val = DIRECT_A2D_ELEM(f2d, i, x);

			// Get the weight
			if (Mweight != NULL)
				my_weight = DIRECT_A2D_ELEM(*Mweight, i, x);
			// else: my_weight was already initialised to 1.

			if (my_weight > 0.)
			{

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

					mfx = 1. - fx;
					mfy = 1. - fy;
					mfz = 1. - fz;

					dd000 = mfz * mfy * mfx;
					dd001 = mfz * mfy *  fx;
					dd010 = mfz *  fy * mfx;
					dd011 = mfz *  fy *  fx;
					dd100 =  fz * mfy * mfx;
					dd101 =  fz * mfy *  fx;
					dd110 =  fz *  fy * mfx;
					dd111 =  fz *  fy *  fx;

					if (is_neg_x)
						my_val = conj(my_val);

					// Store slice in 3D weighted sum
					DIRECT_A3D_ELEM(data, z0, y0, x0) += dd000 * my_val;
					DIRECT_A3D_ELEM(data, z0, y0, x1) += dd001 * my_val;
					DIRECT_A3D_ELEM(data, z0, y1, x0) += dd010 * my_val;
					DIRECT_A3D_ELEM(data, z0, y1, x1) += dd011 * my_val;
					DIRECT_A3D_ELEM(data, z1, y0, x0) += dd100 * my_val;
					DIRECT_A3D_ELEM(data, z1, y0, x1) += dd101 * my_val;
					DIRECT_A3D_ELEM(data, z1, y1, x0) += dd110 * my_val;
					DIRECT_A3D_ELEM(data, z1, y1, x1) += dd111 * my_val;
					// Store corresponding weights
					DIRECT_A3D_ELEM(weight, z0, y0, x0) += dd000 * my_weight;
					DIRECT_A3D_ELEM(weight, z0, y0, x1) += dd001 * my_weight;
					DIRECT_A3D_ELEM(weight, z0, y1, x0) += dd010 * my_weight;
					DIRECT_A3D_ELEM(weight, z0, y1, x1) += dd011 * my_weight;
					DIRECT_A3D_ELEM(weight, z1, y0, x0) += dd100 * my_weight;
					DIRECT_A3D_ELEM(weight, z1, y0, x1) += dd101 * my_weight;
					DIRECT_A3D_ELEM(weight, z1, y1, x0) += dd110 * my_weight;
					DIRECT_A3D_ELEM(weight, z1, y1, x1) += dd111 * my_weight;

				} // endif TRILINEAR
				else if (interpolator == NEAREST_NEIGHBOUR )
				{

					x0 = ROUND(xp);
					y0 = ROUND(yp);
					z0 = ROUND(zp);

					if (x0 < 0)
					{
						A3D_ELEM(data, -z0, -y0, -x0) += conj(my_val);
						A3D_ELEM(weight, -z0, -y0, -x0) += my_weight;
					}
					else
					{
						A3D_ELEM(data, z0, y0, x0) += my_val;
						A3D_ELEM(weight, z0, y0, x0) += my_weight;
					}

				} // endif NEAREST_NEIGHBOUR
				else if (interpolator == CUBIC )
				{
					// std::cerr << " Entering cubic interpolation " << std::endl;

					if (xp <0)
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

#ifdef DEBUG_BACKP
    std::cerr << " u0= "<< u0 << std::endl;
    std::cerr << " u1= "<< u1 << std::endl;
    std::cerr << " u2= "<< u2 << std::endl;
    std::cerr << " u3= "<< u3 << std::endl;
	std::cerr << " v0= "<< v0 << std::endl;
    std::cerr << " v1= "<< v1 << std::endl;
    std::cerr << " v2= "<< v2 << std::endl;
    std::cerr << " v3= "<< v3 << std::endl;
	std::cerr << " w0= "<< w0 << std::endl;
    std::cerr << " w1= "<< w1 << std::endl;
    std::cerr << " w2= "<< w2 << std::endl;
    std::cerr << " w3= "<< w3 << std::endl;
#endif

					// std::cerr << " intXX to w3 value equal to " << intXX << std::endl << intYY << std::endl << intZZ << std::endl;
					// std::cerr << u0 << " " << u1 << " " << u2 << " " << u3 << " " << v0 << " " << v1 << " " << v2 << " " << v3 << 
					//			" " << w0 << " " << w1 << " " << w2 << " " << w3 << std::endl;

					df000 = u1 * u2 * u3 * v1 * v2 * v3 * w1 * w2 * w3;
					df001 = u0 * u2 * u3 * v1 * v2 * v3 * w1 * w2 * w3;
					df002 = u1 * u0 * u3 * v1 * v2 * v3 * w1 * w2 * w3;
					df003 = u1 * u2 * u0 * v1 * v2 * v3 * w1 * w2 * w3;
					df010 = u1 * u2 * u3 * v0 * v2 * v3 * w1 * w2 * w3;
					df011 = u0 * u2 * u3 * v0 * v2 * v3 * w1 * w2 * w3;
					df012 = u1 * u0 * u3 * v0 * v2 * v3 * w1 * w2 * w3;
					df013 = u1 * u2 * u0 * v0 * v2 * v3 * w1 * w2 * w3;
					df020 = u1 * u2 * u3 * v1 * v0 * v3 * w1 * w2 * w3;
					df021 = u0 * u2 * u3 * v1 * v0 * v3 * w1 * w2 * w3;
					df022 = u1 * u0 * u3 * v1 * v0 * v3 * w1 * w2 * w3;
					df023 = u1 * u2 * u0 * v1 * v0 * v3 * w1 * w2 * w3;
					df030 = u1 * u2 * u3 * v1 * v2 * v0 * w1 * w2 * w3;
					df031 = u0 * u2 * u3 * v1 * v2 * v0 * w1 * w2 * w3;
					df032 = u1 * u0 * u3 * v1 * v2 * v0 * w1 * w2 * w3;
					df033 = u1 * u2 * u0 * v1 * v2 * v0 * w1 * w2 * w3;
					df100 = u1 * u2 * u3 * v1 * v2 * v3 * w0 * w2 * w3;
					df101 = u0 * u2 * u3 * v1 * v2 * v3 * w0 * w2 * w3;
					df102 = u1 * u0 * u3 * v1 * v2 * v3 * w0 * w2 * w3;
					df103 = u1 * u2 * u0 * v1 * v2 * v3 * w0 * w2 * w3;
					df110 = u1 * u2 * u3 * v0 * v2 * v3 * w0 * w2 * w3;
					df111 = u0 * u2 * u3 * v0 * v2 * v3 * w0 * w2 * w3;
					df112 = u1 * u0 * u3 * v0 * v2 * v3 * w0 * w2 * w3;
					df113 = u1 * u2 * u0 * v0 * v2 * v3 * w0 * w2 * w3;
					df120 = u1 * u2 * u3 * v1 * v0 * v3 * w0 * w2 * w3;
					df121 = u0 * u2 * u3 * v1 * v0 * v3 * w0 * w2 * w3;
					df122 = u1 * u0 * u3 * v1 * v0 * v3 * w0 * w2 * w3;
					df123 = u1 * u2 * u0 * v1 * v0 * v3 * w0 * w2 * w3;
					df130 = u1 * u2 * u3 * v1 * v2 * v0 * w0 * w2 * w3;
					df131 = u0 * u2 * u3 * v1 * v2 * v0 * w0 * w2 * w3;
					df132 = u1 * u0 * u3 * v1 * v2 * v0 * w0 * w2 * w3;
					df133 = u1 * u2 * u0 * v1 * v2 * v0 * w0 * w2 * w3;
					df200 = u1 * u2 * u3 * v1 * v2 * v3 * w1 * w0 * w3;
					df201 = u0 * u2 * u3 * v1 * v2 * v3 * w1 * w0 * w3;
					df202 = u1 * u0 * u3 * v1 * v2 * v3 * w1 * w0 * w3;
					df203 = u1 * u2 * u0 * v1 * v2 * v3 * w1 * w0 * w3;
					df210 = u1 * u2 * u3 * v0 * v2 * v3 * w1 * w0 * w3;
					df211 = u0 * u2 * u3 * v0 * v2 * v3 * w1 * w0 * w3;
					df212 = u1 * u0 * u3 * v0 * v2 * v3 * w1 * w0 * w3;
					df213 = u1 * u2 * u0 * v0 * v2 * v3 * w1 * w0 * w3;
					df220 = u1 * u2 * u3 * v1 * v0 * v3 * w1 * w0 * w3;
					df221 = u0 * u2 * u3 * v1 * v0 * v3 * w1 * w0 * w3;
					df222 = u1 * u0 * u3 * v1 * v0 * v3 * w1 * w0 * w3;
					df223 = u1 * u2 * u0 * v1 * v0 * v3 * w1 * w0 * w3;
					df230 = u1 * u2 * u3 * v1 * v2 * v0 * w1 * w0 * w3;
					df231 = u0 * u2 * u3 * v1 * v2 * v0 * w1 * w0 * w3;
					df232 = u1 * u0 * u3 * v1 * v2 * v0 * w1 * w0 * w3;
					df233 = u1 * u2 * u0 * v1 * v2 * v0 * w1 * w0 * w3;
					df300 = u1 * u2 * u3 * v1 * v2 * v3 * w1 * w2 * w0;
					df301 = u0 * u2 * u3 * v1 * v2 * v3 * w1 * w2 * w0;
					df302 = u1 * u0 * u3 * v1 * v2 * v3 * w1 * w2 * w0;
					df303 = u1 * u2 * u0 * v1 * v2 * v3 * w1 * w2 * w0;
					df310 = u1 * u2 * u3 * v0 * v2 * v3 * w1 * w2 * w0;
					df311 = u0 * u2 * u3 * v0 * v2 * v3 * w1 * w2 * w0;
					df312 = u1 * u0 * u3 * v0 * v2 * v3 * w1 * w2 * w0;
					df313 = u1 * u2 * u0 * v0 * v2 * v3 * w1 * w2 * w0;
					df320 = u1 * u2 * u3 * v1 * v0 * v3 * w1 * w2 * w0;
					df321 = u0 * u2 * u3 * v1 * v0 * v3 * w1 * w2 * w0;
					df322 = u1 * u0 * u3 * v1 * v0 * v3 * w1 * w2 * w0;
					df323 = u1 * u2 * u0 * v1 * v0 * v3 * w1 * w2 * w0;
					df330 = u1 * u2 * u3 * v1 * v2 * v0 * w1 * w2 * w0;
					df331 = u0 * u2 * u3 * v1 * v2 * v0 * w1 * w2 * w0;
					df332 = u1 * u0 * u3 * v1 * v2 * v0 * w1 * w2 * w0;
					df333 = u1 * u2 * u0 * v1 * v2 * v0 * w1 * w2 * w0;

#ifdef DEBUG_BACKP
    std::cerr << " df000= "<< df000 << std::endl;
    std::cerr << " df001= "<< df001 << std::endl;
    std::cerr << " df002= "<< df002 << std::endl;
    std::cerr << " df003= "<< df003 << std::endl;
	std::cerr << " df010= "<< df010 << std::endl;
    std::cerr << " df011= "<< df011 << std::endl;
    std::cerr << " df012= "<< df012 << std::endl;
    std::cerr << " df013= "<< df013 << std::endl;
	std::cerr << " df020= "<< df020 << std::endl;
    std::cerr << " df021= "<< df021 << std::endl;
    std::cerr << " df022= "<< df022 << std::endl;
    std::cerr << " df023= "<< df023 << std::endl;
	std::cerr << " df030= "<< df030 << std::endl;
    std::cerr << " df031= "<< df031 << std::endl;
    std::cerr << " df032= "<< df032 << std::endl;
    std::cerr << " df033= "<< df033 << std::endl;
#endif

					if (is_neg_x)
						my_val = conj(my_val);

					// std::cerr << " Before changing the matrix data " << std::endl;
					// std::cerr << " Origin data matrix is " << data << std::endl;

					// Store slice in 3D weighted sum
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX - 1) += df000 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX ) += df001 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX + 1) += df002 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX + 2) += df003 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX - 1) += df010 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX ) += df011 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX + 1) += df012 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX + 2) += df013 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX - 1) += df020 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX ) += df021 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX + 1) += df022 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX + 2) += df023 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX - 1) += df030 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX ) += df031 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX + 1) += df032 * my_val;
					DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX + 2) += df033 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX - 1) += df100 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX ) += df101 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX + 1) += df102 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX + 2) += df103 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY , intXX - 1) += df110 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY , intXX ) += df111 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY , intXX + 1) += df112 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY , intXX + 2) += df113 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX - 1) += df120 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX ) += df121 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX + 1) += df122 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX + 2) += df123 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX - 1) += df130 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX ) += df131 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX + 1) += df132 * my_val;
					DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX + 2) += df133 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX - 1) += df200 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX ) += df201 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX + 1) += df202 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX + 2) += df203 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX - 1) += df210 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX ) += df211 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX + 1) += df212 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX + 2) += df213 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX - 1) += df220 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX ) += df221 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX + 1) += df222 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX + 2) += df223 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX - 1) += df230 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX ) += df231 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX + 1) += df232 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX + 2) += df233 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX - 1) += df300 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX ) += df301 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX + 1) += df302 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX + 2) += df303 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX - 1) += df310 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX ) += df311 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX + 1) += df312 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX + 2) += df313 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX - 1) += df320 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX ) += df321 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX + 1) += df322 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX + 2) += df323 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX - 1) += df330 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX ) += df331 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX + 1) += df332 * my_val;
					DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX + 2) += df333 * my_val;

					// std::cerr << " After changing the data matrix " << std::endl;
					// std::cerr << " The new data matrix is " << data << std::endl;

					// std::cerr << " Before changing the weight " << std::endl;
					// std::cerr << " The origin weight matrix is " << weight << std::endl;

					// Store corresponding weights

#ifdef DEBUG_BACKP_NN  // test for what will happen if weight use NN interpolator
					x0 = ROUND(xp);
					y0 = ROUND(yp);
					z0 = ROUND(zp);
					A3D_ELEM(weight, z0, y0, x0) += my_weight;
#endif


// #ifdef DEBUG_ZSY_WEIGHT
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY - 1, intXX - 1) += df000 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY - 1, intXX ) += df001 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY - 1, intXX + 1) += df002 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY - 1, intXX + 2) += df003 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY , intXX - 1) += df010 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY , intXX ) += df011 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY , intXX + 1) += df012 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY , intXX + 2) += df013 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 1, intXX - 1) += df020 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 1, intXX ) += df021 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 1, intXX + 1) += df022 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 1, intXX + 2) += df023 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 2, intXX - 1) += df030 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 2, intXX ) += df031 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 2, intXX + 1) += df032 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 2, intXX + 2) += df033 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY - 1, intXX - 1) += df100 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY - 1, intXX ) += df101 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY - 1, intXX + 1) += df102 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY - 1, intXX + 2) += df103 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY , intXX - 1) += df110 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY , intXX ) += df111 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY , intXX + 1) += df112 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY , intXX + 2) += df113 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY + 1, intXX - 1) += df120 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY + 1, intXX ) += df121 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY + 1, intXX + 1) += df122 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY + 1, intXX + 2) += df123 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY + 2, intXX - 1) += df130 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY + 2, intXX ) += df131 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY + 2, intXX + 1) += df132 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ , intYY + 2, intXX + 2) += df133 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY - 1, intXX - 1) += df200 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY - 1, intXX ) += df201 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY - 1, intXX + 1) += df202 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY - 1, intXX + 2) += df203 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY , intXX - 1) += df210 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY , intXX ) += df211 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY , intXX + 1) += df212 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY , intXX + 2) += df213 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 1, intXX - 1) += df220 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 1, intXX ) += df221 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 1, intXX + 1) += df222 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 1, intXX + 2) += df223 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 2, intXX - 1) += df230 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 2, intXX ) += df231 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 2, intXX + 1) += df232 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 2, intXX + 2) += df233 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY - 1, intXX - 1) += df300 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY - 1, intXX ) += df301 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY - 1, intXX + 1) += df302 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY - 1, intXX + 2) += df303 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY , intXX - 1) += df310 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY , intXX ) += df311 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY , intXX + 1) += df312 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY , intXX + 2) += df313 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 1, intXX - 1) += df320 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 1, intXX ) += df321 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 1, intXX + 1) += df322 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 1, intXX + 2) += df323 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 2, intXX - 1) += df330 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 2, intXX ) += df331 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 2, intXX + 1) += df332 * my_weight;
					DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 2, intXX + 2) += df333 * my_weight;
// #endif
					// std::cerr << " After changing the weight matrix " << std::endl;
					// std::cerr << " The new weight matrix is " << weight << std::endl;

				} // endif CUBIC
				else
				{
					REPORT_ERROR("FourierInterpolator::backproject%%ERROR: unrecognized interpolator ");
				}
			} // endif weight>0.
		} // endif x-loop
	} // endif y-loop

	std::cerr << " Leaving backprojection " << std::endl;

}

void BackProjector::backrotate2D(const MultidimArray<Complex > &f2d,
		                         const Matrix2D<DOUBLE> &A, bool inv,
		                         const MultidimArray<DOUBLE> *Mweight)
{
	DOUBLE fx, fy, mfx, mfy, xp, yp;
	int first_x, x0, x1, y0, y1, y, y2, r2;
	bool is_neg_x;
	DOUBLE dd00, dd01, dd10, dd11;
	Complex my_val;
	Matrix2D<DOUBLE> Ainv;
	DOUBLE my_weight = 1.;

	int intXX, intYY;
	DOUBLE u0, u1, u2, u3, v0, v1, v2, v3;
	DOUBLE df00, df01, df02, df03, df10, df11, df12, df13, df20, df21, df22, df23, df30, df31, df32, df33;

	// cubic_factor set -0.5 for test
	DOUBLE cubic_factor = -0.5;

	// f2d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero...

	// Use the inverse matrix
    if (inv)
    	Ainv = A;
    else
    	Ainv = A.transpose();

    // Go from the 2D slice coordinates to the data-array coordinates
    Ainv *= (DOUBLE)padding_factor;  // take scaling into account directly
    int max_r2 = r_max * r_max;
    int min_r2_nn = r_min_nn * r_min_nn;

//#define DEBUG_BACKROTATE
#ifdef DEBUG_BACKROTATE
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
		if (i <= r_max)
		{
			y = i;
			first_x = 0;
		}
		else if (i >= YSIZE(f2d) - r_max)
		{
			y = i - YSIZE(f2d);
			// x==0 plane is stored twice in the FFTW format. Dont set it twice in BACKPROJECTION!
			first_x = 1;
		}
		else
			continue;

		y2 = y * y;
		for (int x=first_x; x <= r_max; x++)
		{
	    	// Only include points with radius < max_r (exclude points outside circle in square)
			r2 = x * x + y2;
			if (r2 > max_r2)
				continue;

			// Get the relevant value in the input image
			my_val = DIRECT_A2D_ELEM(f2d, i, x);

			// Get the weight
			if (Mweight != NULL)
				my_weight = DIRECT_A2D_ELEM(*Mweight, i, x);
			// else: my_weight was already initialised to 1.

			if (my_weight > 0.)
			{
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
					// In that way use DIRECT_A2D_ELEM, rather than A2D_ELEM
					x0 = FLOOR(xp);
					fx = xp - x0;
					x1 = x0 + 1;

					y0 = FLOOR(yp);
					fy = yp - y0;
					y0 -=  STARTINGY(data);
					y1 = y0 + 1;

					mfx = 1. - fx;
					mfy = 1. - fy;

					dd00 = mfy * mfx;
					dd01 = mfy *  fx;
					dd10 =  fy * mfx;
					dd11 =  fy *  fx;

					if (is_neg_x)
						my_val = conj(my_val);

					// Store slice in 3D weighted sum
					DIRECT_A2D_ELEM(data, y0, x0) += dd00 * my_val;
					DIRECT_A2D_ELEM(data, y0, x1) += dd01 * my_val;
					DIRECT_A2D_ELEM(data, y1, x0) += dd10 * my_val;
					DIRECT_A2D_ELEM(data, y1, x1) += dd11 * my_val;

					// Store corresponding weights
					DIRECT_A2D_ELEM(weight, y0, x0) += dd00 * my_weight;
					DIRECT_A2D_ELEM(weight, y0, x1) += dd01 * my_weight;
					DIRECT_A2D_ELEM(weight, y1, x0) += dd10 * my_weight;
					DIRECT_A2D_ELEM(weight, y1, x1) += dd11 * my_weight;

				} // endif TRILINEAR
				else if (interpolator == NEAREST_NEIGHBOUR )
				{
					x0 = ROUND(xp);
					y0 = ROUND(yp);
					if (x0 < 0)
					{
						A2D_ELEM(data, -y0, -x0) += conj(my_val);
						A2D_ELEM(weight, -y0, -x0) += my_weight;
					}
					else
					{
						A2D_ELEM(data, y0, x0) += my_val;
						A2D_ELEM(weight, y0, x0) += my_weight;
					}
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

					df00 = u1 * u2 * u3 * v1 * v2 * v3;
					df01 = u0 * u2 * u3 * v1 * v2 * v3;
					df02 = u1 * u0 * u3 * v1 * v2 * v3;
					df03 = u1 * u2 * u0 * v1 * v2 * v3;
					df10 = u1 * u2 * u3 * v0 * v2 * v3;
					df11 = u0 * u2 * u3 * v0 * v2 * v3;
					df12 = u1 * u0 * u3 * v0 * v2 * v3;
					df13 = u1 * u2 * u0 * v0 * v2 * v3;
					df20 = u1 * u2 * u3 * v1 * v0 * v3;
					df21 = u0 * u2 * u3 * v1 * v0 * v3;
					df22 = u1 * u0 * u3 * v1 * v0 * v3;
					df23 = u1 * u2 * u0 * v1 * v0 * v3;
					df30 = u1 * u2 * u3 * v1 * v2 * v0;
					df31 = u0 * u2 * u3 * v1 * v2 * v0;
					df32 = u1 * u0 * u3 * v1 * v2 * v0;
					df33 = u1 * u2 * u0 * v1 * v2 * v0;

					if (is_neg_x)
						my_val = conj(my_val);

					// Store slice in 3D weighted sum
					DIRECT_A2D_ELEM(data, intYY - 1, intXX - 1) += df00 * my_val;
					DIRECT_A2D_ELEM(data, intYY - 1, intXX ) += df01 * my_val;
					DIRECT_A2D_ELEM(data, intYY - 1, intXX + 1) += df02 * my_val;
					DIRECT_A2D_ELEM(data, intYY - 1, intXX + 2) += df03 * my_val;
					DIRECT_A2D_ELEM(data, intYY , intXX - 1) += df10 * my_val;
					DIRECT_A2D_ELEM(data, intYY , intXX ) += df11 * my_val;
					DIRECT_A2D_ELEM(data, intYY , intXX + 1) += df12 * my_val;
					DIRECT_A2D_ELEM(data, intYY , intXX + 2) += df13 * my_val;
					DIRECT_A2D_ELEM(data, intYY + 1, intXX - 1) += df20 * my_val;
					DIRECT_A2D_ELEM(data, intYY + 1, intXX ) += df21 * my_val;
					DIRECT_A2D_ELEM(data, intYY + 1, intXX + 1) += df22 * my_val;
					DIRECT_A2D_ELEM(data, intYY + 1, intXX + 2) += df23 * my_val;
					DIRECT_A2D_ELEM(data, intYY + 2, intXX - 1) += df30 * my_val;
					DIRECT_A2D_ELEM(data, intYY + 2, intXX ) += df31 * my_val;
					DIRECT_A2D_ELEM(data, intYY + 2, intXX + 1) += df32 * my_val;
					DIRECT_A2D_ELEM(data, intYY + 2, intXX + 2) += df33 * my_val;

					// Store corresponding weights
					DIRECT_A2D_ELEM(weight, intYY - 1, intXX - 1) += df00 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY - 1, intXX ) += df01 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY - 1, intXX + 1) += df02 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY - 1, intXX + 2) += df03 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY , intXX - 1) += df10 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY , intXX ) += df11 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY , intXX + 1) += df12 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY , intXX + 2) += df13 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY + 1, intXX - 1) += df20 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY + 1, intXX ) += df21 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY + 1, intXX + 1) += df22 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY + 1, intXX + 2) += df23 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY + 2, intXX - 1) += df30 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY + 2, intXX ) += df31 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY + 2, intXX + 1) += df32 * my_weight;
					DIRECT_A2D_ELEM(weight, intYY + 2, intXX + 2) += df33 * my_weight;

				} // endif CUBIC
				else
				{
					REPORT_ERROR("FourierInterpolator::backrotate2D%%ERROR: unrecognized interpolator ");
				}
			} // endif weight > 0.
		} // endif x-loop
	} // endif y-loop
}

void BackProjector::backrotate3D(const MultidimArray<Complex > &f3d,
		                         const Matrix2D<DOUBLE> &A, bool inv,
		                         const MultidimArray<DOUBLE> *Mweight)
{
	DOUBLE fx, fy, fz, mfx, mfy, mfz, xp, yp, zp;
	int first_x, x0, x1, y0, y1, z0, z1, y, y2, z, z2, r2;
	bool is_neg_x;
	DOUBLE dd000, dd010, dd100, dd110, dd001, dd011, dd101, dd111;
	Complex my_val;
	Matrix2D<DOUBLE> Ainv;
	DOUBLE my_weight = 1.;

	int intXX, intYY, intZZ;
	DOUBLE u0, u1, u2, u3, v0, v1, v2, v3, w0, w1, w2, w3;
	DOUBLE df000, df001, df002, df003, df010, df011, df012, df013, df020, df021, df022, df023, df030, df031, df032, df033,
		   df100, df101, df102, df103, df110, df111, df112, df113, df120, df121, df122, df123, df130, df131, df132, df133,
		   df200, df201, df202, df203, df210, df211, df212, df213, df220, df221, df222, df223, df230, df231, df232, df233,
		   df300, df301, df302, df303, df310, df311, df312, df313, df320, df321, df322, df323, df330, df331, df332, df333;

	// cubic_factor set -0.5 for test
	DOUBLE cubic_factor = -0.5;

	// f3d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero...

	// Use the inverse matrix
    if (inv)
    	Ainv = A;
    else
    	Ainv = A.transpose();

    // Go from the 2D slice coordinates to the data-array coordinates
    Ainv *= (DOUBLE)padding_factor;  // take scaling into account directly
    int max_r2 = r_max * r_max;
    int min_r2_nn = r_min_nn * r_min_nn;

//#define DEBUG_BACKROTATE
#ifdef DEBUG_BACKROTATE
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

    for (int k=0; k < ZSIZE(f3d); k++)
	{
		// Don't search beyond square with side max_r
		if (k <= r_max)
		{
			z = k;
			first_x = 0;
		}
		else if (k >= YSIZE(f3d) - r_max)
		{
			z = k - YSIZE(f3d);
			/// TODO: still check this better in the 3D case!!!
			// x==0 (y,z)-plane is stored twice in the FFTW format. Don't set it twice in BACKPROJECTION!
			first_x = 1;
		}
		else
			continue;

		z2 = z * z;
		for (int i=0; i < YSIZE(f3d); i++)
		{
			// Don't search beyond square with side max_r
			if (i <= r_max)
			{
				y = i;
			}
			else if (i >= YSIZE(f3d) - r_max)
			{
				y = i - YSIZE(f3d);
			}
			else
				continue;

			y2 = y * y;
			for (int x = first_x; x <= r_max; x++)
			{
				// Only include points with radius < max_r (exclude points outside circle in square)
				r2 = x * x + y2 + z2;
				if (r2 > max_r2)
					continue;

				// Get the relevant value in the input image
				my_val = DIRECT_A3D_ELEM(f3d, k, i, x);

				// Get the weight
				if (Mweight != NULL)
					my_weight = DIRECT_A3D_ELEM(*Mweight, k, i, x);
				// else: my_weight was already initialised to 1.

				if (my_weight > 0.)
				{
					// Get logical coordinates in the 3D map
					xp = Ainv(0,0) * x + Ainv(0,1) * y + Ainv(0,2) * z;
					yp = Ainv(1,0) * x + Ainv(1,1) * y + Ainv(1,2) * z;
					zp = Ainv(2,0) * x + Ainv(2,1) * y + Ainv(2,2) * z;

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

						mfx = 1. - fx;
						mfy = 1. - fy;
						mfz = 1. - fz;

						dd000 = mfz * mfy * mfx;
						dd001 = mfz * mfy *  fx;
						dd010 = mfz *  fy * mfx;
						dd011 = mfz *  fy *  fx;
						dd100 =  fz * mfy * mfx;
						dd101 =  fz * mfy *  fx;
						dd110 =  fz *  fy * mfx;
						dd111 =  fz *  fy *  fx;

						if (is_neg_x)
							my_val = conj(my_val);

						// Store slice in 3D weighted sum
						DIRECT_A3D_ELEM(data, z0, y0, x0) += dd000 * my_val;
						DIRECT_A3D_ELEM(data, z0, y0, x1) += dd001 * my_val;
						DIRECT_A3D_ELEM(data, z0, y1, x0) += dd010 * my_val;
						DIRECT_A3D_ELEM(data, z0, y1, x1) += dd011 * my_val;
						DIRECT_A3D_ELEM(data, z1, y0, x0) += dd100 * my_val;
						DIRECT_A3D_ELEM(data, z1, y0, x1) += dd101 * my_val;
						DIRECT_A3D_ELEM(data, z1, y1, x0) += dd110 * my_val;
						DIRECT_A3D_ELEM(data, z1, y1, x1) += dd111 * my_val;
						// Store corresponding weights
						DIRECT_A3D_ELEM(weight, z0, y0, x0) += dd000 * my_weight;
						DIRECT_A3D_ELEM(weight, z0, y0, x1) += dd001 * my_weight;
						DIRECT_A3D_ELEM(weight, z0, y1, x0) += dd010 * my_weight;
						DIRECT_A3D_ELEM(weight, z0, y1, x1) += dd011 * my_weight;
						DIRECT_A3D_ELEM(weight, z1, y0, x0) += dd100 * my_weight;
						DIRECT_A3D_ELEM(weight, z1, y0, x1) += dd101 * my_weight;
						DIRECT_A3D_ELEM(weight, z1, y1, x0) += dd110 * my_weight;
						DIRECT_A3D_ELEM(weight, z1, y1, x1) += dd111 * my_weight;


					} // endif TRILINEAR
					else if (interpolator == NEAREST_NEIGHBOUR )
					{
						x0 = ROUND(xp);
						y0 = ROUND(yp);
						z0 = ROUND(zp);

						if (x0 < 0)
						{
							A3D_ELEM(data, -z0, -y0, -x0) += conj(my_val);
							A3D_ELEM(weight, -z0, -y0, -x0) += my_weight;
						}
						else
						{
							A3D_ELEM(data, z0, y0, x0) += my_val;
							A3D_ELEM(weight, z0, y0, x0) += my_weight;
						}

					} // endif NEAREST_NEIGHBOUR
					else if (interpolator == CUBIC )
					{
						if (xp <0)
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

						df000 = u1 * u2 * u3 * v1 * v2 * v3 * w1 * w2 * w3;
						df001 = u0 * u2 * u3 * v1 * v2 * v3 * w1 * w2 * w3;
						df002 = u1 * u0 * u3 * v1 * v2 * v3 * w1 * w2 * w3;
						df003 = u1 * u2 * u0 * v1 * v2 * v3 * w1 * w2 * w3;
						df010 = u1 * u2 * u3 * v0 * v2 * v3 * w1 * w2 * w3;
						df011 = u0 * u2 * u3 * v0 * v2 * v3 * w1 * w2 * w3;
						df012 = u1 * u0 * u3 * v0 * v2 * v3 * w1 * w2 * w3;
						df013 = u1 * u2 * u0 * v0 * v2 * v3 * w1 * w2 * w3;
						df020 = u1 * u2 * u3 * v1 * v0 * v3 * w1 * w2 * w3;
						df021 = u0 * u2 * u3 * v1 * v0 * v3 * w1 * w2 * w3;
						df022 = u1 * u0 * u3 * v1 * v0 * v3 * w1 * w2 * w3;
						df023 = u1 * u2 * u0 * v1 * v0 * v3 * w1 * w2 * w3;
						df030 = u1 * u2 * u3 * v1 * v2 * v0 * w1 * w2 * w3;
						df031 = u0 * u2 * u3 * v1 * v2 * v0 * w1 * w2 * w3;
						df032 = u1 * u0 * u3 * v1 * v2 * v0 * w1 * w2 * w3;
						df033 = u1 * u2 * u0 * v1 * v2 * v0 * w1 * w2 * w3;
						df100 = u1 * u2 * u3 * v1 * v2 * v3 * w0 * w2 * w3;
						df101 = u0 * u2 * u3 * v1 * v2 * v3 * w0 * w2 * w3;
						df102 = u1 * u0 * u3 * v1 * v2 * v3 * w0 * w2 * w3;
						df103 = u1 * u2 * u0 * v1 * v2 * v3 * w0 * w2 * w3;
						df110 = u1 * u2 * u3 * v0 * v2 * v3 * w0 * w2 * w3;
						df111 = u0 * u2 * u3 * v0 * v2 * v3 * w0 * w2 * w3;
						df112 = u1 * u0 * u3 * v0 * v2 * v3 * w0 * w2 * w3;
						df113 = u1 * u2 * u0 * v0 * v2 * v3 * w0 * w2 * w3;
						df120 = u1 * u2 * u3 * v1 * v0 * v3 * w0 * w2 * w3;
						df121 = u0 * u2 * u3 * v1 * v0 * v3 * w0 * w2 * w3;
						df122 = u1 * u0 * u3 * v1 * v0 * v3 * w0 * w2 * w3;
						df123 = u1 * u2 * u0 * v1 * v0 * v3 * w0 * w2 * w3;
						df130 = u1 * u2 * u3 * v1 * v2 * v0 * w0 * w2 * w3;
						df131 = u0 * u2 * u3 * v1 * v2 * v0 * w0 * w2 * w3;
						df132 = u1 * u0 * u3 * v1 * v2 * v0 * w0 * w2 * w3;
						df133 = u1 * u2 * u0 * v1 * v2 * v0 * w0 * w2 * w3;
						df200 = u1 * u2 * u3 * v1 * v2 * v3 * w1 * w0 * w3;
						df201 = u0 * u2 * u3 * v1 * v2 * v3 * w1 * w0 * w3;
						df202 = u1 * u0 * u3 * v1 * v2 * v3 * w1 * w0 * w3;
						df203 = u1 * u2 * u0 * v1 * v2 * v3 * w1 * w0 * w3;
						df210 = u1 * u2 * u3 * v0 * v2 * v3 * w1 * w0 * w3;
						df211 = u0 * u2 * u3 * v0 * v2 * v3 * w1 * w0 * w3;
						df212 = u1 * u0 * u3 * v0 * v2 * v3 * w1 * w0 * w3;
						df213 = u1 * u2 * u0 * v0 * v2 * v3 * w1 * w0 * w3;
						df220 = u1 * u2 * u3 * v1 * v0 * v3 * w1 * w0 * w3;
						df221 = u0 * u2 * u3 * v1 * v0 * v3 * w1 * w0 * w3;
						df222 = u1 * u0 * u3 * v1 * v0 * v3 * w1 * w0 * w3;
						df223 = u1 * u2 * u0 * v1 * v0 * v3 * w1 * w0 * w3;
						df230 = u1 * u2 * u3 * v1 * v2 * v0 * w1 * w0 * w3;
						df231 = u0 * u2 * u3 * v1 * v2 * v0 * w1 * w0 * w3;
						df232 = u1 * u0 * u3 * v1 * v2 * v0 * w1 * w0 * w3;
						df233 = u1 * u2 * u0 * v1 * v2 * v0 * w1 * w0 * w3;
						df300 = u1 * u2 * u3 * v1 * v2 * v3 * w1 * w2 * w0;
						df301 = u0 * u2 * u3 * v1 * v2 * v3 * w1 * w2 * w0;
						df302 = u1 * u0 * u3 * v1 * v2 * v3 * w1 * w2 * w0;
						df303 = u1 * u2 * u0 * v1 * v2 * v3 * w1 * w2 * w0;
						df310 = u1 * u2 * u3 * v0 * v2 * v3 * w1 * w2 * w0;
						df311 = u0 * u2 * u3 * v0 * v2 * v3 * w1 * w2 * w0;
						df312 = u1 * u0 * u3 * v0 * v2 * v3 * w1 * w2 * w0;
						df313 = u1 * u2 * u0 * v0 * v2 * v3 * w1 * w2 * w0;
						df320 = u1 * u2 * u3 * v1 * v0 * v3 * w1 * w2 * w0;
						df321 = u0 * u2 * u3 * v1 * v0 * v3 * w1 * w2 * w0;
						df322 = u1 * u0 * u3 * v1 * v0 * v3 * w1 * w2 * w0;
						df323 = u1 * u2 * u0 * v1 * v0 * v3 * w1 * w2 * w0;
						df330 = u1 * u2 * u3 * v1 * v2 * v0 * w1 * w2 * w0;
						df331 = u0 * u2 * u3 * v1 * v2 * v0 * w1 * w2 * w0;
						df332 = u1 * u0 * u3 * v1 * v2 * v0 * w1 * w2 * w0;
						df333 = u1 * u2 * u0 * v1 * v2 * v0 * w1 * w2 * w0;
						
						if (is_neg_x)
							my_val = conj(my_val);

						// Store slice in 3D weighted sum
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX - 1) += df000 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX ) += df001 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX + 1) += df002 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY - 1, intXX + 2) += df003 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX - 1) += df010 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX ) += df011 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX + 1) += df012 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY , intXX + 2) += df013 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX - 1) += df020 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX ) += df021 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX + 1) += df022 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 1, intXX + 2) += df023 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX - 1) += df030 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX ) += df031 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX + 1) += df032 * my_val;
						DIRECT_A3D_ELEM(data, intZZ - 1, intYY + 2, intXX + 2) += df033 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX - 1) += df100 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX ) += df101 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX + 1) += df102 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY - 1, intXX + 2) += df103 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY , intXX - 1) += df110 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY , intXX ) += df111 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY , intXX + 1) += df112 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY , intXX + 2) += df113 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX - 1) += df120 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX ) += df121 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX + 1) += df122 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY + 1, intXX + 2) += df123 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX - 1) += df130 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX ) += df131 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX + 1) += df132 * my_val;
						DIRECT_A3D_ELEM(data, intZZ , intYY + 2, intXX + 2) += df133 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX - 1) += df200 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX ) += df201 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX + 1) += df202 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY - 1, intXX + 2) += df203 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX - 1) += df210 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX ) += df211 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX + 1) += df212 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY , intXX + 2) += df213 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX - 1) += df220 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX ) += df221 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX + 1) += df222 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 1, intXX + 2) += df223 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX - 1) += df230 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX ) += df231 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX + 1) += df232 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 1, intYY + 2, intXX + 2) += df233 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX - 1) += df300 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX ) += df301 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX + 1) += df302 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY - 1, intXX + 2) += df303 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX - 1) += df310 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX ) += df311 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX + 1) += df312 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY , intXX + 2) += df313 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX - 1) += df320 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX ) += df321 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX + 1) += df322 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 1, intXX + 2) += df323 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX - 1) += df330 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX ) += df331 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX + 1) += df332 * my_val;
						DIRECT_A3D_ELEM(data, intZZ + 2, intYY + 2, intXX + 2) += df333 * my_val;

						// Store corresponding weights
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY - 1, intXX - 1) += df000 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY - 1, intXX ) += df001 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY - 1, intXX + 1) += df002 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY - 1, intXX + 2) += df003 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY , intXX - 1) += df010 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY , intXX ) += df011 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY , intXX + 1) += df012 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY , intXX + 2) += df013 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 1, intXX - 1) += df020 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 1, intXX ) += df021 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 1, intXX + 1) += df022 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 1, intXX + 2) += df023 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 2, intXX - 1) += df030 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 2, intXX ) += df031 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 2, intXX + 1) += df032 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ - 1, intYY + 2, intXX + 2) += df033 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY - 1, intXX - 1) += df100 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY - 1, intXX ) += df101 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY - 1, intXX + 1) += df102 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY - 1, intXX + 2) += df103 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY , intXX - 1) += df110 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY , intXX ) += df111 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY , intXX + 1) += df112 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY , intXX + 2) += df113 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY + 1, intXX - 1) += df120 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY + 1, intXX ) += df121 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY + 1, intXX + 1) += df122 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY + 1, intXX + 2) += df123 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY + 2, intXX - 1) += df130 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY + 2, intXX ) += df131 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY + 2, intXX + 1) += df132 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ , intYY + 2, intXX + 2) += df133 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY - 1, intXX - 1) += df200 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY - 1, intXX ) += df201 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY - 1, intXX + 1) += df202 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY - 1, intXX + 2) += df203 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY , intXX - 1) += df210 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY , intXX ) += df211 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY , intXX + 1) += df212 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY , intXX + 2) += df213 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 1, intXX - 1) += df220 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 1, intXX ) += df221 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 1, intXX + 1) += df222 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 1, intXX + 2) += df223 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 2, intXX - 1) += df230 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 2, intXX ) += df231 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 2, intXX + 1) += df232 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 1, intYY + 2, intXX + 2) += df233 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY - 1, intXX - 1) += df300 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY - 1, intXX ) += df301 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY - 1, intXX + 1) += df302 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY - 1, intXX + 2) += df303 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY , intXX - 1) += df310 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY , intXX ) += df311 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY , intXX + 1) += df312 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY , intXX + 2) += df313 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 1, intXX - 1) += df320 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 1, intXX ) += df321 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 1, intXX + 1) += df322 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 1, intXX + 2) += df323 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 2, intXX - 1) += df330 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 2, intXX ) += df331 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 2, intXX + 1) += df332 * my_weight;
						DIRECT_A3D_ELEM(weight, intZZ + 2, intYY + 2, intXX + 2) += df333 * my_weight;

					} // endif CUBIC
					else
					{
						REPORT_ERROR("BackProjector::backrotate3D%%ERROR: unrecognized interpolator ");
					}
				} // endif weight > 0.
			} // endif x-loop
		} // endif y-loop
	} // endif z-loop
}

void BackProjector::getLowResDataAndWeight(MultidimArray<Complex > &lowres_data, MultidimArray<DOUBLE> &lowres_weight,
		int lowres_r_max)
{

	int lowres_r2_max = padding_factor * padding_factor * lowres_r_max * lowres_r_max;
	int lowres_pad_size = 2 * (padding_factor * lowres_r_max + 1) + 1;

	// Check for dimension
	if (ref_dim != 3)
		REPORT_ERROR("BackProjector::getLowResDataAndWeight%%ERROR: only implemented for 3D case....");

	// Check lowres_r_max is not too big
	if (lowres_r_max > r_max)
		REPORT_ERROR("BackProjector::getLowResDataAndWeight%%ERROR: lowres_r_max is bigger than r_max");

	// Initialize lowres_data and low_res_weight arrays
	lowres_data.clear();
	lowres_data.resize(lowres_pad_size, lowres_pad_size, lowres_pad_size / 2 + 1);
	lowres_data.setXmippOrigin();
	lowres_data.xinit=0;
	lowres_weight.clear();
	lowres_weight.resize(lowres_pad_size, lowres_pad_size, lowres_pad_size / 2 + 1);
	lowres_weight.setXmippOrigin();
	lowres_weight.xinit=0;

	// fill lowres arrays with relevant values
	FOR_ALL_ELEMENTS_IN_ARRAY3D(lowres_data)
	{
		if (k*k + i*i + j*j <= lowres_r2_max)
		{
			A3D_ELEM(lowres_data, k, i, j) = A3D_ELEM(data, k , i, j);
			A3D_ELEM(lowres_weight, k, i, j) = A3D_ELEM(weight, k , i, j);
		}
	}

}

void BackProjector::setLowResDataAndWeight(MultidimArray<Complex > &lowres_data, MultidimArray<DOUBLE> &lowres_weight,
		int lowres_r_max)
{

	int lowres_r2_max = padding_factor * padding_factor * lowres_r_max * lowres_r_max;
	int lowres_pad_size = 2 * (padding_factor * lowres_r_max + 1) + 1;

	// Check for dimension
	if (ref_dim != 3)
		REPORT_ERROR("BackProjector::getLowResDataAndWeight%%ERROR: only implemented for 3D case....");

	// Check lowres_r_max is not too big
	if (lowres_r_max > r_max)
		REPORT_ERROR("BackProjector::getLowResDataAndWeight%%ERROR: lowres_r_max is bigger than r_max");

	// Check sizes of lowres_data and lowres_weight
	if (ZSIZE(lowres_data) != lowres_pad_size || YSIZE(lowres_data) != lowres_pad_size || XSIZE(lowres_data) != lowres_pad_size / 2 + 1)
		REPORT_ERROR("BackProjector::setLowResDataAndWeight%%ERROR: lowres_data is not of expected size...");
	if (ZSIZE(lowres_weight) != lowres_pad_size || YSIZE(lowres_weight) != lowres_pad_size || XSIZE(lowres_weight) != lowres_pad_size / 2 + 1)
		REPORT_ERROR("BackProjector::setLowResDataAndWeight%%ERROR: lowres_weight is not of expected size...");

	// Re-set origin to the expected place
	lowres_data.setXmippOrigin();
	lowres_data.xinit=0;
	lowres_weight.setXmippOrigin();
	lowres_weight.xinit=0;

	// Overwrite data and weight with the lowres arrays
	FOR_ALL_ELEMENTS_IN_ARRAY3D(lowres_data)
	{
		if (k*k + i*i + j*j <= lowres_r2_max)
		{
			A3D_ELEM(data, k, i, j) = A3D_ELEM(lowres_data, k , i, j);
			A3D_ELEM(weight, k, i, j) = A3D_ELEM(lowres_weight, k , i, j);
		}
	}

}


void BackProjector::getDownsampledAverage(MultidimArray<Complex > &avg)
{
	MultidimArray<DOUBLE> down_weight;

	// Pre-set down_data and down_weight sizes
	int down_size = 2 * (r_max + 1) + 1;
	int r2_max = r_max * r_max;
	// Short side of data array
	switch (ref_dim)
	{
	case 2:
	   avg.initZeros(down_size, down_size / 2 + 1);
	   break;
	case 3:
	   avg.initZeros(down_size, down_size, down_size / 2 + 1);
	   break;
	default:
	   REPORT_ERROR("BackProjector::getDownsampledAverage%%ERROR: Dimension of the data array should be 2 or 3");
	}
	// Set origin in the y.z-center, but on the left side for x.
	avg.setXmippOrigin();
	avg.xinit=0;
	// Resize down_weight the same as down_data
	down_weight.initZeros(avg);

	// Now calculate the down-sized sum
	int kp, ip, jp;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(data)
	{
		kp = ROUND((DOUBLE)k/padding_factor);
		ip = ROUND((DOUBLE)i/padding_factor);
		jp = ROUND((DOUBLE)j/padding_factor);

// TMP
//#define CHECK_SIZE
#ifdef CHECK_SIZE
		if (kp > FINISHINGZ(avg) || ip > FINISHINGY(avg) || jp > FINISHINGX(avg) ||
				kp < STARTINGZ(avg) || ip < STARTINGY(avg) || jp < STARTINGX(avg))
		{
			std::cerr << " kp= " << kp << " ip= " << ip << " jp= " << jp << std::endl;
			avg.printShape();
			REPORT_ERROR("BackProjector::getDownsampledAverage: indices out of range");
		}
#endif
		A3D_ELEM(avg, kp, ip, jp) += A3D_ELEM(data, k , i, j);
		A3D_ELEM(down_weight, kp, ip, jp) += A3D_ELEM(weight, k , i, j);
	}

	// Then enforce Hermitian symmetry in the downsampled arrays
	// We already took the average.... so not completely correct, but does not really matter for FSC calculation anyway
	// enforceHermitianSymmetry(avg, down_weight);

	// And enforce symmetry in the downsampled arrays
	symmetrise(avg, down_weight, r2_max);

	// Calculate the straightforward average in the downsampled arrays
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(avg)
	{
		if (DIRECT_MULTIDIM_ELEM(down_weight, n) > 0.)
			DIRECT_MULTIDIM_ELEM(avg, n) /= DIRECT_MULTIDIM_ELEM(down_weight, n);
		else
			DIRECT_MULTIDIM_ELEM(avg, n) = 0.;
	}


}

void BackProjector::calculateDownSampledFourierShellCorrelation(MultidimArray<Complex > &avg1,
																MultidimArray<Complex > &avg2,
																MultidimArray<DOUBLE> &fsc)
{

    if (!avg1.sameShape(avg2))
    	REPORT_ERROR("ERROR BackProjector::calculateDownSampledFourierShellCorrelation: two arrays have different sizes");

    MultidimArray<DOUBLE> num, den1, den2;
    num.initZeros(ori_size/2 + 1);
    den1.initZeros(num);
    den2.initZeros(num);
    fsc.initZeros(num);

    FOR_ALL_ELEMENTS_IN_ARRAY3D(avg1)
    {
    	DOUBLE R = sqrt(k*k + i*i + j*j);
        if (R > r_max)
            continue;
        int idx=ROUND(R);
        Complex z1=A3D_ELEM(avg1, k, i, j);
        Complex z2=A3D_ELEM(avg2, k, i, j);
        DOUBLE absz1=abs(z1);
        DOUBLE absz2=abs(z2);
        num(idx)+=(conj(z1) * z2).real;
        den1(idx)+= absz1*absz1;
        den2(idx)+= absz2*absz2;
    }

    FOR_ALL_ELEMENTS_IN_ARRAY1D(fsc)
    {
    	if (den1(i)*den2(i) > 0.)
    		fsc(i) = num(i)/sqrt(den1(i)*den2(i));
    }

    // Always set zero-resolution shell to FSC=1
    // Raimond Ravelli reported a problem with FSC=1 at res=0 on 13feb2013...
    // (because of a suboptimal normalisation scheme, but anyway)
    fsc(0) = 1.;

}


void BackProjector::reconstruct(MultidimArray<DOUBLE> &vol_out,
                                int max_iter_preweight,
                                bool do_map,
                                DOUBLE tau2_fudge,
                                MultidimArray<DOUBLE> &tau2,
                                MultidimArray<DOUBLE> &sigma2,
                                MultidimArray<DOUBLE> &data_vs_prior,
                                MultidimArray<DOUBLE> fsc, // only input
                                DOUBLE normalise,
                                bool update_tau2_with_fsc,
                                bool is_whole_instead_of_half,
                                int nr_threads,
                                int minres_map)

{


    FourierTransformer transformer;
	// The threads are giving me a headache. Let's switch them off
    // Somehow I get lots of bad/non-reproducible errors when having these...
    //transformer.setThreadsNumber(nr_threads);
	MultidimArray<Complex > Fconv;
	MultidimArray<DOUBLE> Fweight;
        // Fnewweight can become too large for a float: always keep this one in double-precision
        MultidimArray<double> Fnewweight;

	int max_r2 = r_max * r_max * padding_factor * padding_factor;

//#define DEBUG_RECONSTRUCT
#ifdef DEBUG_RECONSTRUCT
	Image<DOUBLE> ttt;
	FileName fnttt;
	ttt()=weight;
	ttt.write("reconstruct_initial_weight.spi");
#endif

	// At the x=0 line, we have collected either the positive y-z coordinate, or its negative Friedel pair.
	// Sum these two together for both the data and the weight arrays
	enforceHermitianSymmetry(data, weight);

#ifdef DEBUG_RECONSTRUCT
	ttt()=weight;
	ttt.write("reconstruct_hermitian_weight.spi");
#endif

	// First enforce Hermitian symmetry, then symmetry!
	// This way the redundancy at the x=0 plane is handled correctly
	symmetrise(data, weight, max_r2);
#ifdef DEBUG_RECONSTRUCT
	ttt()=weight;
	ttt.write("reconstruct_symmetrised_weight.spi");
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(data)
	{
		DIRECT_MULTIDIM_ELEM(ttt(), n) = DIRECT_MULTIDIM_ELEM(data, n).real;
	}
	ttt.write("reconstruct_symmetrised_data_real.spi");

	std::cerr << " pad_size= " << pad_size << " padding_factor= " << padding_factor << " max_r2= " << max_r2 << std::endl;
#endif


	// Set Fweight, Fnewweight and Fconv to the right size
	if (ref_dim == 2)
		vol_out.resize(pad_size, pad_size);
	else
		vol_out.resize(pad_size, pad_size, pad_size);

	transformer.setReal(vol_out);
	transformer.getFourierAlias(Fconv);

	// clear vol_out to save memory!
	vol_out.clear();

	Fweight.resize(Fconv);
	Fnewweight.resize(Fconv);
	// Go from projector-centered to FFTW-uncentered
	decenter(weight, Fweight, max_r2);

	// Take oversampling into account
	DOUBLE oversampling_correction = (ref_dim == 3) ? (padding_factor * padding_factor * padding_factor) : (padding_factor * padding_factor);
	MultidimArray<DOUBLE> counter;


	// First calculate the radial average of the (inverse of the) power of the noise in the reconstruction
	// This is the left-hand side term in the nominator of the Wiener-filter-like update formula
	// and it is stored inside the weight vector
	// Then, if (do_map) add the inverse of tau2-spectrum values to the weight
	sigma2.initZeros(ori_size/2 + 1);
	counter.initZeros(ori_size/2 + 1);
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv)
	{
		int r2 = kp * kp + ip * ip + jp * jp;
		if (r2 < max_r2)
		{

#define DEBUG_ZSY
#ifdef DEBUG_ZSY
	std::cerr << " kp= " << kp << " ip= " << ip << " jp= " << jp << std::endl;
	std::cerr << " padding_factor= " << padding_factor << " Fweight= " << DIRECT_A3D_ELEM(Fweight, k, i, j) << std::endl;
#endif

			int ires = ROUND( sqrt((DOUBLE)r2) / padding_factor );
			DOUBLE invw = oversampling_correction * DIRECT_A3D_ELEM(Fweight, k, i, j);
			DIRECT_A1D_ELEM(sigma2, ires) += invw;
			DIRECT_A1D_ELEM(counter, ires) += 1.;
		}
	}

#define DEBUG_ZSY
#define DEBUG_ZSY_WEIGHT     // for the weight which is negative 
#ifdef DEBUG_ZSY
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(sigma2)
	{
	std::cerr << " i= " << i << " sigma2= " << DIRECT_A1D_ELEM(sigma2, i) << std::endl;
	}
#endif	

	// Average (inverse of) sigma2 in reconstruction
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(sigma2)
	{
		if (DIRECT_A1D_ELEM(sigma2, i) > 1e-10)
			DIRECT_A1D_ELEM(sigma2, i) = DIRECT_A1D_ELEM(counter, i) / DIRECT_A1D_ELEM(sigma2, i);
		else if (DIRECT_A1D_ELEM(sigma2, i) == 0)
			DIRECT_A1D_ELEM(sigma2, i) = 0.;

#ifdef DEBUG_ZSY_WEIGHT
		else if (DIRECT_A1D_ELEM(sigma2, i) < 0)
			DIRECT_A1D_ELEM(sigma2, i) = - DIRECT_A1D_ELEM(sigma2, i);
#endif

		else
		{
			std::cerr << " DIRECT_A1D_ELEM(sigma2, i)= " << DIRECT_A1D_ELEM(sigma2, i) << std::endl;
			REPORT_ERROR("BackProjector::reconstruct: ERROR: unexpectedly small, yet non-zero sigma2 value, this should not happen...a");
		}
	}

	if (update_tau2_with_fsc)
	{
		tau2.resize(ori_size/2 + 1);
		data_vs_prior.initZeros(ori_size/2 + 1);
		// Then calculate new tau2 values, based on the FSC
		if (!fsc.sameShape(sigma2) || !fsc.sameShape(tau2))
		{
			fsc.printShape(std::cerr);
			tau2.printShape(std::cerr);
			sigma2.printShape(std::cerr);
			REPORT_ERROR("ERROR BackProjector::reconstruct: sigma2, tau2 and fsc have different sizes");
		}
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(sigma2)
		{
			// FSC cannot be negative or zero for conversion into tau2
			DOUBLE myfsc = XMIPP_MAX(0.001, DIRECT_A1D_ELEM(fsc, i));
			if (is_whole_instead_of_half)
			{
				// Factor two because of twice as many particles
				// Sqrt-term to get 60-degree phase errors....
				myfsc = sqrt(2. * myfsc / (myfsc + 1.));
			}
			myfsc = XMIPP_MIN(0.999, myfsc);
			DOUBLE myssnr = myfsc / (1. - myfsc);
			DOUBLE fsc_based_tau = myssnr * DIRECT_A1D_ELEM(sigma2, i);
			DIRECT_A1D_ELEM(tau2, i) = fsc_based_tau;
			// data_vs_prior is merely for reporting: it is not used for anything in the reconstruction
			DIRECT_A1D_ELEM(data_vs_prior, i) = myssnr;

		}
	}

	// Apply MAP-additional term to the Fnewweight array
	// This will regularise the actual reconstruction
	if (do_map)
	{
		// Then, add the inverse of tau2-spectrum values to the weight
		// and also calculate spherical average of data_vs_prior ratios
		if (!update_tau2_with_fsc)
			data_vs_prior.initZeros(ori_size/2 + 1);
		counter.initZeros(ori_size/2 + 1);
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv)
 		{
			int r2 = kp * kp + ip * ip + jp * jp;
			if (r2 < max_r2)
			{
				int ires = ROUND( sqrt((DOUBLE)r2) / padding_factor );
				DOUBLE invw = DIRECT_A3D_ELEM(Fweight, k, i, j);

				DOUBLE invtau2;
				if (DIRECT_A1D_ELEM(tau2, ires) > 0.)
				{
					// Calculate inverse of tau2
					invtau2 = 1. / (oversampling_correction * tau2_fudge * DIRECT_A1D_ELEM(tau2, ires));
				}
				else if (DIRECT_A1D_ELEM(tau2, ires) == 0.)
				{
					// If tau2 is zero, use small value instead
					invtau2 = 1./ ( 0.001 * invw);
				}
				else
				{
					std::cerr << " sigma2= " << sigma2 << std::endl;
					std::cerr << " fsc= " << fsc << std::endl;
					std::cerr << " tau2= " << tau2 << std::endl;
					REPORT_ERROR("ERROR BackProjector::reconstruct: Negative or zero values encountered for tau2 spectrum!");
				}

				// Keep track of spectral evidence-to-prior ratio and remaining noise in the reconstruction
				if (!update_tau2_with_fsc)
					DIRECT_A1D_ELEM(data_vs_prior, ires) += invw / invtau2;
				DIRECT_A1D_ELEM(counter, ires) += 1.;

				// Only for (ires >= minres_map) add Wiener-filter like term
				if (ires >= minres_map)
				{
					// Now add the inverse-of-tau2_class term
					invw += invtau2;
					// Store the new weight again in Fweight
					DIRECT_A3D_ELEM(Fweight, k, i, j) = invw;
				}
			}
		}

		// Average data_vs_prior
		if (!update_tau2_with_fsc)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(data_vs_prior)
			{
				if (i > r_max)
					DIRECT_A1D_ELEM(data_vs_prior, i) = 0.;
				else if (DIRECT_A1D_ELEM(counter, i) < 0.001)
					DIRECT_A1D_ELEM(data_vs_prior, i) = 999.;
				else
					DIRECT_A1D_ELEM(data_vs_prior, i) /= DIRECT_A1D_ELEM(counter, i);
			}
		}

	} //end if do_map

	// Divide both data and Fweight by normalisation factor to prevent FFT's with very large values....
#ifdef DEBUG_RECONSTRUCT
	std::cerr << " normalise= " << normalise << std::endl;
#endif
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fweight)
	{
		DIRECT_MULTIDIM_ELEM(Fweight, n) /= normalise;
	}
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(data)
	{
		DIRECT_MULTIDIM_ELEM(data, n) /= normalise;
	}

	// Initialise Fnewweight with 1's and 0's. (also see comments below)
	FOR_ALL_ELEMENTS_IN_ARRAY3D(weight)
	{
		if (k * k + i * i + j * j < max_r2)
			A3D_ELEM(weight, k, i, j) = 1.;
		else
			A3D_ELEM(weight, k, i, j) = 0.;
	}
	decenter(weight, Fnewweight, max_r2);

	// Iterative algorithm as in  Eq. [14] in Pipe & Menon (1999)
	// or Eq. (4) in Matej (2001)
	for (int iter = 0; iter < max_iter_preweight; iter++)
	{

		// Set Fnewweight * Fweight in the transformer
		// In Matej et al (2001), weights w_P^i are convoluted with the kernel,
		// and the initial w_P^0 are 1 at each sampling point
		// Here the initial weights are also 1 (see initialisation Fnewweight above),
		// but each "sampling point" counts "Fweight" times!
		// That is why Fnewweight is multiplied by Fweight prior to the convolution
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv)
		{
			DIRECT_MULTIDIM_ELEM(Fconv, n) = DIRECT_MULTIDIM_ELEM(Fnewweight, n) * DIRECT_MULTIDIM_ELEM(Fweight, n);
		}

                // convolute through Fourier-transform (as both grids are rectangular)
                // Note that convoluteRealSpace acts on the complex array inside the transformer
                convoluteBlobRealSpace(transformer);

                DOUBLE w, corr_min = 99.e99, corr_max = -99.e99, corr_avg=0., corr_nn=0.;
                FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv)
                {
                    if (kp * kp + ip * ip + jp * jp < max_r2)
                    {

        		// Make sure no division by zero can occur....
        		w = XMIPP_MAX(1e-6, abs(DIRECT_A3D_ELEM(Fconv, k, i, j)));
        		// Monitor min, max and avg conv_weight
        		corr_min = XMIPP_MIN(corr_min, w);
        		corr_max = XMIPP_MAX(corr_max, w);
        		corr_avg += w;
        		corr_nn += 1.;
        		// Apply division of Eq. [14] in Pipe & Menon (1999)
        		DIRECT_A3D_ELEM(Fnewweight, k, i, j) /= w;
                    }
                }

#ifdef DEBUG_RECONSTRUCT
        std::cerr << " PREWEIGHTING ITERATION: "<< iter + 1 << " OF " << max_iter_preweight << std::endl;
        // report of maximum and minimum values of current conv_weight
        std::cerr << " corr_avg= " << corr_avg / corr_nn << std::endl;
        std::cerr << " corr_min= " << corr_min << std::endl;
        std::cerr << " corr_max= " << corr_max << std::endl;
#endif
	}

#ifdef DEBUG_RECONSTRUCT
	Image<double> tttt;
	tttt()=Fnewweight;
	tttt.write("reconstruct_gridding_weight.spi");
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv)
	{
		DIRECT_MULTIDIM_ELEM(ttt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fconv, n));
	}
	ttt.write("reconstruct_gridding_correction_term.spi");
#endif

	// Clear memory
	Fweight.clear();

	// Note that Fnewweight now holds the approximation of the inverse of the weights on a regular grid

	// Now do the actual reconstruction with the data array
	// Apply the iteratively determined weight
	Fconv.initZeros(); // to remove any stuff from the input volume
	decenter(data, Fconv, max_r2);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv)
	{
#ifdef  FLOAT_PRECISION
            // Prevent numerical instabilities in single-precision reconstruction with very unevenly sampled orientations
            if (DIRECT_MULTIDIM_ELEM(Fnewweight, n) > 1e20)
                DIRECT_MULTIDIM_ELEM(Fnewweight, n) = 1e20;
#endif
		DIRECT_MULTIDIM_ELEM(Fconv, n) *= DIRECT_MULTIDIM_ELEM(Fnewweight, n);
	}

	// Clear memory
	Fnewweight.clear();

// Gridding theory says one now has to interpolate the fine grid onto the coarse one using a blob kernel
// and then do the inverse transform and divide by the FT of the blob (i.e. do the gridding correction)
// In practice, this gives all types of artefacts (perhaps I never found the right implementation?!)
// Therefore, window the Fourier transform and then do the inverse transform
//#define RECONSTRUCT_CONVOLUTE_BLOB
#ifdef RECONSTRUCT_CONVOLUTE_BLOB

	// Apply the same blob-convolution as above to the data array
	// Mask real-space map beyond its original size to prevent aliasing in the downsampling step below
	convoluteBlobRealSpace(transformer, true);

	// Now just pick every 3rd pixel in Fourier-space (i.e. down-sample)
	// and do a final inverse FT
	if (ref_dim == 2)
		vol_out.resize(ori_size, ori_size);
	else
		vol_out.resize(ori_size, ori_size, ori_size);

	FourierTransformer transformer2;
	MultidimArray<Complex > Ftmp;
	transformer2.setReal(vol_out); // cannot use the first transformer because Fconv is inside there!!
	transformer2.getFourierAlias(Ftmp);
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Ftmp)
	{
		if (kp * kp + ip * ip + jp * jp < r_max * r_max)
		{
			DIRECT_A3D_ELEM(Ftmp, k, i, j) = FFTW_ELEM(Fconv, kp * padding_factor, ip * padding_factor, jp * padding_factor);
		}
		else
		{
			DIRECT_A3D_ELEM(Ftmp, k, i, j) = 0.;
		}
	}

	// inverse FFT leaves result in vol_out
	transformer2.inverseFourierTransform();

	// Shift the map back to its origin
	CenterFFT(vol_out, false);

	// Un-normalize FFTW (because original FFTs were done with the size of 2D FFTs)
	if (ref_dim==3)
		vol_out /= ori_size;

	// Mask out corners to prevent aliasing artefacts
	softMaskOutsideMap(vol_out);

	// Gridding correction for the blob
	DOUBLE normftblob = tab_ftblob(0.);
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_out)
	{

		DOUBLE r = sqrt((DOUBLE)(k*k+i*i+j*j));
		DOUBLE rval = r / (ori_size * padding_factor);
		A3D_ELEM(vol_out, k, i, j) /= tab_ftblob(rval) / normftblob;
		//if (k==0 && i==0)
		//	std::cerr << " j= " << j << " rval= " << rval << " tab_ftblob(rval) / normftblob= " << tab_ftblob(rval) / normftblob << std::endl;
	}


#else



	// rather than doing the blob-convolution to downsample the data array, do a windowing operation:
	// This is the same as convolution with a SINC. It seems to give better maps.
	// Then just make the blob look as much as a SINC as possible....
	// The "standard" r1.9, m2 and a15 blob looks quite like a sinc until the first zero (perhaps that's why it is standard?)
	//for (DOUBLE r = 0.1; r < 10.; r+=0.01)
	//{
	//	DOUBLE sinc = sin(PI * r / padding_factor ) / ( PI * r / padding_factor);
	//	std::cout << " r= " << r << " sinc= " << sinc << " blob= " << blob_val(r, blob) << std::endl;
	//}

	// Now do inverse FFT and window to original size in real-space
	// Pass the transformer to prevent making and clearing a new one before clearing the one declared above....
	// The latter may give memory problems as detected by electric fence....
	windowToOridimRealSpace(transformer, Fconv, vol_out, nr_threads);

#endif

#ifdef DEBUG_RECONSTRUCT
	ttt()=vol_out;
	ttt.write("reconstruct_before_gridding_correction.spi");
#endif

	// Correct for the linear/nearest-neighbour interpolation that led to the data array
	griddingCorrect(vol_out);


	// If the tau-values were calculated based on the FSC, then now re-calculate the power spectrum of the actual reconstruction
	if (update_tau2_with_fsc)
	{

		// New tau2 will be the power spectrum of the new map
		MultidimArray<DOUBLE> spectrum, count;

		// Calculate this map's power spectrum
		// Don't call getSpectrum() because we want to use the same transformer object to prevent memory trouble....
		spectrum.initZeros(XSIZE(vol_out));
	    count.initZeros(XSIZE(vol_out));
	    // recycle the same transformer for all images
	    transformer.FourierTransform(vol_out, Fconv, false);
	    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv)
	    {
	    	long int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
	    	spectrum(idx) += norm(dAkij(Fconv, k, i, j));
	        count(idx) += 1.;
	    }
	    spectrum /= count;

		// Factor two because of two-dimensionality of the complex plane
		// (just like sigma2_noise estimates, the power spectra should be divided by 2)
		DOUBLE normfft = (ref_dim == 3 && data_dim == 2) ? (DOUBLE)(ori_size * ori_size) : 1.;
		spectrum *= normfft / 2.;

		// New SNR^MAP will be power spectrum divided by the noise in the reconstruction (i.e. sigma2)
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(data_vs_prior)
		{
			DIRECT_MULTIDIM_ELEM(tau2, n) =  tau2_fudge * DIRECT_MULTIDIM_ELEM(spectrum, n);
		}

	}

	// Completely empty the transformer object
	transformer.cleanup();

#ifdef DEBUG_RECONSTRUCT
    std::cerr<<"done with reconstruct"<<std::endl;
#endif

}

void BackProjector::enforceHermitianSymmetry(MultidimArray<Complex > &my_data,
											 MultidimArray<DOUBLE> &my_weight)
{

	for (int iz = STARTINGZ(my_data); iz <=FINISHINGZ(my_data); iz++)
	{
		// Make sure all points are only included once.
		int starty = (iz < 0) ? 0 : 1;
		for (int iy = starty; iy <= FINISHINGY(my_data); iy++)
		{
			// I just need to sum the two points, not divide by 2!
			Complex fsum = (A3D_ELEM(my_data, iz, iy, 0) + conj(A3D_ELEM(my_data, -iz, -iy, 0)));
			A3D_ELEM(my_data, iz, iy, 0) = fsum;
			A3D_ELEM(my_data, -iz, -iy, 0) = conj(fsum);
			DOUBLE sum = (A3D_ELEM(my_weight, iz, iy, 0) + A3D_ELEM(my_weight, -iz, -iy, 0));
			A3D_ELEM(my_weight, iz, iy, 0) = sum;
			A3D_ELEM(my_weight, -iz, -iy, 0) = sum;
		}
	}

}

void BackProjector::symmetrise(MultidimArray<Complex > &my_data,
		 MultidimArray<DOUBLE> &my_weight, int my_rmax2)
{

//#define DEBUG_SYMM
#ifdef DEBUG_SYMM
	std::cerr << " SL.SymsNo()= " << SL.SymsNo() << std::endl;
	std::cerr << " SL.true_symNo= " << SL.true_symNo << std::endl;
#endif

	if (SL.SymsNo() > 0 && ref_dim == 3)
	{
		Matrix2D<DOUBLE> L(4, 4), R(4, 4); // A matrix from the list
		MultidimArray<DOUBLE> sum_weight;
		MultidimArray<Complex > sum_data;
        DOUBLE x, y, z, fx, fy, fz, xp, yp, zp, r2;
        bool is_neg_x;
        int x0, x1, y0, y1, z0, z1;
    	Complex d000, d001, d010, d011, d100, d101, d110, d111;
    	Complex dx00, dx01, dx10, dx11, dxy0, dxy1;
    	DOUBLE dd000, dd001, dd010, dd011, dd100, dd101, dd110, dd111;
    	DOUBLE ddx00, ddx01, ddx10, ddx11, ddxy0, ddxy1;

        // First symmetry operator (not stored in SL) is the identity matrix
		sum_weight = my_weight;
		sum_data = my_data;
		// Loop over all other symmetry operators
	    for (int isym = 0; isym < SL.SymsNo(); isym++)
	    {
	        SL.get_matrices(isym, L, R);
#ifdef DEBUG_SYMM
	        std::cerr << " isym= " << isym << " R= " << R << std::endl;
#endif

	        // Loop over all points in the output (i.e. rotated, or summed) array
	        FOR_ALL_ELEMENTS_IN_ARRAY3D(sum_weight)
	        {

	        	x = (DOUBLE)j; // STARTINGX(sum_weight) is zero!
	        	y = (DOUBLE)i;
	        	z = (DOUBLE)k;
	        	r2 = x*x + y*y + z*z;
	        	if (r2 <= my_rmax2)
	        	{
	        		// coords_output(x,y) = A * coords_input (xp,yp)
					xp = x * R(0, 0) + y * R(0, 1) + z * R(0, 2);
					yp = x * R(1, 0) + y * R(1, 1) + z * R(1, 2);
					zp = x * R(2, 0) + y * R(2, 1) + z * R(2, 2);

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
					y0 -=  STARTINGY(my_data);
					y1 = y0 + 1;

					z0 = FLOOR(zp);
					fz = zp - z0;
					z0 -= STARTINGZ(my_data);
					z1 = z0 + 1;

#ifdef CHECK_SIZE
					if (x0 < 0 || y0 < 0 || z0 < 0 ||
						x1 < 0 || y1 < 0 || z1 < 0 ||
						x0 >= XSIZE(my_data) || y0  >= YSIZE(my_data) || z0 >= ZSIZE(my_data) ||
						x1 >= XSIZE(my_data) || y1  >= YSIZE(my_data)  || z1 >= ZSIZE(my_data) 	)
					{
						std::cerr << " x0= " << x0 << " y0= " << y0 << " z0= " << z0 << std::endl;
						std::cerr << " x1= " << x1 << " y1= " << y1 << " z1= " << z1 << std::endl;
						my_data.printShape();
						REPORT_ERROR("BackProjector::symmetrise: checksize!!!");
					}
#endif
					// First interpolate (complex) data
					d000 = DIRECT_A3D_ELEM(my_data, z0, y0, x0);
					d001 = DIRECT_A3D_ELEM(my_data, z0, y0, x1);
					d010 = DIRECT_A3D_ELEM(my_data, z0, y1, x0);
					d011 = DIRECT_A3D_ELEM(my_data, z0, y1, x1);
					d100 = DIRECT_A3D_ELEM(my_data, z1, y0, x0);
					d101 = DIRECT_A3D_ELEM(my_data, z1, y0, x1);
					d110 = DIRECT_A3D_ELEM(my_data, z1, y1, x0);
					d111 = DIRECT_A3D_ELEM(my_data, z1, y1, x1);

					dx00 = LIN_INTERP(fx, d000, d001);
					dx01 = LIN_INTERP(fx, d100, d101);
					dx10 = LIN_INTERP(fx, d010, d011);
					dx11 = LIN_INTERP(fx, d110, d111);
					dxy0 = LIN_INTERP(fy, dx00, dx10);
					dxy1 = LIN_INTERP(fy, dx01, dx11);

					// Take complex conjugated for half with negative x
					if (is_neg_x)
						A3D_ELEM(sum_data, k, i, j) += conj(LIN_INTERP(fz, dxy0, dxy1));
					else
						A3D_ELEM(sum_data, k, i, j) += LIN_INTERP(fz, dxy0, dxy1);

					// Then interpolate (real) weight
					dd000 = DIRECT_A3D_ELEM(my_weight, z0, y0, x0);
					dd001 = DIRECT_A3D_ELEM(my_weight, z0, y0, x1);
					dd010 = DIRECT_A3D_ELEM(my_weight, z0, y1, x0);
					dd011 = DIRECT_A3D_ELEM(my_weight, z0, y1, x1);
					dd100 = DIRECT_A3D_ELEM(my_weight, z1, y0, x0);
					dd101 = DIRECT_A3D_ELEM(my_weight, z1, y0, x1);
					dd110 = DIRECT_A3D_ELEM(my_weight, z1, y1, x0);
					dd111 = DIRECT_A3D_ELEM(my_weight, z1, y1, x1);

					ddx00 = LIN_INTERP(fx, dd000, dd001);
					ddx01 = LIN_INTERP(fx, dd100, dd101);
					ddx10 = LIN_INTERP(fx, dd010, dd011);
					ddx11 = LIN_INTERP(fx, dd110, dd111);
					ddxy0 = LIN_INTERP(fy, ddx00, ddx10);
					ddxy1 = LIN_INTERP(fy, ddx01, ddx11);

					A3D_ELEM(sum_weight, k, i, j) +=  LIN_INTERP(fz, ddxy0, ddxy1);

	        	} // end if r2 <= my_rmax2

	        } // end loop over all elements of sum_weight

	    } // end loop over symmetry operators

	    my_data = sum_data;
	    my_weight = sum_weight;
	    // Average
	    // The division should only be done if we would search all (C1) directions, not if we restrict the angular search!
	    /*
	    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(data)
	    {
	    	DIRECT_MULTIDIM_ELEM(data, n) = DIRECT_MULTIDIM_ELEM(sum_data, n) / (DOUBLE)(SL.SymsNo() + 1);
	    	DIRECT_MULTIDIM_ELEM(weight, n) = DIRECT_MULTIDIM_ELEM(sum_weight, n) / (DOUBLE)(SL.SymsNo() + 1);
	    }
	    */
	}

}

void BackProjector::convoluteBlobRealSpace(FourierTransformer &transformer, bool do_mask)
{

	MultidimArray<DOUBLE> Mconv;
	int padhdim = pad_size / 2;

	// Set up right dimension of real-space array
	// TODO: resize this according to r_max!!!
	if (ref_dim==2)
		Mconv.resize(pad_size, pad_size);
	else
		Mconv.resize(pad_size, pad_size, pad_size);

	// inverse FFT
	transformer.setReal(Mconv);
	transformer.inverseFourierTransform();

	// Blob normalisation in Fourier space
	DOUBLE normftblob = tab_ftblob(0.);

	// TMP DEBUGGING
	//struct blobtype blob;
	//blob.order = 0;
	//blob.radius = 1.9 * padding_factor;
	//blob.alpha = 15;

	// Multiply with FT of the blob kernel
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Mconv)
    {
		int kp = (k < padhdim) ? k : k - pad_size;
		int ip = (i < padhdim) ? i : i - pad_size;
		int jp = (j < padhdim) ? j : j - pad_size;
    	DOUBLE rval = sqrt ( (DOUBLE)(kp * kp + ip * ip + jp * jp) ) / (ori_size * padding_factor);
    	//if (kp==0 && ip==0 && jp > 0)
		//	std::cerr << " jp= " << jp << " rval= " << rval << " tab_ftblob(rval) / normftblob= " << tab_ftblob(rval) / normftblob << " ori_size/2= " << ori_size/2 << std::endl;
    	// In the final reconstruction: mask the real-space map beyond its original size to prevent aliasing ghosts
    	// Note that rval goes until 1/2 in the oversampled map
    	if (do_mask && rval > 1./(2. * padding_factor))
    		DIRECT_A3D_ELEM(Mconv, k, i, j) = 0.;
    	else
    		DIRECT_A3D_ELEM(Mconv, k, i, j) *= (tab_ftblob(rval) / normftblob);
    }

    // forward FFT to go back to Fourier-space
    transformer.FourierTransform();

}

void BackProjector::windowToOridimRealSpace(FourierTransformer &transformer, MultidimArray<Complex > &Fin, MultidimArray<DOUBLE> &Mout, int nr_threads)
{

	MultidimArray<Complex > Ftmp;
	int padoridim = padding_factor * ori_size;
	DOUBLE normfft;

//#define DEBUG_WINDOWORIDIMREALSPACE
#ifdef DEBUG_WINDOWORIDIMREALSPACE
	Image<DOUBLE> tt;
	tt().resize(ZSIZE(Fin), YSIZE(Fin), XSIZE(Fin));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fin)
	{
		DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fin, n));
	}
	tt.write("windoworidim_Fin.spi");
#endif

	if (ref_dim == 2)
	{
		Mout.resize(padoridim, padoridim);
		normfft = (DOUBLE)(padding_factor * padding_factor);
	}
	else
	{
		Mout.resize(padoridim, padoridim, padoridim);
		if (data_dim == 3)
			normfft = (DOUBLE)(padding_factor * padding_factor * padding_factor);
		else
			normfft = (DOUBLE)(padding_factor * padding_factor * padding_factor * ori_size);
	}
	Mout.setXmippOrigin();

	// Resize incoming complex array to the correct size
	windowFourierTransform(Fin, Ftmp, padoridim);

#ifdef DEBUG_WINDOWORIDIMREALSPACE
	tt().resize(ZSIZE(Ftmp), YSIZE(Ftmp), XSIZE(Ftmp));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Ftmp)
	{
		DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Ftmp, n));
	}
	tt.write("windoworidim_Fresized.spi");
#endif

	// Do the inverse FFT
	transformer.inverseFourierTransform(Ftmp, Mout);
	Mout.setXmippOrigin();

	// Shift the map back to its origin
	CenterFFT(Mout,true);

#ifdef DEBUG_WINDOWORIDIMREALSPACE
	tt()=Mout;
	tt.write("windoworidim_Munwindowed.spi");
#endif

	// Window in real-space
	if (ref_dim==2)
	{
		Mout.window(FIRST_XMIPP_INDEX(ori_size), FIRST_XMIPP_INDEX(ori_size),
				       LAST_XMIPP_INDEX(ori_size), LAST_XMIPP_INDEX(ori_size));
	}
	else
	{
		Mout.window(FIRST_XMIPP_INDEX(ori_size), FIRST_XMIPP_INDEX(ori_size), FIRST_XMIPP_INDEX(ori_size),
				       LAST_XMIPP_INDEX(ori_size), LAST_XMIPP_INDEX(ori_size), LAST_XMIPP_INDEX(ori_size));
	}
	Mout.setXmippOrigin();

	// Normalisation factor of FFTW
	// The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
	Mout /= normfft;

#ifdef DEBUG_WINDOWORIDIMREALSPACE
	tt()=Mout;
	tt.write("windoworidim_Mwindowed.spi");
#endif

	// Mask out corners to prevent aliasing artefacts
	softMaskOutsideMap(Mout);

#ifdef DEBUG_WINDOWORIDIMREALSPACE
	tt()=Mout;
	tt.write("windoworidim_Mwindowed_masked.spi");
	FourierTransformer ttf;
	ttf.FourierTransform(Mout, Ftmp);
	tt().resize(ZSIZE(Ftmp), YSIZE(Ftmp), XSIZE(Ftmp));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Ftmp)
	{
		DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Ftmp, n));
	}
	tt.write("windoworidim_Fnew.spi");
#endif


}
