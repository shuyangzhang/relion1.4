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

#ifndef PARTICLE_POLISHER_MPI_H_
#define PARTICLE_POLISHER_MPI_H_

#include "src/mpi.h"
#include "src/parallel.h"
#include "src/particle_polisher.h"

class ParticlePolisherMpi : public ParticlePolisher
{
private:
	MpiNode *node;

public:
	/** Destructor, calls MPI_Finalize */
    ~ParticlePolisherMpi()
    {
        delete node;
    }

    /** Read
     * This could take care of mpi-parallelisation-dependent variables
     */
    void read(int argc, char **argv);

	// Parallelized fit the beam-induced translations for all average micrographs
	void fitMovementsAllMicrographs();

	// Parallelized calculation of B-factors for single-frame reconstructions
	void calculateAllSingleFrameReconstructionsAndBfactors();

	// Parallelized movie frame re-alignment for all micrographs
	void polishParticlesAllMicrographs();

	// Calculate two half-reconstructions from shiny particles and calculate FSC-weighted average map, store that in refvol
	void reconstructShinyParticlesAndFscWeight(int ipass);

	// Parallelized optimisation of beamtilt for all micrographs
	void optimiseBeamTilt();

	// Parallelized run function
    void run();

};



#endif /* PARTICLE_POLISHER_MPI_H_ */
