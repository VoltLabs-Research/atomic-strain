#pragma once

#include <memory>
#include <vector>
#include <atomic>

#include <volt/core/volt.h>
#include <volt/core/simulation_cell.h>
#include <volt/core/particle_property.h>

namespace Volt{

class CutoffNeighborFinder;
class AtomicStrainLog;

class AtomicStrainModifier{
public:
	// The strain over the whole frame, for the summary table and the log.
	struct StrainSummary{
		int count = 0;
		double totalShear = 0.0;
		double maxShear = 0.0;
		double totalVolumetric = 0.0;

		double meanShear() const{ return count > 0 ? totalShear / count : 0.0; }
		double meanVolumetric() const{ return count > 0 ? totalVolumetric / count : 0.0; }
	};

	class AtomicStrainEngine{
	public:
		AtomicStrainEngine(
			Particles::ParticleProperty* positions,
			const SimulationCell& cell,
			Particles::ParticleProperty* refPositions,
			const SimulationCell& refCell,
			Particles::ParticleProperty* identifiers,
			Particles::ParticleProperty* refIdentifiers,
			double cutoff,
			bool eliminateCellDeformation,
			bool assumeUnwrappedCoordinates,
			bool calculateDeformationGradients,
			bool calculateStrainTensors,
			bool calculateNonaffineSquaredDisplacements,
			bool calculatePolarDecomposition = false
		);

		void perform(AtomicStrainLog& log);

		std::shared_ptr<Particles::ParticleProperty> rotationTensors() const{
			return _rotationTensors;
		}

		std::shared_ptr<Particles::ParticleProperty> stretchTensors() const{
			return _stretchTensors;
		}

		std::shared_ptr<Particles::ParticleProperty> shearStrains() const{
			return _shearStrains;
		}

		std::shared_ptr<Particles::ParticleProperty> volumetricStrains() const{
			return _volumetricStrains;
		}

		std::shared_ptr<Particles::ParticleProperty> nonaffineSquaredDisplacements() const{
			return _nonaffineSquaredDisplacements;
		}

		std::shared_ptr<Particles::ParticleProperty> invalidParticles() const{
			return _invalidParticles;
		}

		std::shared_ptr<Particles::ParticleProperty> strainTensors() const{
			return _strainTensors;
		}

		std::shared_ptr<Particles::ParticleProperty> deformationGradients() const{
			return _deformationGradients;
		}

		std::size_t numInvalidParticles() const{
			return _numInvalidParticles.load(std::memory_order_relaxed);
		}

	private:
		Particles::ParticleProperty* positions() const{
			return _positions;
		}

		Particles::ParticleProperty* refPositions() const{
			return _refPositions;
		}

		Particles::ParticleProperty* identifiers() const{
			return _identifiers;
		}

		Particles::ParticleProperty* refIdentifiers() const{
			return _refIdentifiers;
		}

		const SimulationCell& refCell() const{
			return _simCellRef;
		}

		bool computeStrain(
			std::size_t particleIndex,
			CutoffNeighborFinder& neighborFinder,
			const std::vector<std::size_t>& refToCurrentIndexMap,
			const std::vector<std::size_t>& currentToRefIndexMap
		);

		Particles::ParticleProperty* _positions;
		Particles::ParticleProperty* _refPositions;
		Particles::ParticleProperty* _identifiers;
		Particles::ParticleProperty* _refIdentifiers;

		SimulationCell _simCell;
		SimulationCell _simCellRef; 

		AffineTransformation _currentSimCellInv;
		AffineTransformation _reducedToAbsolute;

		double _cutoff;
		bool _eliminateCellDeformation;
		bool _assumeUnwrappedCoordinates;
		bool _calculateDeformationGradients;
		bool _calculateStrainTensors;
		bool _calculateNonaffineSquaredDisplacements;
		bool _calculatePolarDecomposition;

		std::shared_ptr<Particles::ParticleProperty> _shearStrains;
		std::shared_ptr<Particles::ParticleProperty> _volumetricStrains;
		std::shared_ptr<Particles::ParticleProperty> _nonaffineSquaredDisplacements;
		std::shared_ptr<Particles::ParticleProperty> _invalidParticles;
		std::shared_ptr<Particles::ParticleProperty> _strainTensors;
		std::shared_ptr<Particles::ParticleProperty> _deformationGradients;
		std::shared_ptr<Particles::ParticleProperty> _rotationTensors;
		std::shared_ptr<Particles::ParticleProperty> _stretchTensors;

		std::atomic<std::size_t> _numInvalidParticles{0};
	};
};

}
