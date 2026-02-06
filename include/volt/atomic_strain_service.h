#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <nlohmann/json.hpp>
#include <string>

namespace Volt{

using json = nlohmann::json;

class AtomicStrainService{
public:
	AtomicStrainService();

	void setCutoff(double cutoff);
	void setReferenceFrame(const LammpsParser::Frame &ref);
	void setOptions(
		bool eliminateCellDeformation,
		bool assumeUnwrappedCoordinates,
		bool calculateDeformationGradient,
		bool calculateStrainTensors,
		bool calculateD2min
	);

	json compute(
		const LammpsParser::Frame& currentFrame,
		const std::string& outputFilename = ""
	);

private:
	double _cutoff;
	bool _eliminateCellDeformation;
	bool _assumeUnwrappedCoordinates;
	bool _calculateDeformationGradient;
	bool _calculateStrainTensors;
	bool _calculateD2min;

	bool _hasReference;
	LammpsParser::Frame _referenceFrame;

	json computeAtomicStrain(
		const LammpsParser::Frame& currentFrame,
		const LammpsParser::Frame& refFrame,
		Particles::ParticleProperty* positions,
		const std::string& outputFilename
	);
};
    
}
