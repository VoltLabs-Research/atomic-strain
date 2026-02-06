#include <volt/atomic_strain_service.h>
#include <volt/atomic_strain_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/concurrence/parallel_system.h>
#include <spdlog/spdlog.h>

namespace Volt{

using namespace Volt::Particles;

AtomicStrainService::AtomicStrainService()
    : _cutoff(0.10),
      _eliminateCellDeformation(false),
      _assumeUnwrappedCoordinates(false),
      _calculateDeformationGradient(true),
      _calculateStrainTensors(true),
      _calculateD2min(true),
      _hasReference(false){}


void AtomicStrainService::setCutoff(double cutoff){
    _cutoff = cutoff;
}

void AtomicStrainService::setReferenceFrame(const LammpsParser::Frame &ref){
    _referenceFrame = ref;
    _hasReference = true;
}

void AtomicStrainService::setOptions(
    bool eliminateCellDeformation,
    bool assumeUnwrappedCoordinates,
    bool calculateDeformationGradient,
    bool calculateStrainTensors,
    bool calculateD2min
){
    _eliminateCellDeformation = eliminateCellDeformation;
    _assumeUnwrappedCoordinates = assumeUnwrappedCoordinates;
    _calculateDeformationGradient = calculateDeformationGradient;
    _calculateStrainTensors = calculateStrainTensors;
    _calculateD2min = calculateD2min;
}

json AtomicStrainService::compute(const LammpsParser::Frame& currentFrame, const std::string &outputFilename){
    const LammpsParser::Frame &refFrame = _hasReference ? _referenceFrame : currentFrame;

    auto positions = FrameAdapter::createPositionProperty(currentFrame);
    if(!positions){
        return AnalysisResult::failure("Failed to create position property");
    }

    json result = computeAtomicStrain(currentFrame, refFrame, positions.get(), outputFilename);
    result["is_failed"] = false;
    return result;
}

json AtomicStrainService::computeAtomicStrain(
    const LammpsParser::Frame& currentFrame,
    const LammpsParser::Frame& refFrame,
    ParticleProperty* positions,
    const std::string& outputFilename
){
    if(currentFrame.natoms != refFrame.natoms){
        throw std::runtime_error("Cannot calculate atomic strain. Number of atoms in current and reference frames does not match.");
    }

    auto refPositions = std::make_shared<ParticleProperty>(
        refFrame.positions.size(),
        ParticleProperty::PositionProperty,
        3,
        false
    );

    for(std::size_t i = 0; i < refFrame.positions.size(); i++){
        refPositions->setPoint3(i, refFrame.positions[i]);
    }

    auto identifiers = FrameAdapter::createIdentifierProperty(currentFrame);
    auto refIdentifiers = FrameAdapter::createIdentifierProperty(refFrame);

    AtomicStrainModifier::AtomicStrainEngine engine(
        positions,
        currentFrame.simulationCell,
        refPositions.get(),
        refFrame.simulationCell,
        identifiers.get(),
        refIdentifiers.get(),
        _cutoff,
        _eliminateCellDeformation,
        _assumeUnwrappedCoordinates,
        _calculateDeformationGradient,
        _calculateStrainTensors,
        _calculateD2min
    );

    engine.perform();

    // calculate summary stats
    double totalShear = 0.0;
    double totalVolumetric = 0.0;
    double maxShear = 0.0;
    int count = 0;
    
    auto shear = engine.shearStrains();
    auto volumetric = engine.volumetricStrains();

    size_t n = currentFrame.positions.size();
    for(size_t i = 0; i < n; i++){
        if(shear){
            double s = shear->getDouble(i);
            totalShear += s;
            if(s > maxShear) maxShear = s;
        }

        if(volumetric){
            totalVolumetric += volumetric->getDouble(i);
        }
        count++;
    }

    json root;
    root["cutoff"] = _cutoff;
    root["num_invalid_particles"] = engine.numInvalidParticles();
    root["summary"] = {
        { "average_shear_strain", count > 0 ? totalShear / count : 0.0 },
        { "average_volumetric_strain", count > 0 ? totalVolumetric / count : 0.0 },
        { "max_shear_strain", maxShear }
    };

    if(!outputFilename.empty()){
        // TODO: Implement msgpack export in standalone package
        spdlog::warn("File output not yet implemented in standalone package. Returning inline JSON data.");
    }
    {
        root["atomic_strain"] = json::array();
        auto strainProp = engine.strainTensors();
        auto defgrad = engine.deformationGradients();
        auto D2minProp = engine.nonaffineSquaredDisplacements();
        auto invalid = engine.invalidParticles();

        // per atom properties
        for(std::size_t i = 0; i < currentFrame.positions.size(); i++){
            json a;
            a["id"] = currentFrame.ids[i];
            a["shear_strain"] = shear ? shear->getDouble(i) : 0.0;
            a["volumetric_strain"] = volumetric ? volumetric->getDouble(i) : 0.0;

            if(strainProp){
                double xx = strainProp->getDoubleComponent(i, 0);
                double yy = strainProp->getDoubleComponent(i, 1);
                double zz = strainProp->getDoubleComponent(i, 2);
                double yz = strainProp->getDoubleComponent(i, 3);
                double xz = strainProp->getDoubleComponent(i, 4);
                double xy = strainProp->getDoubleComponent(i, 5);
                a["strain_tensor"] = { xx, yy, zz, xy, xz, yz };
            }

            if(defgrad){
                double xx = defgrad->getDoubleComponent(i, 0);
                double yx = defgrad->getDoubleComponent(i, 1);
                double zx = defgrad->getDoubleComponent(i, 2);
                double xy = defgrad->getDoubleComponent(i, 3);
                double yy = defgrad->getDoubleComponent(i, 4);
                double zy = defgrad->getDoubleComponent(i, 5);
                double xz = defgrad->getDoubleComponent(i, 6);
                double yz = defgrad->getDoubleComponent(i, 7);
                double zz = defgrad->getDoubleComponent(i, 8);
                a["deformation_gradient"] = { xx, yx, zx, xy, yy, zy, xz, yz, zz };
            }

            if(D2minProp){
                a["D2min"] = D2minProp->getDouble(i);
            }else{
                a["D2min"] = nullptr;
            }

            if(invalid){
                a["invalid"] = (invalid->getInt(i) != 0);
            }else{
                a["invalid"] = false;
            }

            root["atomic_strain"].push_back(a);
        }
    }

    return root;
}

}
