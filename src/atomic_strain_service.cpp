#include <volt/atomic_strain_service.h>
#include <volt/atomic_strain_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/json_utils.h>
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

    auto positions = FrameAdapter::createPositionPropertyShared(currentFrame);
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

    auto refPositions = FrameAdapter::createPositionPropertyShared(refFrame);

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

    auto shear = engine.shearStrains();
    auto volumetric = engine.volumetricStrains();
    auto strainProp = engine.strainTensors();
    auto defgrad = engine.deformationGradients();
    auto D2minProp = engine.nonaffineSquaredDisplacements();
    auto invalid = engine.invalidParticles();

    size_t n = currentFrame.positions.size();

    double totalShear = 0.0;
    double totalVolumetric = 0.0;
    double maxShear = 0.0;
    int count = 0;

    for(size_t i = 0; i < n; i++){
        if(shear){
            double s = shear->getDouble(i);
            totalShear += s;
            if(s > maxShear) maxShear = s;
        }
        if(volumetric) totalVolumetric += volumetric->getDouble(i);
        count++;
    }

    json root;
    root["main_listing"] = {
        { "cutoff", _cutoff },
        { "num_invalid_particles", engine.numInvalidParticles() },
        { "average_shear_strain", count > 0 ? totalShear / count : 0.0 },
        { "average_volumetric_strain", count > 0 ? totalVolumetric / count : 0.0 },
        { "max_shear_strain", maxShear }
    };

    json perAtom = json::array();
    for(std::size_t i = 0; i < n; i++){
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
        } else {
            a["D2min"] = nullptr;
        }
        a["invalid"] = invalid ? (invalid->getInt(i) != 0) : false;

        perAtom.push_back(a);
    }

    root["per-atom-properties"] = perAtom;

    if(!outputFilename.empty()){
        const std::string outputPath = outputFilename + "_atomic_strain.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(root, outputPath, false)){
            spdlog::info("Atomic strain msgpack written to {}", outputPath);
        }else{
            spdlog::warn("Could not write atomic strain msgpack: {}", outputPath);
        }

        // --- atoms.msgpack (AtomisticExporter) ---
        // Canonical per-atom envelope grouped by "valid"/"invalid" so the
        // Volt viewport can render the strain-annotated atoms. Mirrors
        // OVITO's AtomicStrainModifier which publishes a SelectionProperty
        // for invalid particles plus StrainTensor / DeformationGradient /
        // D2min user properties.
        json atomsByBucket = json::object();
        json validAtoms = json::array();
        json invalidAtoms = json::array();
        int validCount = 0;
        int invalidCount = 0;
        for(std::size_t i = 0; i < n; i++){
            const Point3& pos = currentFrame.positions[i];
            const bool isInvalid = invalid ? (invalid->getInt(i) != 0) : false;
            const int structureId = isInvalid ? 1 : 0;
            const char* structureName = isInvalid ? "INVALID" : "VALID";
            json atom = {
                {"id", currentFrame.ids[i]},
                {"pos", {pos.x(), pos.y(), pos.z()}},
                {"structure_id", structureId},
                {"structure_name", structureName},
                {"cluster_id", 0},
                {"shear_strain", shear ? shear->getDouble(i) : 0.0},
                {"volumetric_strain", volumetric ? volumetric->getDouble(i) : 0.0},
                {"invalid", isInvalid}
            };
            if(strainProp){
                atom["strain_tensor"] = {
                    strainProp->getDoubleComponent(i, 0),
                    strainProp->getDoubleComponent(i, 1),
                    strainProp->getDoubleComponent(i, 2),
                    strainProp->getDoubleComponent(i, 5),
                    strainProp->getDoubleComponent(i, 4),
                    strainProp->getDoubleComponent(i, 3)
                };
            }
            if(defgrad){
                json grad = json::array();
                for(int c = 0; c < 9; c++)
                    grad.push_back(defgrad->getDoubleComponent(i, c));
                atom["deformation_gradient"] = grad;
            }
            if(D2minProp){
                atom["D2min"] = D2minProp->getDouble(i);
            }
            if(isInvalid){
                invalidAtoms.push_back(std::move(atom));
                ++invalidCount;
            }else{
                validAtoms.push_back(std::move(atom));
                ++validCount;
            }
        }
        if(validCount > 0) atomsByBucket["VALID"] = std::move(validAtoms);
        if(invalidCount > 0) atomsByBucket["INVALID"] = std::move(invalidAtoms);

        json structuresListing = json::array();
        if(validCount > 0){
            structuresListing.push_back({
                {"structure_id", 0}, {"structure_name", "VALID"}, {"atom_count", validCount}
            });
        }
        if(invalidCount > 0){
            structuresListing.push_back({
                {"structure_id", 1}, {"structure_name", "INVALID"}, {"atom_count", invalidCount}
            });
        }

        json exportWrapper;
        exportWrapper["main_listing"] = {
            {"total_atoms", static_cast<int>(n)},
            {"structure_count", static_cast<int>(structuresListing.size())},
            {"valid_atoms", validCount},
            {"invalid_atoms", invalidCount}
        };
        exportWrapper["sub_listings"] = { {"structures", structuresListing} };
        exportWrapper["export"] = json::object();
        exportWrapper["export"]["AtomisticExporter"] = atomsByBucket;
        const std::string atomsPath = outputFilename + "_atoms.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(exportWrapper, atomsPath, false)){
            spdlog::info("Exported atoms data to: {}", atomsPath);
        }else{
            spdlog::warn("Could not write atoms msgpack: {}", atomsPath);
        }
    }

    return root;
}


}
