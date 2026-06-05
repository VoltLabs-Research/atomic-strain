#include <volt/plugin/plugin_main.h>
#include <volt/plugin/output_serializer.h>
#include <volt/atomic_strain_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>

static const Volt::Plugin::PluginDescriptor descriptor{
    .name = "volt-atomic-strain",
    .description = "Atomic Strain Analysis",
    .options = {
        {"--cutoff", "float", "Cutoff radius for neighbor search", "3.0"},
        {"--eliminate_cell_deformation", "bool", "Eliminate cell deformation", "false"},
        {"--assume_unwrapped", "bool", "Assume unwrapped coordinates", "false"},
        {"--calc_deformation_gradient", "bool", "Compute deformation gradient F", "true"},
        {"--calc_strain_tensors", "bool", "Compute strain tensors", "true"},
        {"--calc_d2min", "bool", "Compute D2min (nonaffine displacement)", "true"},
    },
    .needsReferenceFrame = true
};

VOLT_PLUGIN_MAIN(descriptor,
    [](const auto& opts, const Volt::LammpsParser::Frame& frame,
       const Volt::LammpsParser::Frame* refFramePtr,
       const std::string& outputBase) -> Volt::Plugin::json {
        using namespace Volt;
        using namespace Volt::Particles;

        const LammpsParser::Frame& refFrame = refFramePtr ? *refFramePtr : frame;

        auto positions = FrameAdapter::createPositionPropertyShared(frame);
        if (!positions) return AnalysisResult::failure("Failed to create position property");

        if (frame.natoms != refFrame.natoms)
            return AnalysisResult::failure("Atom count mismatch between current and reference frames");

        auto refPositions = FrameAdapter::createPositionPropertyShared(refFrame);
        auto identifiers = FrameAdapter::createIdentifierProperty(frame);
        auto refIdentifiers = FrameAdapter::createIdentifierProperty(refFrame);

        double cutoff = CLI::getDouble(opts, "--cutoff", 3.0);
        bool eliminateCell = CLI::getBool(opts, "--eliminate_cell_deformation", false);
        bool assumeUnwrapped = CLI::getBool(opts, "--assume_unwrapped", false);
        bool calcDefGrad = CLI::getBool(opts, "--calc_deformation_gradient", true);
        bool calcStrain = CLI::getBool(opts, "--calc_strain_tensors", true);
        bool calcD2min = CLI::getBool(opts, "--calc_d2min", true);

        AtomicStrainModifier::AtomicStrainEngine engine(
            positions.get(), frame.simulationCell,
            refPositions.get(), refFrame.simulationCell,
            identifiers.get(), refIdentifiers.get(),
            cutoff, eliminateCell, assumeUnwrapped,
            calcDefGrad, calcStrain, calcD2min
        );
        engine.perform();

        auto shear = engine.shearStrains();
        auto volumetric = engine.volumetricStrains();
        auto strainProp = engine.strainTensors();
        auto defgrad = engine.deformationGradients();
        auto D2minProp = engine.nonaffineSquaredDisplacements();
        auto invalid = engine.invalidParticles();

        const std::size_t n = frame.positions.size();
        double totalShear = 0.0, totalVolumetric = 0.0, maxShear = 0.0;
        int count = 0;

        for (std::size_t i = 0; i < n; ++i) {
            if (shear) { double s = shear->getDouble(i); totalShear += s; if (s > maxShear) maxShear = s; }
            if (volumetric) totalVolumetric += volumetric->getDouble(i);
            ++count;
        }

        nlohmann::json result;
        result["main_listing"] = {
            {"cutoff", cutoff},
            {"num_invalid_particles", engine.numInvalidParticles()},
            {"average_shear_strain", count > 0 ? totalShear / count : 0.0},
            {"average_volumetric_strain", count > 0 ? totalVolumetric / count : 0.0},
            {"max_shear_strain", maxShear}
        };

        if (!outputBase.empty()) {
            const int fixedExtra = 3;
            const int strainExtra = strainProp ? 1 : 0;
            const int defgradExtra = defgrad ? 1 : 0;
            const int d2minExtra = D2minProp ? 1 : 0;

            Plugin::serializePluginOutput(outputBase, frame, result, {
                .summaryFileSuffix = "_atomic_strain",
                .bucketResolver = [&invalid](std::size_t i) {
                    return (invalid && invalid->getInt(i) != 0)
                        ? std::string("INVALID") : std::string("VALID");
                },
                .atomFieldWriter = [&](MsgpackWriter& w, std::size_t i, int& extraCount) {
                    extraCount = fixedExtra + strainExtra + defgradExtra + d2minExtra;
                    w.write_key("shear_strain"); w.write_double(shear ? shear->getDouble(i) : 0.0);
                    w.write_key("volumetric_strain"); w.write_double(volumetric ? volumetric->getDouble(i) : 0.0);
                    w.write_key("invalid"); w.write_bool(invalid ? (invalid->getInt(i) != 0) : false);
                    if (strainProp) {
                        w.write_key("strain_tensor"); w.write_array_header(6);
                        for (int c = 0; c < 6; ++c) w.write_double(strainProp->getDoubleComponent(i, c));
                    }
                    if (defgrad) {
                        w.write_key("deformation_gradient"); w.write_array_header(9);
                        for (int c = 0; c < 9; ++c) w.write_double(defgrad->getDoubleComponent(i, c));
                    }
                    if (D2minProp) {
                        w.write_key("D2min"); w.write_double(D2minProp->getDouble(i));
                    }
                }
            });
        }

        result["is_failed"] = false;
        return result;
    })
