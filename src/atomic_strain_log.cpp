#include <volt/atomic_strain_log.h>

#include <spdlog/spdlog.h>

#include <cstddef>

namespace Volt{

void AtomicStrainLog::dump(const AtomicStrainModifier::AtomicStrainEngine& engine){
    const std::size_t n = engine.shearStrains()->size();
    const std::size_t invalid = engine.numInvalidParticles();
    spdlog::info("Atomic strain computed for {} atoms in {}", n, phaseElapsed());
    if(invalid > 0){
        spdlog::warn(
            "{} atoms ({:.1f}%) have fewer than 3 reference neighbors within the cutoff or a degenerate neighborhood; their strain is zero",
            invalid, percentOf(invalid, n));
    }
}

void AtomicStrainLog::dump(const AtomicStrainModifier::StrainSummary& summary){
    if(summary.count > 0)
        spdlog::info("Shear strain (von Mises): mean {:.6f}, max {:.6f}; volumetric strain: mean {:.6f}",
            summary.meanShear(), summary.maxShear, summary.meanVolumetric());
}

}
