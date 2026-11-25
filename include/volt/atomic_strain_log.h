#pragma once

#include <volt/atomic_strain_engine.h>
#include <volt/core/run_log.h>

namespace Volt{

// The execution log of one atomic strain analysis. Each phase hands it what it built, and the log writes what that
// holds and how long the phase took. Scripts parse these lines, so their wording stays as it is.
class AtomicStrainLog : public RunLog{
public:
    using RunLog::dump;

    void dump(const AtomicStrainModifier::AtomicStrainEngine& engine);
    void dump(const AtomicStrainModifier::StrainSummary& summary);
};

}
