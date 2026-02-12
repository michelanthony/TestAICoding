#include "monod_model.hpp"

#include <stdexcept>

namespace monod {

namespace {
constexpr double kMuMax = 0.5;
constexpr double kKs = 0.1;
}  // namespace

double growth_rate(double substrate) {
    return kMuMax * substrate / (kKs + substrate);
}

double dX_dt(double biomass, double substrate) {
    return growth_rate(substrate) * biomass;
}

double dS_dt(double biomass, double substrate) {
    return -growth_rate(substrate) * biomass;
}

SimulationResult integrate(double x0, double s0, double t_max, double dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument("dt must be strictly positive");
    }
    if (t_max <= 0.0) {
        throw std::invalid_argument("t_max must be strictly positive");
    }

    const std::size_t steps = static_cast<std::size_t>(t_max / dt);
    if (steps == 0) {
        throw std::invalid_argument("t_max / dt must be >= 1");
    }

    SimulationResult result;
    result.time_points.resize(steps);
    result.biomass.resize(steps);
    result.substrate.resize(steps);

    result.time_points[0] = 0.0;
    result.biomass[0] = x0;
    result.substrate[0] = s0;

    for (std::size_t i = 1; i < steps; ++i) {
        result.time_points[i] = i * dt;
        result.biomass[i] =
            result.biomass[i - 1] + dt * dX_dt(result.biomass[i - 1], result.substrate[i - 1]);
        result.substrate[i] =
            result.substrate[i - 1] + dt * dS_dt(result.biomass[i - 1], result.substrate[i - 1]);
    }

    return result;
}

}  // namespace monod
