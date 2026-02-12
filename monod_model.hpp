#pragma once

#include <vector>

namespace monod {

struct SimulationResult {
    std::vector<double> time_points;
    std::vector<double> biomass;
    std::vector<double> substrate;
};

double growth_rate(double substrate);
double dX_dt(double biomass, double substrate);
double dS_dt(double biomass, double substrate);
SimulationResult integrate(double x0, double s0, double t_max, double dt);

}  // namespace monod
