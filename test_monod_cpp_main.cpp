#include "monod_model.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

bool almost_equal(double a, double b, double tolerance = 1e-12) {
    return std::fabs(a - b) <= tolerance;
}

bool compare_vectors(const std::vector<double>& got,
                     const std::vector<double>& expected,
                     const char* label) {
    if (got.size() != expected.size()) {
        std::cerr << label << " size mismatch. got=" << got.size()
                  << " expected=" << expected.size() << '\n';
        return false;
    }

    for (std::size_t i = 0; i < got.size(); ++i) {
        if (!almost_equal(got[i], expected[i])) {
            std::cerr << label << " mismatch at index " << i << ": got=" << got[i]
                      << " expected=" << expected[i] << '\n';
            return false;
        }
    }

    return true;
}

}  // namespace

int main() {
    // Python reference values computed with the same equations and parameters:
    // X0=0.1, S0=1.0, t_max=1.0, dt=0.1.
    const std::vector<double> expected_t = {
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    };
    const std::vector<double> expected_x = {
        0.1,
        0.104545454545,
        0.109295548849,
        0.114259294373,
        0.119446077777,
        0.124865674044,
        0.13052825973,
        0.136444426298,
        0.142625193468,
        0.149082022501,
    };
    const std::vector<double> expected_s = {
        1.0,
        0.995454545455,
        0.990704451151,
        0.985740705627,
        0.980553922223,
        0.975134325956,
        0.96947174027,
        0.963555573702,
        0.957374806532,
        0.950917977499,
    };

    const monod::SimulationResult result = monod::integrate(0.1, 1.0, 1.0, 0.1);

    const bool ok_t = compare_vectors(result.time_points, expected_t, "time_points");
    const bool ok_x = compare_vectors(result.biomass, expected_x, "biomass");
    const bool ok_s = compare_vectors(result.substrate, expected_s, "substrate");

    if (ok_t && ok_x && ok_s) {
        std::cout << "C++ Monod integration matches Python reference values." << std::endl;
        return 0;
    }

    return 1;
}
