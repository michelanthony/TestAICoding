#include "monod_model.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
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

bool ensure_results_dir() {
#ifdef _WIN32
    const int code = std::system("if not exist results mkdir results");
#else
    const int code = std::system("mkdir -p results");
#endif
    if (code != 0) {
        std::cerr << "Unable to create results directory." << '\n';
        return false;
    }
    return true;
}

void print_table(const std::vector<double>& expected_x,
                 const std::vector<double>& expected_s,
                 const monod::SimulationResult& result) {
    std::cout << "\nComparison table (Python reference vs C++):\n";
    std::cout << std::left << std::setw(8) << "step" << std::setw(12) << "time"
              << std::setw(18) << "x_python" << std::setw(18) << "x_cpp"
              << std::setw(18) << "s_python" << std::setw(18) << "s_cpp" << '\n';

    std::cout << std::string(92, '-') << '\n';

    for (std::size_t i = 0; i < result.time_points.size(); ++i) {
        std::cout << std::fixed << std::setprecision(6) << std::left << std::setw(8) << i
                  << std::setw(12) << result.time_points[i] << std::setw(18) << expected_x[i]
                  << std::setw(18) << result.biomass[i] << std::setw(18) << expected_s[i]
                  << std::setw(18) << result.substrate[i] << '\n';
    }
}

bool write_tsv(const std::string& path,
               const std::vector<double>& expected_t,
               const std::vector<double>& expected_x,
               const std::vector<double>& expected_s,
               const monod::SimulationResult& result) {
    std::ofstream out(path.c_str());
    if (!out) {
        std::cerr << "Unable to write TSV file: " << path << '\n';
        return false;
    }

    out << "step\ttime\tx_python\tx_cpp\ts_python\ts_cpp\n";
    out << std::setprecision(15);
    for (std::size_t i = 0; i < result.time_points.size(); ++i) {
        out << i << '\t' << expected_t[i] << '\t' << expected_x[i] << '\t' << result.biomass[i]
            << '\t' << expected_s[i] << '\t' << result.substrate[i] << '\n';
    }

    return true;
}

}  // namespace

int main() {
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

    print_table(expected_x, expected_s, result);

    if (!ensure_results_dir()) {
        return 1;
    }

    const std::string tsv_path = "results/monod_test_results.tsv";
    const bool tsv_ok = write_tsv(tsv_path, expected_t, expected_x, expected_s, result);
    if (!tsv_ok) {
        return 1;
    }

    std::cout << "\nGenerated table file: " << tsv_path << '\n';

    if (ok_t && ok_x && ok_s) {
        std::cout << "C++ Monod integration matches Python reference values." << std::endl;
        return 0;
    }

    return 1;
}
