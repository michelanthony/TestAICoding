#include "monod_model.hpp"

#include <algorithm>
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

bool write_csv(const std::string& path,
               const std::vector<double>& expected_t,
               const std::vector<double>& expected_x,
               const std::vector<double>& expected_s,
               const monod::SimulationResult& result) {
    std::ofstream out(path.c_str());
    if (!out) {
        std::cerr << "Unable to write CSV file: " << path << '\n';
        return false;
    }

    out << "step,time,x_python,x_cpp,s_python,s_cpp\n";
    out << std::setprecision(15);
    for (std::size_t i = 0; i < result.time_points.size(); ++i) {
        out << i << ',' << expected_t[i] << ',' << expected_x[i] << ',' << result.biomass[i]
            << ',' << expected_s[i] << ',' << result.substrate[i] << '\n';
    }

    return true;
}

double scale_x(double value, double min_value, double max_value, double width, double margin) {
    const double range = (max_value - min_value);
    if (range <= 0.0) {
        return margin;
    }
    return margin + ((value - min_value) / range) * width;
}

double scale_y(double value, double min_value, double max_value, double height, double margin) {
    const double range = (max_value - min_value);
    if (range <= 0.0) {
        return margin + height;
    }
    return margin + height - ((value - min_value) / range) * height;
}

std::string polyline_points(const std::vector<double>& x_values,
                            const std::vector<double>& y_values,
                            double min_x,
                            double max_x,
                            double min_y,
                            double max_y,
                            double chart_width,
                            double chart_height,
                            double margin) {
    std::string points;

    for (std::size_t i = 0; i < x_values.size(); ++i) {
        const double x = scale_x(x_values[i], min_x, max_x, chart_width, margin);
        const double y = scale_y(y_values[i], min_y, max_y, chart_height, margin);
        points += std::to_string(x) + "," + std::to_string(y);
        if (i + 1U < x_values.size()) {
            points += " ";
        }
    }

    return points;
}

bool write_svg(const std::string& path,
               const std::vector<double>& expected_t,
               const std::vector<double>& expected_x,
               const std::vector<double>& expected_s,
               const monod::SimulationResult& result) {
    std::ofstream out(path.c_str());
    if (!out) {
        std::cerr << "Unable to write SVG file: " << path << '\n';
        return false;
    }

    const double svg_width = 980.0;
    const double svg_height = 430.0;
    const double margin = 40.0;
    const double chart_width = (svg_width - (3.0 * margin)) / 2.0;
    const double chart_height = svg_height - (2.0 * margin);

    const double min_t = expected_t.front();
    const double max_t = expected_t.back();

    const double min_x = std::min(*std::min_element(expected_x.begin(), expected_x.end()),
                                  *std::min_element(result.biomass.begin(), result.biomass.end()));
    const double max_x = std::max(*std::max_element(expected_x.begin(), expected_x.end()),
                                  *std::max_element(result.biomass.begin(), result.biomass.end()));

    const double min_s = std::min(*std::min_element(expected_s.begin(), expected_s.end()),
                                  *std::min_element(result.substrate.begin(), result.substrate.end()));
    const double max_s = std::max(*std::max_element(expected_s.begin(), expected_s.end()),
                                  *std::max_element(result.substrate.begin(), result.substrate.end()));

    const std::string x_py_points = polyline_points(expected_t, expected_x, min_t, max_t, min_x, max_x,
                                                    chart_width, chart_height, margin);
    const std::string x_cpp_points = polyline_points(result.time_points, result.biomass, min_t, max_t,
                                                     min_x, max_x, chart_width, chart_height, margin);

    std::vector<double> shifted_t(result.time_points.size());
    for (std::size_t i = 0; i < shifted_t.size(); ++i) {
        shifted_t[i] = result.time_points[i] + chart_width + margin;
    }

    const std::string s_py_points = polyline_points(shifted_t, expected_s,
                                                    min_t + chart_width + margin,
                                                    max_t + chart_width + margin,
                                                    min_s, max_s, chart_width, chart_height, margin);
    const std::string s_cpp_points = polyline_points(shifted_t, result.substrate,
                                                     min_t + chart_width + margin,
                                                     max_t + chart_width + margin,
                                                     min_s, max_s, chart_width, chart_height, margin);

    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << svg_width
        << "\" height=\"" << svg_height << "\" viewBox=\"0 0 " << svg_width << " "
        << svg_height << "\">\n";
    out << "  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    out << "  <text x=\"" << margin << "\" y=\"22\" font-size=\"16\">Biomass</text>\n";
    out << "  <text x=\"" << (2.0 * margin + chart_width)
        << "\" y=\"22\" font-size=\"16\">Substrate</text>\n";

    out << "  <rect x=\"" << margin << "\" y=\"" << margin << "\" width=\"" << chart_width
        << "\" height=\"" << chart_height << "\" fill=\"none\" stroke=\"#808080\"/>\n";
    out << "  <rect x=\"" << (2.0 * margin + chart_width) << "\" y=\"" << margin
        << "\" width=\"" << chart_width << "\" height=\"" << chart_height
        << "\" fill=\"none\" stroke=\"#808080\"/>\n";

    out << "  <polyline fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"2\" points=\""
        << x_py_points << "\"/>\n";
    out << "  <polyline fill=\"none\" stroke=\"#d62728\" stroke-width=\"2\" points=\""
        << x_cpp_points << "\"/>\n";

    out << "  <polyline fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"2\" points=\""
        << s_py_points << "\"/>\n";
    out << "  <polyline fill=\"none\" stroke=\"#d62728\" stroke-width=\"2\" points=\""
        << s_cpp_points << "\"/>\n";

    out << "  <text x=\"" << margin << "\" y=\"" << (svg_height - 10.0)
        << "\" font-size=\"12\" fill=\"#1f77b4\">Blue: Python reference</text>\n";
    out << "  <text x=\"" << (margin + 190.0) << "\" y=\"" << (svg_height - 10.0)
        << "\" font-size=\"12\" fill=\"#d62728\">Red: C++ result</text>\n";
    out << "</svg>\n";

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

    print_table(expected_x, expected_s, result);

    if (!ensure_results_dir()) {
        return 1;
    }

    const std::string csv_path = "results/monod_test_results.csv";
    const std::string svg_path = "results/monod_test_plot.svg";
    const bool csv_ok = write_csv(csv_path, expected_t, expected_x, expected_s, result);
    const bool svg_ok = write_svg(svg_path, expected_t, expected_x, expected_s, result);

    if (!(csv_ok && svg_ok)) {
        return 1;
    }

    std::cout << "\nGenerated table file: " << csv_path << '\n';
    std::cout << "Generated graph file: " << svg_path << '\n';

    if (ok_t && ok_x && ok_s) {
        std::cout << "C++ Monod integration matches Python reference values." << std::endl;
        return 0;
    }

    return 1;
}
