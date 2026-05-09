#include "reactive_transport_2d.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "=== 2D Reactive Transport Model Test ===" << std::endl;

    // Grid setup: 20x20 cells on 10m x 10m domain
    int nx = 21;  // 21 nodes = 20 cells
    int ny = 21;
    double Lx = 10.0;  // meters
    double Ly = 10.0;  // meters

    ReactiveTransport2D model(nx, ny, Lx, Ly);
    std::cout << "Grid: " << nx << "x" << ny << " nodes on " << Lx << "x" << Ly << " m domain" << std::endl;

    // Physical parameters
    double ux = 0.1;      // Velocity in x: 0.1 m/s (advection direction)
    double uy = 0.0;      // Velocity in y: 0 m/s
    double alphaL = 0.1;  // Longitudinal dispersivity: 0.1 m
    double alphaT = 0.01; // Transverse dispersivity: 0.01 m
    double Dm = 1e-4;     // Molecular diffusion: 1e-4 m²/s
    double k_decay = 0.01; // First-order decay: 0.01 1/s

    model.setVelocity(ux, uy);
    model.setDispersion(alphaL, alphaT, Dm);
    model.setReaction(k_decay);

    std::cout << "\nPhysical Parameters:" << std::endl;
    std::cout << "  Velocity: (" << ux << ", " << uy << ") m/s" << std::endl;
    std::cout << "  Dispersivity (L, T): (" << alphaL << ", " << alphaT << ") m" << std::endl;
    std::cout << "  Diffusion: " << Dm << " m²/s" << std::endl;
    std::cout << "  Decay rate: " << k_decay << " 1/s" << std::endl;

    // Initialize with Gaussian plume at center-left
    double x0 = 2.0;      // Plume center at x = 2m
    double y0 = 5.0;      // Plume center at y = 5m (domain center)
    double amplitude = 1.0;  // Peak concentration: 1.0 mg/L
    double sigma = 0.5;      // Plume width: 0.5 m

    model.initializeGaussianPlume(x0, y0, amplitude, sigma);
    std::cout << "\nInitial Gaussian Plume:" << std::endl;
    std::cout << "  Center: (" << x0 << ", " << y0 << ") m" << std::endl;
    std::cout << "  Amplitude: " << amplitude << " mg/L" << std::endl;
    std::cout << "  Sigma: " << sigma << " m" << std::endl;

    // Time stepping parameters
    double dt = 0.1;      // Time step: 0.1 seconds
    int n_steps = 100;    // 100 steps = 10 seconds total
    int output_interval = 10;  // Output every 10 steps

    std::cout << "\nTime Stepping:" << std::endl;
    std::cout << "  dt = " << dt << " s" << std::endl;
    std::cout << "  Total steps: " << n_steps << std::endl;
    std::cout << "  Total time: " << n_steps * dt << " s" << std::endl;

    // Prepare output file
    std::ofstream outfile("results/reactive_transport_2d_results.tsv");
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open results/reactive_transport_2d_results.tsv" << std::endl;
        return 1;
    }

    // Write header
    outfile << "Time(s)\tMaxConc(mg/L)\tTotalMass(mg)\tMass_Initial_Ratio" << std::endl;

    double initial_mass = model.computeTotalMass();
    std::cout << "\nInitial total mass: " << initial_mass << " mg" << std::endl;

    // Time-stepping loop
    std::cout << "\nTime stepping..." << std::endl;
    for (int step = 0; step <= n_steps; ++step) {
        double current_time = model.getCurrentTime();
        double max_conc = model.getMaxConcentration();
        double total_mass = model.computeTotalMass();
        double mass_ratio = (initial_mass > 1e-12) ? total_mass / initial_mass : 0.0;

        // Output to file
        outfile << std::scientific << std::setprecision(6);
        outfile << current_time << "\t" << max_conc << "\t" << total_mass 
                << "\t" << mass_ratio << std::endl;

        if (step % output_interval == 0) {
            std::cout << "Step " << std::setw(3) << step << " | t=" << std::setw(7) << current_time 
                      << "s | c_max=" << std::setw(10) << max_conc 
                      << " mg/L | mass=" << std::setw(10) << total_mass << " mg" << std::endl;
        }

        // Perform time step
        if (step < n_steps) {
            model.timeStepBackwardEuler(dt, 20, 1e-6);
        }
    }

    outfile.close();
    std::cout << "\nResults written to: results/reactive_transport_2d_results.tsv" << std::endl;

    // Final diagnostics
    std::cout << "\n=== Final Diagnostics ===" << std::endl;
    std::cout << "Final time: " << model.getCurrentTime() << " s" << std::endl;
    std::cout << "Max concentration: " << model.getMaxConcentration() << " mg/L" << std::endl;
    std::cout << "Final total mass: " << model.computeTotalMass() << " mg" << std::endl;
    std::cout << "Mass balance (final/initial): " 
              << (model.computeTotalMass() / initial_mass) << std::endl;

    // Mass balance check
    double final_mass = model.computeTotalMass();
    double mass_loss_percent = (1.0 - final_mass / initial_mass) * 100.0;
    std::cout << "\nMass loss due to decay: " << mass_loss_percent << "%" << std::endl;

    // Expected mass loss from first-order decay
    double expected_remaining = initial_mass * std::exp(-k_decay * n_steps * dt);
    std::cout << "Expected mass (exponential decay): " << expected_remaining << " mg" << std::endl;
    std::cout << "Actual final mass: " << final_mass << " mg" << std::endl;
    std::cout << "Relative error: " 
              << std::abs(final_mass - expected_remaining) / expected_remaining * 100.0 
              << "%" << std::endl;

    if (mass_loss_percent > 0 && mass_loss_percent < 50) {
        std::cout << "\n✓ Test PASSED: Mass balance is reasonable" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗ Test WARNING: Mass balance may be incorrect" << std::endl;
        return 1;
    }
}
