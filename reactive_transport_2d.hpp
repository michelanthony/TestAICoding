#ifndef REACTIVE_TRANSPORT_2D_HPP
#define REACTIVE_TRANSPORT_2D_HPP

#include <vector>
#include <Eigen/Sparse>

/**
 * @class ReactiveTransport2D
 * @brief Solves 2D reactive transport equations using finite volume method with backward Euler time stepping.
 * 
 * Equation solved:
 * dc/dt + u·∇c - ∇·(D∇c) + kc = 0
 * 
 * where:
 * - c: concentration field
 * - u = (ux, uy): velocity field
 * - D: anisotropic dispersion tensor (includes molecular diffusion)
 * - k: decay constant
 */
class ReactiveTransport2D {
public:
    /**
     * @brief Constructor
     * @param nx Number of grid points in x-direction
     * @param ny Number of grid points in y-direction
     * @param Lx Domain length in x-direction
     * @param Ly Domain length in y-direction
     */
    ReactiveTransport2D(int nx, int ny, double Lx, double Ly);

    /**
     * @brief Set velocity field components
     * @param ux Velocity in x-direction
     * @param uy Velocity in y-direction
     */
    void setVelocity(double ux, double uy);

    /**
     * @brief Set dispersion parameters
     * @param alphaL Longitudinal dispersivity
     * @param alphaT Transverse dispersivity
     * @param Dm Molecular diffusion coefficient
     */
    void setDispersion(double alphaL, double alphaT, double Dm);

    /**
     * @brief Set reaction (decay) constant
     * @param k First-order decay constant
     */
    void setReaction(double k);

    /**
     * @brief Initialize concentration field as a 2D Gaussian plume
     * @param x0 x-coordinate of plume center
     * @param y0 y-coordinate of plume center
     * @param amplitude Maximum concentration
     * @param sigma Standard deviation of Gaussian
     */
    void initializeGaussianPlume(double x0, double y0, double amplitude, double sigma);

    /**
     * @brief Perform one time step using backward Euler with Newton iteration
     * @param dt Time step size
     * @param max_newton Maximum Newton iterations
     * @param tol Convergence tolerance for Newton residual
     */
    void timeStepBackwardEuler(double dt, int max_newton = 10, double tol = 1e-6);

    /**
     * @brief Get maximum concentration value
     * @return Maximum concentration in the domain
     */
    double getMaxConcentration() const;

    /**
     * @brief Compute total mass (integral of concentration)
     * @return Total mass × area = Σ c_i × dx × dy
     */
    double computeTotalMass() const;

    /**
     * @brief Get current simulation time
     */
    double getCurrentTime() const { return current_time_; }

    /**
     * @brief Get concentration at grid point (i, j)
     */
    double getConcentration(int i, int j) const { return c_[nodeIndex(i, j)]; }

    /**
     * @brief Get grid spacing in x-direction
     */
    double getDx() const { return dx_; }

    /**
     * @brief Get grid spacing in y-direction
     */
    double getDy() const { return dy_; }

private:
    // Grid parameters
    int nx_, ny_;           // Grid dimensions
    double Lx_, Ly_;        // Domain dimensions
    double dx_, dy_;        // Grid spacing
    int n_dofs_;            // Total number of degrees of freedom

    // Flow parameters
    double ux_, uy_;        // Velocity components

    // Dispersion parameters
    double alphaL_;         // Longitudinal dispersivity
    double alphaT_;         // Transverse dispersivity
    double Dm_;             // Molecular diffusion coefficient

    // Reaction parameter
    double k_decay_;        // First-order decay constant

    // Solution storage
    std::vector<double> c_;         // Concentration field
    double current_time_;           // Current simulation time
    double dt_;                     // Time step (set in timeStepBackwardEuler)

    /**
     * @brief Convert 2D grid indices (i, j) to 1D DOF index (row-major ordering)
     */
    int nodeIndex(int i, int j) const;

    /**
     * @brief Compute anisotropic dispersion tensor components
     * @param Dxx [out] xx-component of dispersion tensor
     * @param Dyy [out] yy-component of dispersion tensor
     * @param Dxy [out] xy-component of dispersion tensor
     */
    void computeDispersionTensor(double& Dxx, double& Dyy, double& Dxy) const;

    /**
     * @brief Assemble sparse Jacobian matrix for backward Euler: dF/dc_new = I + dt*A
     */
    void assembleJacobianMatrix(double dt, Eigen::SparseMatrix<double>& jacobian);

    /**
     * @brief Compute residual: F(c_new) = (c_new - c_old)/dt + advection + diffusion + reaction
     */
    std::vector<double> computeResidual(const std::vector<double>& c_new,
                                        const std::vector<double>& c_old,
                                        double dt);
};

#endif // REACTIVE_TRANSPORT_2D_HPP
