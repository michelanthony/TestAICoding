#include "reactive_transport_2d.hpp"
#include <Eigen/SparseLU>
#include <iostream>
#include <algorithm>
#include <cmath>

ReactiveTransport2D::ReactiveTransport2D(int nx, int ny, double Lx, double Ly)
    : nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly),
      dx_(Lx / (nx - 1)), dy_(Ly / (ny - 1)),
      n_dofs_(nx * ny),
      ux_(0.0), uy_(0.0),
      alphaL_(0.0), alphaT_(0.0), Dm_(0.0),
      k_decay_(0.0),
      c_(n_dofs_, 0.0),
      current_time_(0.0) {
}

int ReactiveTransport2D::nodeIndex(int i, int j) const {
    // Row-major: index = i * ny + j
    return i * ny_ + j;
}

void ReactiveTransport2D::setVelocity(double ux, double uy) {
    ux_ = ux;
    uy_ = uy;
}

void ReactiveTransport2D::setDispersion(double alphaL, double alphaT, double Dm) {
    alphaL_ = alphaL;
    alphaT_ = alphaT;
    Dm_ = Dm;
}

void ReactiveTransport2D::setReaction(double k) {
    k_decay_ = k;
}

void ReactiveTransport2D::initializeGaussianPlume(double x0, double y0, 
                                                   double amplitude, double sigma) {
    std::fill(c_.begin(), c_.end(), 0.0);
    
    for (int i = 0; i < nx_; ++i) {
        for (int j = 0; j < ny_; ++j) {
            double x = i * dx_;
            double y = j * dy_;
            double dist_sq = (x - x0) * (x - x0) + (y - y0) * (y - y0);
            c_[nodeIndex(i, j)] = amplitude * std::exp(-dist_sq / (2.0 * sigma * sigma));
        }
    }
}

void ReactiveTransport2D::computeDispersionTensor(double& Dxx, double& Dyy, double& Dxy) const {
    double u_mag = std::sqrt(ux_ * ux_ + uy_ * uy_);
    
    if (u_mag < 1e-12) {
        // No velocity: pure diffusion
        Dxx = Dm_;
        Dyy = Dm_;
        Dxy = 0.0;
    } else {
        double u_ratio = 1.0 / u_mag;
        double ux_norm = ux_ * u_ratio;
        double uy_norm = uy_ * u_ratio;
        
        // D_ij = (alphaL - alphaT)*u_i*u_j/|u| + alphaT*|u|*delta_ij + Dm*delta_ij
        double coeff_aniso = (alphaL_ - alphaT_) * u_ratio;
        double coeff_iso = alphaT_ * u_mag + Dm_;
        
        Dxx = coeff_aniso * ux_norm * ux_norm + coeff_iso;
        Dyy = coeff_aniso * uy_norm * uy_norm + coeff_iso;
        Dxy = coeff_aniso * ux_norm * uy_norm;
    }
}

double ReactiveTransport2D::getMaxConcentration() const {
    return *std::max_element(c_.begin(), c_.end());
}

double ReactiveTransport2D::computeTotalMass() const {
    double total_mass = 0.0;
    for (int i = 0; i < nx_; ++i) {
        for (int j = 0; j < ny_; ++j) {
            total_mass += c_[nodeIndex(i, j)];
        }
    }
    return total_mass * dx_ * dy_;  // Grid cell area
}

void ReactiveTransport2D::assembleJacobianMatrix(double dt, 
                                                  Eigen::SparseMatrix<double>& jacobian) {
    // The Jacobian for backward Euler implicit: dF/dc_new = I + dt*A
    // where F(c_new) = (I + dt*A)*c_new - c_old
    // A is the finite volume operator (advection + dispersion + reaction)
    
    jacobian.resize(n_dofs_, n_dofs_);
    std::vector<Eigen::Triplet<double>> triplets;
    
    double Dxx, Dyy, Dxy;
    computeDispersionTensor(Dxx, Dyy, Dxy);
    
    for (int i = 0; i < nx_; ++i) {
        for (int j = 0; j < ny_; ++j) {
            int idx = nodeIndex(i, j);
            
            // Diagonal term: 1 + dt*(k + advection_div + dispersion_div)
            double diag_coeff = 1.0 + dt_ * k_decay_;
            
            // Upwind advection divergence contribution
            // For cell (i,j), use upwind scheme
            double adv_coeff = 0.0;
            
            // x-direction advection
            if (ux_ > 0 && i > 0) {
                // Upwind: use c[i-1,j], contributes to divergence at (i,j)
                adv_coeff += ux_ / dx_;
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i-1, j), -dt_ * ux_ / dx_));
            } else if (ux_ < 0 && i < nx_ - 1) {
                // Upwind: use c[i+1,j], contributes to divergence at (i,j)
                adv_coeff += -ux_ / dx_;
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i+1, j), -dt_ * (-ux_) / dx_));
            }
            
            // y-direction advection
            if (uy_ > 0 && j > 0) {
                // Upwind: use c[i,j-1]
                adv_coeff += uy_ / dy_;
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i, j-1), -dt_ * uy_ / dy_));
            } else if (uy_ < 0 && j < ny_ - 1) {
                // Upwind: use c[i,j+1]
                adv_coeff += -uy_ / dy_;
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i, j+1), -dt_ * (-uy_) / dy_));
            }
            
            // Dispersion: -div(D*grad(c))
            // Using centered differences for grad(c)
            // div_x = d/dx(Dxx*dc/dx + Dxy*dc/dy)
            // div_y = d/dy(Dxy*dc/dx + Dyy*dc/dy)
            
            double D_coeff = 0.0;
            
            if (i > 0 && i < nx_ - 1) {
                // d²c/dx² coefficient
                double coeff_xx = -dt_ * Dxx / (dx_ * dx_);
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i-1, j), coeff_xx));
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i+1, j), coeff_xx));
                D_coeff += 2.0 * Dxx / (dx_ * dx_);
            }
            
            if (j > 0 && j < ny_ - 1) {
                // d²c/dy² coefficient
                double coeff_yy = -dt_ * Dyy / (dy_ * dy_);
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i, j-1), coeff_yy));
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i, j+1), coeff_yy));
                D_coeff += 2.0 * Dyy / (dy_ * dy_);
            }
            
            // Mixed derivatives (Dxy terms)
            if (i > 0 && i < nx_ - 1 && j > 0 && j < ny_ - 1) {
                double coeff_xy = -dt_ * Dxy / (4.0 * dx_ * dy_);
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i-1, j-1), coeff_xy));
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i-1, j+1), -coeff_xy));
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i+1, j-1), -coeff_xy));
                triplets.push_back(Eigen::Triplet<double>(idx, nodeIndex(i+1, j+1), coeff_xy));
            }
            
            diag_coeff += dt_ * D_coeff;
            
            triplets.push_back(Eigen::Triplet<double>(idx, idx, diag_coeff));
        }
    }
    
    jacobian.setFromTriplets(triplets.begin(), triplets.end());
    jacobian.makeCompressed();
}

std::vector<double> ReactiveTransport2D::computeResidual(const std::vector<double>& c_new,
                                                          const std::vector<double>& c_old,
                                                          double dt) {
    std::vector<double> residual(n_dofs_, 0.0);
    
    double Dxx, Dyy, Dxy;
    computeDispersionTensor(Dxx, Dyy, Dxy);
    
    for (int i = 0; i < nx_; ++i) {
        for (int j = 0; j < ny_; ++j) {
            int idx = nodeIndex(i, j);
            
            // Residual: (c_new - c_old)/dt + advection + diffusion + reaction
            double res = (c_new[idx] - c_old[idx]) / dt + k_decay_ * c_new[idx];
            
            // Advection (upwind)
            if (ux_ > 0 && i > 0) {
                res += ux_ * (c_new[idx] - c_new[nodeIndex(i-1, j)]) / dx_;
            } else if (ux_ < 0 && i < nx_ - 1) {
                res += ux_ * (c_new[nodeIndex(i+1, j)] - c_new[idx]) / dx_;
            }
            
            if (uy_ > 0 && j > 0) {
                res += uy_ * (c_new[idx] - c_new[nodeIndex(i, j-1)]) / dy_;
            } else if (uy_ < 0 && j < ny_ - 1) {
                res += uy_ * (c_new[nodeIndex(i, j+1)] - c_new[idx]) / dy_;
            }
            
            // Dispersion (centered differences)
            if (i > 0 && i < nx_ - 1) {
                res -= Dxx * (c_new[nodeIndex(i+1, j)] - 2.0 * c_new[idx] + c_new[nodeIndex(i-1, j)]) / (dx_ * dx_);
            }
            
            if (j > 0 && j < ny_ - 1) {
                res -= Dyy * (c_new[nodeIndex(i, j+1)] - 2.0 * c_new[idx] + c_new[nodeIndex(i, j-1)]) / (dy_ * dy_);
            }
            
            // Mixed derivatives
            if (i > 0 && i < nx_ - 1 && j > 0 && j < ny_ - 1) {
                double d2c_dxdy = (c_new[nodeIndex(i+1, j+1)] - c_new[nodeIndex(i-1, j+1)] 
                                 - c_new[nodeIndex(i+1, j-1)] + c_new[nodeIndex(i-1, j-1)]) / (4.0 * dx_ * dy_);
                res -= Dxy * d2c_dxdy;
            }
            
            residual[idx] = res;
        }
    }
    
    return residual;
}

void ReactiveTransport2D::timeStepBackwardEuler(double dt, int max_newton, double tol) {
    // Save old solution
    std::vector<double> c_old = c_;
    
    // Newton iterations for: F(c_new) = (c_new - c_old)/dt + advection(c_new) + diffusion(c_new) + reaction(c_new) = 0
    for (int newton_iter = 0; newton_iter < max_newton; ++newton_iter) {
        // Compute residual
        std::vector<double> residual = computeResidual(c_, c_old, dt);
        
        // Compute residual norm
        double residual_norm = 0.0;
        for (double r : residual) {
            residual_norm += r * r;
        }
        residual_norm = std::sqrt(residual_norm);
        
        if (newton_iter % 5 == 0) {
            // Optionally print convergence (commented out for cleaner output)
            // std::cout << "  Newton iter " << newton_iter << ": ||F|| = " << residual_norm << std::endl;
        }
        
        if (residual_norm < tol) {
            break;  // Converged
        }
        
        // Assemble Jacobian
        Eigen::SparseMatrix<double> jacobian;
        assembleJacobianMatrix(dt, jacobian);
        
        // Convert residual to Eigen vector
        Eigen::VectorXd residual_eigen(n_dofs_);
        for (int i = 0; i < n_dofs_; ++i) {
            residual_eigen(i) = residual[i];
        }
        
        // Solve jacobian * delta_c = -residual
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(jacobian);
        Eigen::VectorXd delta_c = solver.solve(-residual_eigen);
        
        // Update solution
        for (int i = 0; i < n_dofs_; ++i) {
            c_[i] += delta_c(i);
        }
    }
    
    // Enforce non-negativity (optional physical constraint)
    for (double& conc : c_) {
        if (conc < 0.0) {
            conc = 0.0;
        }
    }
    
    current_time_ += dt;
}
