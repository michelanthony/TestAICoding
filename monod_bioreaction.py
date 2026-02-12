import numpy as np
import matplotlib.pyplot as plt

# Monod kinetics parameters
# S: substrate concentration, X: biomass concentration
# mu_max: maximum specific growth rate, K_s: half-saturation constant
mu_max = 0.5  # 1/h
K_s = 0.1  # g/L

try:
    import monod_cpp
except ImportError:
    monod_cpp = None


# Differential equations for Monod kinetics
def dX_dt(X, S):
    mu = mu_max * S / (K_s + S)
    return mu * X


def dS_dt(X, S):
    mu = mu_max * S / (K_s + S)
    return -1 * mu * X


# Numerical integration using the Euler method
def integrate(X0, S0, t_max, dt):
    time_points = np.arange(0, t_max, dt)
    X = np.zeros(len(time_points))
    S = np.zeros(len(time_points))
    X[0] = X0
    S[0] = S0

    for i in range(1, len(time_points)):
        X[i] = X[i - 1] + dt * dX_dt(X[i - 1], S[i - 1])
        S[i] = S[i - 1] + dt * dS_dt(X[i - 1], S[i - 1])

    return time_points, X, S


def integrate_with_pybind11(X0, S0, t_max, dt):
    if monod_cpp is None:
        raise RuntimeError(
            "pybind11 extension unavailable. Build it with `pip install -e .`."
        )

    time_points, X, S = monod_cpp.integrate(X0, S0, t_max, dt)
    return np.array(time_points), np.array(X), np.array(S)


def plot_results(time_points, X, S):
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(time_points, X, label="Biomass Concentration (X)")
    plt.xlabel("Time (hours)")
    plt.ylabel("Concentration (g/L)")
    plt.title("Biomass Growth Over Time")
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(time_points, S, label="Substrate Concentration (S)", color="orange")
    plt.xlabel("Time (hours)")
    plt.ylabel("Concentration (g/L)")
    plt.title("Substrate Depletion Over Time")
    plt.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Initial conditions
    X0 = 0.1  # g/L
    S0 = 1.0  # g/L
    # Time parameters
    t_max = 10  # hours
    dt = 0.1  # hours

    # Integrate the system (prefer pybind11 module when available)
    if monod_cpp is not None:
        time_points, X, S = integrate_with_pybind11(X0, S0, t_max, dt)
    else:
        time_points, X, S = integrate(X0, S0, t_max, dt)

    plot_results(time_points, X, S)
