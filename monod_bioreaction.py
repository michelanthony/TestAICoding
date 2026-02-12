import csv
import os

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
    if dt <= 0:
        raise ValueError("dt must be > 0")

    time_points = []
    X = []
    S = []

    t = 0.0
    x_val = float(X0)
    s_val = float(S0)

    while t < t_max:
        time_points.append(t)
        X.append(x_val)
        S.append(s_val)

        x_next = x_val + dt * dX_dt(x_val, s_val)
        s_next = s_val + dt * dS_dt(x_val, s_val)

        t += dt
        x_val = x_next
        s_val = s_next

    return time_points, X, S


def integrate_with_pybind11(X0, S0, t_max, dt):
    if monod_cpp is None:
        raise RuntimeError(
            "pybind11 extension unavailable. Build it with `pip install -e .`."
        )

    time_points, X, S = monod_cpp.integrate(X0, S0, t_max, dt)
    return list(time_points), list(X), list(S)


def write_tsv(path, time_points, X, S):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["step", "time", "x", "s"])
        for i, (t, x_val, s_val) in enumerate(zip(time_points, X, S)):
            writer.writerow([i, float(t), float(x_val), float(s_val)])


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

    output_path = "results/monod_python_results.tsv"
    write_tsv(output_path, time_points, X, S)
    print("Generated table file:", output_path)
