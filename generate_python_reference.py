"""Generate Python reference values used by the C++ integration test."""

MU_MAX = 0.5
K_S = 0.1


def dX_dt(x, s):
    mu = MU_MAX * s / (K_S + s)
    return mu * x


def dS_dt(x, s):
    mu = MU_MAX * s / (K_S + s)
    return -mu * x


def integrate(x0, s0, t_max, dt):
    steps = int(t_max / dt)
    t = [0.0] * steps
    x = [0.0] * steps
    s = [0.0] * steps
    x[0] = x0
    s[0] = s0
    for i in range(1, steps):
        t[i] = i * dt
        x[i] = x[i - 1] + dt * dX_dt(x[i - 1], s[i - 1])
        s[i] = s[i - 1] + dt * dS_dt(x[i - 1], s[i - 1])
    return t, x, s


if __name__ == "__main__":
    t_vals, x_vals, s_vals = integrate(0.1, 1.0, 1.0, 0.1)
    print("t=", [round(v, 12) for v in t_vals])
    print("x=", [round(v, 12) for v in x_vals])
    print("s=", [round(v, 12) for v in s_vals])
