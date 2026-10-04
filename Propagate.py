import json
import math

def propagate_body(body: dict, t_days: float = 0.0) -> tuple[float, float, float]:
    """
    Propagates a planet or relay to `t_days` (days since 2026-09-22 00:00:00 TDB).
    If you have elapsed time in hours `h`, pass `t_days = h / 24.0`.
    Returns (x, y, z) in AU in the Sun-centered J2000 ecliptic frame.
    """
    # 0. Extract orbital elements from JSON dict
    a = float(body["a_au"])
    e = float(body["e"])
    i = math.radians(float(body["i_deg"]))
    Omega = math.radians(float(body["node_deg"]))
    omega = math.radians(float(body["arg_peri_deg"]))
    M0_deg = float(body["mean_anomaly_deg"])
    n_deg_day = float(body["mean_motion_deg_day"])

    # 1. Advance mean anomaly M(t) = (M0 + n * t) mod 360°
    M_deg = (M0_deg + n_deg_day * t_days) % 360.0
    M = math.radians(M_deg)

    # 2. Solve Kepler's equation: E - e * sin(E) = M using Newton-Raphson
    E = M
    for _ in range(100):
        f_val = E - e * math.sin(E) - M
        if abs(f_val) <= 1e-15:
            break
        E = E - f_val / (1.0 - e * math.cos(E))

    # 3. Compute 2D orbital-plane coordinates (x', y')
    x_prime = a * (math.cos(E) - e)
    y_prime = a * math.sqrt(1.0 - e * e) * math.sin(E)

    # 4. Rotate into 3D J2000 ecliptic frame by omega, then i, then Omega
    cos_w, sin_w = math.cos(omega), math.sin(omega)
    cos_i, sin_i = math.cos(i), math.sin(i)
    cos_O, sin_O = math.cos(Omega), math.sin(Omega)

    # Rotate by omega
    x1 = x_prime * cos_w - y_prime * sin_w
    y1 = x_prime * sin_w + y_prime * cos_w

    # Rotate by i
    x2 = x1
    y2 = y1 * cos_i
    z2 = y1 * sin_i

    # Rotate by Omega
    x = x2 * cos_O - y2 * sin_O
    y = x2 * sin_O + y2 * cos_O
    z = z2

    return x, y, z


# --- Example Usage ---
if __name__ == "__main__":
    with open("orbital_data/orbital_elements.json") as f:
        orb_data = json.load(f)
    with open("orbital_data/network_model.json") as f:
        net_data = json.load(f)

    # Combine the 9 settlements and 2 relays
    all_bodies = orb_data["bodies"] + net_data["relays"]

    # Propagate at t = 0 days (Epoch)
    print("Positions at t = 0:")
    for body in all_bodies:
        x, y, z = propagate_body(body, t_days=0.0)
        print(f"{body['name']:8s}: x = {x:10.6f}, y = {y:10.6f}, z = {z:10.6f}")