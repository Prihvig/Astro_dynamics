import math
import matplotlib.pyplot as plt
import numpy as np

# === Constants ===
mu_sun = 1.32712440018e11  # Sun's gravitational parameter [km^3/s^2]
AU_km = 149597870.7        # Astronomical Unit in kilometers
g0 = 9.80665               # Standard gravity [m/s^2]

# === Planetary Orbit Radii (circular orbits assumed) ===
r_earth_km = 1 * AU_km     # Earth's average orbit radius
r_mars_km = 1.524 * AU_km  # Mars' average orbit radius

# === Spacecraft Info ===
m0 = 1000.0   # Initial spacecraft mass [kg]
Isp = 320.0   # Specific impulse [s]
Thrust = 500.0  # Thrust [N]

# === Plot flag ===
plot_enabled = False  # Set True to show plot, False to skip plotting

# === Hohmann Transfer Function (non-normalized) ===
def hohmann_transfer(r1_km, r2_km):
    a_transfer = (r1_km + r2_km) / 2

    # Orbital velocities
    v1 = math.sqrt(mu_sun / r1_km)
    v2 = math.sqrt(mu_sun / r2_km)

    v_trans_a = math.sqrt(2 * mu_sun / r1_km - mu_sun / a_transfer)
    v_trans_b = math.sqrt(2 * mu_sun / r2_km - mu_sun / a_transfer)

    delta_v1 = v_trans_a - v1
    delta_v2 = v2 - v_trans_b
    total_delta_v = abs(delta_v1) + abs(delta_v2)

    transfer_time_s = math.pi * math.sqrt(a_transfer**3 / mu_sun)

    return {
        "v1": v1,
        "v2": v2,
        "v_trans_a": v_trans_a,
        "v_trans_b": v_trans_b,
        "delta_v1": delta_v1,
        "delta_v2": delta_v2,
        "total_delta_v": total_delta_v,
        "transfer_time_s": transfer_time_s,
        "a_transfer": a_transfer,
        "r1": r1_km,
        "r2": r2_km
    }

# === Rocket Equation to Calculate Fuel Mass ===
def fuel_mass_required(delta_v_kmps, m0_kg, Isp_s):
    delta_v_mps = delta_v_kmps * 1000
    mf = m0_kg * (1 - math.exp(-delta_v_mps / (g0 * Isp_s)))
    return mf

# === Plotting function ===
def plot_hohmann_transfer(results):
    r1 = results['r1']
    r2 = results['r2']
    a = (r1 + r2) / 2
    e = (r2 - r1) / (r2 + r1)  # eccentricity

    theta = np.linspace(0, 2 * np.pi, 500)
    r_transfer = a * (1 - e**2) / (1 + e * np.cos(theta))

    # Orbits of Earth and Mars as circles
    earth_x = r1 * np.cos(theta)
    earth_y = r1 * np.sin(theta)
    mars_x = r2 * np.cos(theta)
    mars_y = r2 * np.sin(theta)

    transfer_x = r_transfer * np.cos(theta)
    transfer_y = r_transfer * np.sin(theta)

    plt.figure(figsize=(8,8))
    plt.plot(earth_x, earth_y, label="Earth Orbit")
    plt.plot(mars_x, mars_y, label="Mars Orbit")
    plt.plot(transfer_x, transfer_y, label="Hohmann Transfer Orbit", linestyle='--')
    plt.plot(0, 0, 'yo', label="Sun")
    plt.plot(r1, 0, 'bo', label="Earth")
    plt.plot(r2, 0, 'ro', label="Mars")
    plt.gca().set_aspect('equal', adjustable='box')

    margin = 0.1*r2
    plt.xlim(-r2-margin, r2+margin)
    plt.ylim(-r2-margin, r2+margin)

    plt.xlabel("Distance (km x 10^6)")
    plt.ylabel("Distance (km x 10^6)")
    plt.title("Hohmann Transfer from Earth to Mars")
    #plt.legend(loc = 'upper right', bbox_to_anchor = (1,0.5), fontsize=10)
    plt.grid(True)
    plt.show()

# === Run the Calculation ===
results = hohmann_transfer(r_earth_km, r_mars_km)
delta_v_total_kmps = results['total_delta_v']

# Calculate fuel mass required
fuel_mass = fuel_mass_required(delta_v_total_kmps, m0, Isp)

# Calculate acceleration and burn time
acceleration = Thrust / m0       # m/s²
burn_time = (delta_v_total_kmps * 1000) / acceleration  # seconds

# === Print Results ===
print("\nHohmann Transfer (Earth to Mars):")
print(f"  Earth Orbit Radius:     {r_earth_km:,.0f} km")
print(f"  Mars Orbit Radius:      {r_mars_km:,.0f} km\n")

print(f"  === Velocity info ===")
print(f"  Velocity at Earth orbit:       {results['v1']:.4f} km/s")
print(f"  Velocity at Mars orbit:        {results['v2']:.4f} km/s")
print(f"  Transfer velocity at Earth:    {results['v_trans_a']:.4f} km/s")
print(f"  Transfer velocity at Mars:     {results['v_trans_b']:.4f} km/s\n")

print(f"  === Delta-V info ===")
print(f"  Delta-V at Earth orbit:        {results['delta_v1']:.4f} km/s")
print(f"  Delta-V at Mars arrival:       {results['delta_v2']:.4f} km/s")
print(f"  Total Delta-V:                 {delta_v_total_kmps:.4f} km/s\n")

print(f"  === Time of Flight ===")
print(f"  Time of flight:                {results['transfer_time_s']:.2f} seconds")
print(f"  Time of flight:                {results['transfer_time_s'] / 86400:.2f} days")
print(f"  Time of flight:                {results['transfer_time_s'] / (86400*365.25):.2f} years\n")

print(f"  === Propulsion Details ===")
print(f"  Initial Mass (m0):             {m0:.2f} kg")
print(f"  Specific Impulse (Isp):        {Isp:.2f} s")
print(f"  Thrust:                        {Thrust:.2f} N")
print(f"  Fuel Mass Required:            {fuel_mass:.2f} kg")
print(f"  Acceleration:                  {acceleration:.4f} m/s^2")
print(f"  Burn Time:                     {burn_time:.2f} seconds ({burn_time/3600:.2f} hours)")
print(f"  Burn Time:                     {burn_time/86400:.2f} days")

# === Conditionally show plot ===
if plot_enabled:
    plot_hohmann_transfer(results)
# Note: The plot will show the orbits of Earth and Mars, along with the Hohmann transfer orbit.
# The plot will be displayed if the `plot_enabled` variable is set to True.