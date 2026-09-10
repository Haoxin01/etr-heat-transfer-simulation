#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simulate 2-D transient heat conduction in argon, a sample, and carbon paper.

The model uses an explicit finite-difference (FTCS) scheme in conservative
form. Thermal conductivity at each cell face is the arithmetic mean of the
conductivities at the two adjacent nodes. The measured temperature curve
``T_hot(t)`` is read from the supplied CSV file and imposed as a Dirichlet
boundary condition on both carbon-paper heaters.

Outputs
-------
- Temperature-field snapshots at 1, 2, 5, 6, 10, 20, 30, 40, 60, and 64 s.
- Temperatures at points 0.5, 1.0, and 1.5 cm above the upper carbon paper.
- CSV files sampled every 0.1 s.
- PNG and PDF plots of probe temperature versus time.
"""

from pathlib import Path
from time import time

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Plot settings
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


# Input and output settings
INPUT_FILE = Path("raw_Time_T_every_0p1s.csv")
SNAPSHOT_TIMES = np.array([1, 2, 5, 6, 10, 20, 30, 40, 60, 64], dtype=float)
CSV_TIMES = np.round(np.arange(0.1, SNAPSHOT_TIMES[-1] + 1e-12, 0.1), 10)

# Domain and mesh
DOMAIN_LENGTH = 0.05  # m; square domain (5 cm x 5 cm)
GRID_SIZE = 201

# Geometry
CARBON_PAPER_WIDTH = 0.015  # m (1.5 cm)
CARBON_PAPER_THICKNESS = 0.000120  # m (0.12 mm)
SAMPLE_THICKNESS = 0.0010  # m (1 mm)
PROBE_DISTANCES_CM = (0.5, 1.0, 1.5)

# Thermal properties: thermal conductivity [W/(m K)], density [kg/m^3],
# and specific heat capacity [J/(kg K)]
ARGON_PROPERTIES = (0.0177, 1.6, 520.0)
GEO2_PROPERTIES = (0.59, 3650.0, 590.0)
CARBON_PAPER_PROPERTIES = (0.25, 1700.0, 710.0)

# Numerical and boundary-condition settings
FOURIER_NUMBER = 0.4
AMBIENT_TEMPERATURE_C = 25.0
EFFECTIVE_HEAT_TRANSFER_COEFFICIENT = 10.0
PLOT_TEMPERATURE_LIMITS_C = (25.0, 1800.0)


def load_hot_temperature_curve(input_file):
    """Load and validate the measured hot-side temperature curve.

    The supplied reviewer dataset uses the columns ``time_s`` and ``T_hot_C``.
    The earlier Excel layout with columns ``Time`` and ``T`` remains supported.
    """
    input_file = Path(input_file)
    if not input_file.exists():
        raise FileNotFoundError(
            f"Input temperature profile not found: {input_file.resolve()}"
        )

    if input_file.suffix.lower() == ".csv":
        data = pd.read_csv(input_file)
    elif input_file.suffix.lower() in {".xlsx", ".xls"}:
        data = pd.read_excel(input_file, engine="openpyxl")
    else:
        raise ValueError("Input must be a CSV or Excel file.")

    if {"time_s", "T_hot_C"}.issubset(data.columns):
        curve_time_s = pd.to_numeric(data["time_s"], errors="raise").to_numpy()
        curve_temperature_c = pd.to_numeric(
            data["T_hot_C"], errors="raise"
        ).to_numpy()
    elif {"Time", "T"}.issubset(data.columns):
        curve_time_s = (
            pd.to_timedelta(data["Time"]).dt.total_seconds().to_numpy()
        )
        curve_temperature_c = pd.to_numeric(
            data["T"], errors="raise"
        ).to_numpy()
    else:
        raise ValueError(
            "Input must contain either 'time_s' and 'T_hot_C' columns or "
            "'Time' and 'T' columns."
        )

    if len(curve_time_s) < 2:
        raise ValueError("The input temperature profile requires at least two rows.")
    if not np.isfinite(curve_time_s).all() or not np.isfinite(
        curve_temperature_c
    ).all():
        raise ValueError(
            "The input temperature profile contains missing or non-finite values."
        )
    if np.any(np.diff(curve_time_s) <= 0):
        raise ValueError("Input time values must be strictly increasing.")
    if curve_time_s[-1] < SNAPSHOT_TIMES[-1]:
        raise ValueError(
            f"The input profile ends at {curve_time_s[-1]:g} s but the "
            f"simulation requires data through {SNAPSHOT_TIMES[-1]:g} s."
        )
    if curve_time_s[0] > 0:
        print(
            f"Note: the first input time is {curve_time_s[0]:g} s. "
            "Its temperature is held constant from t = 0 to the first record."
        )
    return curve_time_s, curve_temperature_c


def build_geometry():
    """Calculate grid indices for the sample and carbon-paper heaters.

    The 0.12-mm carbon paper is thinner than the 0.25-mm grid spacing. Because
    its measured temperature is imposed as a Dirichlet condition, each heater
    is represented as a subgrid interface on the nearest grid row rather than
    as a 0.25-mm-thick resolved material layer. Its physical 0.12-mm thickness
    is retained for geometry reporting and plotting.
    """
    dx = DOMAIN_LENGTH / (GRID_SIZE - 1)

    x_start = int((DOMAIN_LENGTH / 2 - CARBON_PAPER_WIDTH / 2) / dx)
    x_end = int((DOMAIN_LENGTH / 2 + CARBON_PAPER_WIDTH / 2) / dx)

    # A one-row Dirichlet mask represents each unresolved 0.12-mm heater.
    carbon_layers = 1
    sample_layers = max(1, int(round(SAMPLE_THICKNESS / dx)))

    y_mid = GRID_SIZE // 2
    sample_y_start = y_mid - sample_layers // 2
    sample_y_end = sample_y_start + sample_layers - 1

    lower_cp_y_start = sample_y_start - carbon_layers
    lower_cp_y_end = sample_y_start - 1
    upper_cp_y_start = sample_y_end + 1
    upper_cp_y_end = sample_y_end + carbon_layers

    return {
        "dx": dx,
        "x_start": x_start,
        "x_end": x_end,
        "sample_y_start": sample_y_start,
        "sample_y_end": sample_y_end,
        "lower_cp_y_start": lower_cp_y_start,
        "lower_cp_y_end": lower_cp_y_end,
        "upper_cp_y_start": upper_cp_y_start,
        "upper_cp_y_end": upper_cp_y_end,
    }


def build_material_fields(geometry):
    """Create spatial fields for thermal conductivity, density, and heat capacity."""
    argon_k, argon_rho, argon_cp = ARGON_PROPERTIES
    conductivity = np.full((GRID_SIZE, GRID_SIZE), argon_k)
    density = np.full((GRID_SIZE, GRID_SIZE), argon_rho)
    heat_capacity = np.full((GRID_SIZE, GRID_SIZE), argon_cp)

    x_slice = slice(geometry["x_start"], geometry["x_end"] + 1)
    sample_slice = slice(
        geometry["sample_y_start"], geometry["sample_y_end"] + 1
    )
    lower_cp_slice = slice(
        geometry["lower_cp_y_start"], geometry["lower_cp_y_end"] + 1
    )
    upper_cp_slice = slice(
        geometry["upper_cp_y_start"], geometry["upper_cp_y_end"] + 1
    )

    for field, sample_value, carbon_value in zip(
        (conductivity, density, heat_capacity),
        GEO2_PROPERTIES,
        CARBON_PAPER_PROPERTIES,
    ):
        field[sample_slice, x_slice] = sample_value
        field[lower_cp_slice, x_slice] = carbon_value
        field[upper_cp_slice, x_slice] = carbon_value

    return conductivity, density, heat_capacity


def build_carbon_paper_masks(geometry):
    """Create Boolean masks for the lower and upper carbon-paper layers."""
    lower_mask = np.zeros((GRID_SIZE, GRID_SIZE), dtype=bool)
    upper_mask = np.zeros((GRID_SIZE, GRID_SIZE), dtype=bool)
    x_slice = slice(geometry["x_start"], geometry["x_end"] + 1)

    lower_mask[
        geometry["lower_cp_y_start"] : geometry["lower_cp_y_end"] + 1,
        x_slice,
    ] = True
    upper_mask[
        geometry["upper_cp_y_start"] : geometry["upper_cp_y_end"] + 1,
        x_slice,
    ] = True
    return lower_mask, upper_mask


def build_probe_specs(geometry):
    """Define probe locations above the upper surface of the upper carbon paper."""
    dx = geometry["dx"]
    upper_surface_y_m = (geometry["upper_cp_y_end"] + 1) * dx
    center_x_index = (geometry["x_start"] + geometry["x_end"]) // 2
    probe_specs = []

    for distance_cm in PROBE_DISTANCES_CM:
        probe_y_m = upper_surface_y_m + distance_cm / 100
        if probe_y_m > DOMAIN_LENGTH:
            probe_y_index = None
            print(
                f"Warning: probe {distance_cm} cm above the carbon paper is "
                "outside the computational domain; its result will be NaN."
            )
        else:
            probe_y_index = int(round(probe_y_m / dx))
            probe_y_index = min(max(probe_y_index, 0), GRID_SIZE - 1)

        distance_label = str(distance_cm).replace(".", "p")
        probe_specs.append(
            {
                "distance_cm": distance_cm,
                "column": f"T_above_CP_{distance_label}cm_C",
                "y_index": probe_y_index,
                "x_index": center_x_index,
            }
        )

    return probe_specs


def apply_robin_boundaries(temperature, conductivity, dx):
    """Apply convective Robin conditions to all four outer boundaries."""
    h_dx = EFFECTIVE_HEAT_TRANSFER_COEFFICIENT * dx
    ambient = AMBIENT_TEMPERATURE_C

    temperature[:, 0] = (
        conductivity[:, 0] * temperature[:, 1] + h_dx * ambient
    ) / (conductivity[:, 0] + h_dx)
    temperature[:, -1] = (
        conductivity[:, -1] * temperature[:, -2] + h_dx * ambient
    ) / (conductivity[:, -1] + h_dx)
    temperature[0, :] = (
        conductivity[0, :] * temperature[1, :] + h_dx * ambient
    ) / (conductivity[0, :] + h_dx)
    temperature[-1, :] = (
        conductivity[-1, :] * temperature[-2, :] + h_dx * ambient
    ) / (conductivity[-1, :] + h_dx)


def calculate_face_conductivities(conductivity):
    """Return arithmetic-mean conductivities at the four interior faces.

    For every interior node P, the face conductivity between P and a neighbor
    F is k_f = (k_P + k_F) / 2, as specified in Supplementary Eq. S2.
    """
    center = conductivity[1:-1, 1:-1]
    north = 0.5 * (center + conductivity[2:, 1:-1])
    south = 0.5 * (center + conductivity[:-2, 1:-1])
    east = 0.5 * (center + conductivity[1:-1, 2:])
    west = 0.5 * (center + conductivity[1:-1, :-2])
    return north, south, east, west


def save_temperature_snapshot(
    temperature,
    snapshot_time_s,
    hot_temperature_c,
    geometry,
    probe_specs,
):
    """Save one temperature-field snapshot as a PNG image."""
    dx = geometry["dx"]
    extent_cm = [0, DOMAIN_LENGTH * 100, 0, DOMAIN_LENGTH * 100]
    x_cp_cm = geometry["x_start"] * dx * 100
    carbon_width_cm = CARBON_PAPER_WIDTH * 100
    sample_y_cm = geometry["sample_y_start"] * dx * 100
    sample_height_cm = (
        geometry["sample_y_end"] - geometry["sample_y_start"] + 1
    ) * dx * 100
    carbon_thickness_cm = CARBON_PAPER_THICKNESS * 100
    lower_cp_y_cm = sample_y_cm - carbon_thickness_cm
    upper_cp_y_cm = sample_y_cm + sample_height_cm

    figure, axis = plt.subplots(figsize=(6, 5))
    image = axis.imshow(
        temperature,
        origin="lower",
        cmap="inferno",
        extent=extent_cm,
        vmin=PLOT_TEMPERATURE_LIMITS_C[0],
        vmax=PLOT_TEMPERATURE_LIMITS_C[1],
    )
    figure.colorbar(image, ax=axis, label="Temperature (degC)")

    axis.add_patch(
        patches.Rectangle(
            (x_cp_cm, sample_y_cm),
            carbon_width_cm,
            sample_height_cm,
            fill=False,
            edgecolor="black",
            linewidth=1.2,
        )
    )
    for carbon_y_cm in (lower_cp_y_cm, upper_cp_y_cm):
        axis.add_patch(
            patches.Rectangle(
                (x_cp_cm, carbon_y_cm),
                carbon_width_cm,
                carbon_thickness_cm,
                facecolor="black",
                edgecolor="black",
                linewidth=0.5,
            )
        )

    for probe in probe_specs:
        if probe["y_index"] is None:
            continue
        probe_x_cm = probe["x_index"] * dx * 100
        probe_y_cm = probe["y_index"] * dx * 100
        axis.plot(probe_x_cm, probe_y_cm, marker="o", markersize=3)
        axis.text(
            probe_x_cm + 0.12,
            probe_y_cm,
            f"+{probe['distance_cm']} cm",
            fontsize=8,
            verticalalignment="center",
        )

    axis.set_xlabel("x (cm)")
    axis.set_ylabel("y (cm)")
    axis.set_title(
        f"t = {snapshot_time_s:.0f} s | T_hot = {hot_temperature_c:.1f} degC"
    )
    figure.savefig(
        f"snapshot_{int(snapshot_time_s)}s.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)


def save_probe_plot(simulation_data, probe_specs):
    """Save probe-temperature histories as PNG and PDF files."""
    figure, axis = plt.subplots(figsize=(6, 4))
    for probe in probe_specs:
        axis.plot(
            simulation_data["time_s"],
            simulation_data[probe["column"]],
            label=f"{probe['distance_cm']} cm above carbon paper",
        )

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Temperature (degC)")
    axis.legend(frameon=False)
    axis.set_xlim(0, SNAPSHOT_TIMES[-1])
    figure.savefig(
        "probe_above_carbon_paper_T_vs_time.png",
        dpi=300,
        bbox_inches="tight",
    )
    figure.savefig("probe_above_carbon_paper_T_vs_time.pdf", bbox_inches="tight")
    plt.close(figure)


def run_simulation(curve_time_s, curve_temperature_c):
    """Run the transient simulation and return the sampled results."""
    geometry = build_geometry()
    dx = geometry["dx"]
    conductivity, density, heat_capacity = build_material_fields(geometry)
    k_north, k_south, k_east, k_west = calculate_face_conductivities(
        conductivity
    )
    lower_cp_mask, upper_cp_mask = build_carbon_paper_masks(geometry)
    probe_specs = build_probe_specs(geometry)

    minimum_volumetric_heat_capacity = (density * heat_capacity).min()
    maximum_conductivity = conductivity.max()
    dt = (
        FOURIER_NUMBER
        * dx**2
        * minimum_volumetric_heat_capacity
        / (4 * maximum_conductivity)
    )
    end_time_s = SNAPSHOT_TIMES[-1]
    number_of_steps = int(np.ceil(end_time_s / dt))
    print(f"Time step = {dt * 1e3:.4f} ms; total steps = {number_of_steps}")

    print("\nProbe positions above the upper carbon paper:")
    for probe in probe_specs:
        if probe["y_index"] is None:
            print(f"  +{probe['distance_cm']} cm: outside the domain")
        else:
            print(
                f"  +{probe['distance_cm']} cm: "
                f"x = {probe['x_index'] * dx * 100:.2f} cm, "
                f"y = {probe['y_index'] * dx * 100:.2f} cm"
            )

    print("\nSample upper-surface temperature (transient solution):")
    print(" t (s) | center (degC) | mean (degC)")
    print("-------------------------------------")

    temperature = np.full(
        (GRID_SIZE, GRID_SIZE), AMBIENT_TEMPERATURE_C, dtype=float
    )
    sample_center_x = (geometry["x_start"] + geometry["x_end"]) // 2
    sample_x_slice = slice(geometry["x_start"], geometry["x_end"] + 1)

    snapshot_index = 0
    csv_index = 0
    current_time_s = 0.0
    start_time = time()
    records = []

    for step in range(number_of_steps + 1):
        hot_temperature_c = np.interp(
            current_time_s, curve_time_s, curve_temperature_c
        )
        temperature[lower_cp_mask] = hot_temperature_c
        temperature[upper_cp_mask] = hot_temperature_c

        updated_temperature = temperature.copy()
        updated_temperature[1:-1, 1:-1] += dt / (
            density[1:-1, 1:-1] * heat_capacity[1:-1, 1:-1]
        ) * (
            k_north
            * (temperature[2:, 1:-1] - temperature[1:-1, 1:-1])
            + k_south
            * (temperature[:-2, 1:-1] - temperature[1:-1, 1:-1])
            + k_east
            * (temperature[1:-1, 2:] - temperature[1:-1, 1:-1])
            + k_west
            * (temperature[1:-1, :-2] - temperature[1:-1, 1:-1])
        ) / dx**2

        temperature = updated_temperature
        temperature[lower_cp_mask] = hot_temperature_c
        temperature[upper_cp_mask] = hot_temperature_c
        apply_robin_boundaries(temperature, conductivity, dx)

        while (
            csv_index < len(CSV_TIMES)
            and current_time_s >= CSV_TIMEScsv_index] - 1e-12
        ):
            output_time_s = CSV_TIMES[csv_index]
            sample_top_center_c = temperature[
                geometry["sample_y_end"], sample_center_x
            ]
            sample_top_mean_c = temperature[
                geometry["sample_y_end"], sample_x_slice
            ].mean()
            record = {
                "time_s": output_time_s,
                "T_hot_C": np.interp(
                    output_time_s, curve_time_s, curve_temperature_c
                ),
                "sample_top_center_C": sample_top_center_c,
                "sample_top_mean_C": sample_top_mean_c,
            }

            for probe in probe_specs:
                if probe["y_index"] is None:
                    record[probe["column"]] = np.nan
                else:
                    record[probe["column"]] = temperature[
                        probe["y_index"], probe["x_index"]
                    ]

            records.append(record)
            csv_index += 1

        if (
            snapshot_index < len(SNAPSHOT_TIMES)
            and current_time_s >= SNAPSHOT_TIMES[snapshot_index] - 1e-12
        ):
            snapshot_time_s = SNAPSHOT_TIMES[snapshot_index]
            sample_top_center_c = temperature[
                geometry["sample_y_end"], sample_center_x
            ]
            sample_top_mean_c = temperature[
                geometry["sample_y_end"], sample_x_slice
            ].mean()
            print(
                f"{snapshot_time_s:6.0f} | {sample_top_center_c:13.2f} | "
                f"{sample_top_mean_c:11.2f}"
            )
            save_temperature_snapshot(
                temperature,
                snapshot_time_s,
                hot_temperature_c,
                geometry,
                probe_specs,
            )
            snapshot_index += 1

        current_time_s += dt

    elapsed_time_s = time() - start_time
    print(
        f"\nSimulation completed in {elapsed_time_s:.2f} s "
        f"(dt = {dt * 1e3:.4f} ms; steps = {step})."
    )
    return pd.DataFrame(records), probe_specs


def export_results(simulation_data, probe_specs):
    """Export simulation tables and the probe-temperature plot."""
    simulation_data.to_csv(
        "simulated_sample_and_probe_Time_T_every_0p1s.csv", index=False
    )

    probe_columns = ["time_s", "T_hot_C"] + [
        probe["column"] for probe in probe_specs
    ]
    simulation_data[probe_columns].to_csv(
        "probe_above_carbon_paper_Time_T_every_0p1s.csv", index=False
    )
    save_probe_plot(simulation_data, probe_specs)

    print("Generated files:")
    print("  1. heater_profile_used_every_0p1s.csv")
    print("  2. simulated_sample_and_probe_Time_T_every_0p1s.csv")
    print("  3. probe_above_carbon_paper_Time_T_every_0p1s.csv")
    print("  4. probe_above_carbon_paper_T_vs_time.png")
    print("  5. probe_above_carbon_paper_T_vs_time.pdf")
    print("  6. snapshot_1s.png, snapshot_2s.png, snapshot_5s.png,")
    print("     snapshot_6s.png, snapshot_10s.png, snapshot_20s.png,")
    print("     snapshot_30s.png, snapshot_40s.png, snapshot_60s.png,")
    print("     snapshot_64s.png")


def main():
    """Load input data, run the simulation, and export all results."""
    curve_time_s, curve_temperature_c = load_hot_temperature_curve(INPUT_FILE)

    raw_data = pd.DataFrame(
        {
            "time_s": CSV_TIMES,
            "T_hot_C": np.interp(CSV_TIMES, curve_time_s, curve_temperature_c),
        }
    )
    raw_data.to_csv("heater_profile_used_every_0p1s.csv", index=False)

    simulation_data, probe_specs = run_simulation(
        curve_time_s, curve_temperature_c
    )
    export_results(simulation_data, probe_specs)


if __name__ == "__main__":
    main()
