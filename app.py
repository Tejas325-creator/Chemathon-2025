import streamlit as st
import math
import matplotlib.pyplot as plt
import numpy as np

# --- Constants ---
L_D_RATIO = 2.5
LINING_THICKNESS_M = 0.003
INSULATION_THICKNESS_M = 0.025

DESIGN_PRESSURE_MPA = 1.5
ALLOWABLE_STRESS_MPA = 140.0
WELD_EFFICIENCY_J = 1.0
MIN_SHELL_THICKNESS_M = 0.003
CORROSION_ALLOWANCE_M = 0.0015

COIL_COVERAGE_FRACTION = 0.75
COIL_PIPE_OD_M = 0.025
COIL_PIPE_THICKNESS_M = 0.003
COIL_PITCH_FACTOR = 1.25

PROCESS_TEMP_C = 5.0
AMBIENT_TEMP_C = 25.0
K_LINING = 1.0
K_STEEL = 54.0
K_PUF = 0.022
H_INNER = 1000.0
H_OUTER = 10.0

K_INSULATION_STANDARD = 0.040
MIN_SHELL_THICKNESS_STANDARD_M = 0.005
DENSITY_STEEL_KG_M3 = 7850.0

# --- Functions ---
def calculate_reactor_geometry(volume_L):
    volume_m3 = volume_L / 1000.0
    inner_diameter_m = (4 * volume_m3 / (L_D_RATIO * math.pi)) ** (1/3)
    height_m = L_D_RATIO * inner_diameter_m
    return {"volume_L": volume_L, "inner_diameter_m": inner_diameter_m, "height_m": height_m}

def calculate_shell_thickness(inner_diameter_m, min_thickness_m):
    Di_mm = inner_diameter_m * 1000
    t_calculated_m = (DESIGN_PRESSURE_MPA * Di_mm) / \
                     (2 * ALLOWABLE_STRESS_MPA * WELD_EFFICIENCY_J - 1.2 * DESIGN_PRESSURE_MPA) / 1000
    base_thickness_m = max(t_calculated_m, min_thickness_m)
    total_steel_shell_m = base_thickness_m + CORROSION_ALLOWANCE_M
    return {"total_steel_shell_m": total_steel_shell_m}

def calculate_coil_dimensions(inner_diameter_m, height_m):
    pitch_m = COIL_PITCH_FACTOR * COIL_PIPE_OD_M
    covered_height_m = COIL_COVERAGE_FRACTION * height_m
    num_turns = math.ceil(covered_height_m / pitch_m)
    total_length_m = num_turns * math.pi * inner_diameter_m
    coil_inner_dia_m = COIL_PIPE_OD_M - (2 * COIL_PIPE_THICKNESS_M)
    return {"pipe_outer_diameter_m": COIL_PIPE_OD_M,
            "pipe_inner_diameter_m": coil_inner_dia_m,
            "pitch_m": pitch_m,
            "number_of_turns": num_turns,
            "total_length_m": total_length_m}

def calculate_heat_loss(geometry, shell, insulation_k):
    r1 = geometry['inner_diameter_m'] / 2
    r2 = r1 + LINING_THICKNESS_M
    r3 = r2 + shell['total_steel_shell_m']
    r4 = r3 + INSULATION_THICKNESS_M
    H = geometry['height_m']
    A_in = math.pi * geometry['inner_diameter_m'] * H
    A_out = math.pi * (2 * r4) * H
    R_conv_in = 1 / (H_INNER * A_in)
    R_lining = math.log(r2/r1) / (2 * math.pi * K_LINING * H)
    R_steel = math.log(r3/r2) / (2 * math.pi * K_STEEL * H)
    R_insulation = math.log(r4/r3) / (2 * math.pi * insulation_k * H)
    R_conv_out = 1 / (H_OUTER * A_out)
    R_total = R_conv_in + R_lining + R_steel + R_insulation + R_conv_out
    delta_T = AMBIENT_TEMP_C - PROCESS_TEMP_C
    heat_loss_W = delta_T / R_total
    return heat_loss_W

def calculate_shell_weight(geometry, shell):
    r_inner_shell_m = (geometry['inner_diameter_m'] / 2) + LINING_THICKNESS_M
    r_outer_shell_m = r_inner_shell_m + shell['total_steel_shell_m']
    vol_cyl_shell_m3 = math.pi * (r_outer_shell_m**2 - r_inner_shell_m**2) * geometry['height_m']
    area_end_cap_m2 = math.pi * r_outer_shell_m**2
    vol_end_caps_m3 = 2 * area_end_cap_m2 * shell['total_steel_shell_m']
    total_vol_m3 = vol_cyl_shell_m3 + vol_end_caps_m3
    weight_kg = total_vol_m3 * DENSITY_STEEL_KG_M3
    return weight_kg

def plot_comparison(our_metrics, standard_metrics):
    labels = ['Heat Loss (W)', 'Shell Weight (kg)']
    our_values = [our_metrics['heat_loss'], our_metrics['weight']]
    standard_values = [standard_metrics['heat_loss'], standard_metrics['weight']]
    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    rects1 = ax.bar(x - width/2, our_values, width, label='Our Design', color='#004488')
    rects2 = ax.bar(x + width/2, standard_values, width, label='Industry Standard', color='#DDAA33')
    ax.set_ylabel('Value')
    ax.set_title('Design Performance Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    ax.bar_label(rects1, padding=3, fmt='%.2f')
    ax.bar_label(rects2, padding=3, fmt='%.2f')
    return fig

# --- Streamlit UI ---
st.title("⚙️ Reaction Vessel Design & Comparison Calculator")
st.write("Designs a mild steel, PUF-insulated reaction vessel for the dye and textile industry.")

volume_input = st.number_input("Enter desired reactor volume (Liters):", min_value=1.0, value=500.0, step=10.0)

if st.button("Calculate Design"):
    geometry = calculate_reactor_geometry(volume_input)
    shell = calculate_shell_thickness(geometry['inner_diameter_m'], MIN_SHELL_THICKNESS_M)
    coil = calculate_coil_dimensions(geometry['inner_diameter_m'], geometry['height_m'])
    our_heat_loss = calculate_heat_loss(geometry, shell, K_PUF)
    our_weight = calculate_shell_weight(geometry, shell)

    standard_shell = calculate_shell_thickness(geometry['inner_diameter_m'], MIN_SHELL_THICKNESS_STANDARD_M)
    standard_heat_loss = calculate_heat_loss(geometry, standard_shell, K_INSULATION_STANDARD)
    standard_weight = calculate_shell_weight(geometry, standard_shell)

    st.subheader("Reactor Geometry (Our Model)")
    st.write(f"**Target Volume:** {geometry['volume_L']:.2f} L")
    st.write(f"**Inner Diameter:** {geometry['inner_diameter_m']*1000:.2f} mm")
    st.write(f"**Height:** {geometry['height_m']*1000:.2f} mm")
    st.write(f"**H/Di Ratio:** {L_D_RATIO:.1f}")

    st.subheader("Shell & Lining")
    st.write(f"**Inner Lining:** {LINING_THICKNESS_M*1000:.1f} mm (TiO₂+Epoxy)")
    st.write(f"**Total Steel Thickness:** {shell['total_steel_shell_m']*1000:.2f} mm")

    st.subheader("Cooling Coil (SS-316L Half-Pipe)")
    st.write(f"**Pipe Outer Diameter:** {coil['pipe_outer_diameter_m']*1000:.1f} mm")
    st.write(f"**Total Length:** {coil['total_length_m']:.2f} m")
    st.write(f"**Number of Turns:** {coil['number_of_turns']}")
    st.write(f"**Vertical Pitch:** {coil['pitch_m']*1000:.2f} mm")

    st.subheader("Insulation Performance")
    st.write(f"**Insulation Material:** Polyurethane Foam (PUF)")
    st.write(f"**Insulation Thickness:** {INSULATION_THICKNESS_M*1000:.1f} mm")
    st.write(f"**Thermal Conductivity:** {K_PUF} W/m·K")
    st.write(f"**Calculated Heat Loss:** {our_heat_loss:.2f} W (ΔT = {AMBIENT_TEMP_C - PROCESS_TEMP_C}°C)")

    fig = plot_comparison({'heat_loss': our_heat_loss, 'weight': our_weight},
                          {'heat_loss': standard_heat_loss, 'weight': standard_weight})
    st.pyplot(fig)

    heat_loss_reduction = (standard_heat_loss - our_heat_loss) / standard_heat_loss * 100
    weight_reduction = (standard_weight - our_weight) / standard_weight * 100

    st.subheader("Conclusion")
    st.success(f"Our design reduces heat loss by {heat_loss_reduction:.1f}% and weight by {weight_reduction:.1f}% compared to the sample industry standard.")

