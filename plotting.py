import matplotlib.pyplot as plt
import numpy as np
def plot_results(hist, p):
    plt.rcParams.update({
        'font.family': 'Courier New',
        'font.size': 14,
        'axes.titlesize': 14,
        'axes.labelsize': 14,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 12, 
        'figure.titlesize': 16
    })
    
    t_plot = p["time_steps"]
    theta_grid_profile = p["theta_grid_profile"]
    fig, axs = plt.subplots(2, 2, figsize=(16, 10))
    
    color_blue = '#1f77b4'
    color_orange = '#ff7f0e'
    color_green = '#2ca02c'
    color_red = '#d62728'
    axs[0, 0].plot(t_plot, theta_grid_profile, "--", color='gray', linewidth=2.0, label="Infinite Bus")
    axs[0, 0].plot(t_plot, hist["theta_inv"], color=color_orange, linewidth=2.5, label="Inverter Angle")
    axs[0, 0].plot(t_plot, hist["theta_pll"], color=color_blue, linewidth=2.0, label="PLL Angle")
    axs[0, 0].set_ylabel('Angle (rad)')
    axs[0, 0].set_title('PLL - Angle Tracking')
    axs[0, 1].plot(t_plot, hist["vd"], color=color_blue, linewidth=2.5, label="$v_d$")
    axs[0, 1].plot(t_plot, hist["vq"], color=color_orange, linewidth=2.5, label="$v_q$")
    axs[0, 1].set_ylabel('Voltage (pu)')
    axs[0, 1].set_title('DQ-Axis Voltages')
    Pref = np.full_like(t_plot, p["Pref_const"])
    axs[1, 0].plot(t_plot, Pref, "--", color=color_red, linewidth=2.0, label="Reference")
    axs[1, 0].plot(t_plot, hist["P_expr"], color=color_blue, linewidth=2.5, label="Delivered")
    axs[1, 0].set_ylabel('Active Power (pu)')
    axs[1, 0].set_xlabel('Time (s)')
    axs[1, 0].set_title('Active Power Tracking')
    Qref = np.full_like(t_plot, p["Qref_const"])
    axs[1, 1].plot(t_plot, Qref, "--", color=color_red, linewidth=2.0, label="Reference")
    axs[1, 1].plot(t_plot, hist["Q_expr"], color=color_green, linewidth=2.5, label="Delivered")
    axs[1, 1].set_ylabel('Reactive Power (pu)')
    axs[1, 1].set_xlabel('Time (s)')
    axs[1, 1].set_title('Reactive Power Tracking')
    for ax in axs.flat:
        ax.grid(True, which='major', linestyle='-', linewidth=0.75, alpha=0.25)
        ax.minorticks_on()
        ax.grid(True, which='minor', linestyle='-', linewidth=0.25, alpha=0.15)
        ax.set_axisbelow(True)
        ax.legend(loc='best', frameon=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    output_filename = 'inverter_dynamic_response_grid.png'
    # plt.savefig(output_filename, format='png', dpi=600)
    plt.show()
    # print(f"output saved '{output_filename}'")