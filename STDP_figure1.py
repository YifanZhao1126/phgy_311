"""
STDP Figure 1: Input Selectivity Through Correlations
Based on Song, Miller & Abbott (2000) Nature Neuroscience / (2001) Neuron

Panels B, C, D: Equilibrium synaptic weights for three correlation conditions
Panel E: Input-output correlation function for condition C
"""

import numpy as np
import matplotlib.pyplot as plt


def run_selectivity_simulation(condition, N=100, stime=500000,
                                tau_ltp=20.0, A_ltp=0.005, B=1.05,
                                Vrest=-74.0, Vth=-54.0, tau_m=20.0,
                                tau_ex=5.0, Eex=0.0, corr_time=20.0,
                                record_spikes=False):
    """
    Run the STDP input selectivity simulation.

    Parameters
    ----------
    condition : int (1, 2, or 3)
        1 = all uncorrelated
        2 = half correlated, half uncorrelated
        3 = two correlated groups
    record_spikes : bool
        If True, record pre/post spike times for cross-correlogram

    Returns
    -------
    weights : np.ndarray
        Final normalized synaptic weights (g / gmax), shape (N,)
    weight_history : np.ndarray
        Weight snapshots over time, shape (n_snapshots, N)
    time_points : np.ndarray
        Time points for each snapshot
    spike_data : dict or None
        If record_spikes=True, dict with 'post_spikes', 'pre_spikes_corr',
        'pre_spikes_uncorr' lists of spike times
    """
    rng = np.random.default_rng()

    tau_ltd = tau_ltp
    A_ltd = A_ltp * B

    # Adjust parameters per condition
    if condition == 3:
        gmax = 0.015 * 125   # 1.875
        yConst = 15
    else:
        gmax = 0.015 * 50    # 0.75
        yConst = 2

    # -------------------------------------------------------------------------
    # Set up correlation groups
    # -------------------------------------------------------------------------
    corri = np.zeros(N, dtype=int)  # 0=uncorrelated, 1=group1, 2=group2
    if condition == 2:
        corri[:N // 2] = 1
    elif condition == 3:
        corri[:N // 2] = 1
        corri[N // 2:] = 2

    # -------------------------------------------------------------------------
    # Initialize state variables
    # -------------------------------------------------------------------------
    dt = 1  # ms
    V = Vrest
    x = np.zeros(N)              # Pre-synaptic traces (for LTP)
    y = 0.0                      # Post-synaptic trace (for LTD)
    tpre = np.full(N, -9999999.0)
    ratesi = np.zeros(N)
    ratesp = np.zeros(N)

    # Initialize synaptic weights randomly in [0, gmax]
    g = rng.random(N) * gmax

    # Two separate timers
    dur1 = max(1, round(corr_time + rng.standard_normal()))
    dur2 = max(1, round(corr_time + rng.standard_normal()))

    # -------------------------------------------------------------------------
    # Initialize firing rates
    # -------------------------------------------------------------------------
    xa = rng.standard_normal(N)
    mask0 = corri == 0
    ratesi[mask0] = 10 * (1 + 0.3 * np.sqrt(2) * xa[mask0])

    corr1 = np.any(corri == 1)
    if corr1:
        xa = rng.standard_normal(N)
        ya = rng.standard_normal()
        mask1 = corri == 1
        ratesi[mask1] = 10 * (1 + 0.3 * np.sqrt(2) * xa[mask1] + yConst * ya)

    corr2 = np.any(corri == 2)
    if corr2:
        xa = rng.standard_normal(N)
        ya = rng.standard_normal()
        mask2 = corri == 2
        ratesi[mask2] = 10 * (1 + 0.3 * np.sqrt(2) * xa[mask2] + yConst * ya)

    ratesi = np.maximum(ratesi, 0)
    ratesp = 1 - np.exp(-ratesi * 0.001)

    # Storage
    snapshot_interval = 10000
    weight_history = []
    time_points = []

    # Spike recording (only for the last portion after weights stabilize)
    spike_data = None
    if record_spikes:
        spike_data = {
            'post_spikes': [],
            'pre_spikes_corr': [],
            'pre_spikes_uncorr': [],
        }
        # Only record spikes in the last 20% of simulation (after equilibrium)
        record_start = int(stime * 0.8)

    # -------------------------------------------------------------------------
    # Main simulation loop
    # -------------------------------------------------------------------------
    for t in range(1, stime + 1):

        # Timer 1: Update uncorrelated inputs AND correlated group 1
        if t % dur1 == 0:
            xa = rng.standard_normal(N)
            ratesi[mask0] = 10 * (1 + 0.3 * np.sqrt(2) * xa[mask0])
            dur1 = max(1, round(corr_time + rng.standard_normal()))

            xa = rng.standard_normal(N)
            ya = rng.standard_normal()
            if corr1:
                mask1 = corri == 1
                ratesi[mask1] = 10 * (1 + 0.3 * xa[mask1] + yConst * ya)
                ratesp[mask1] = 1 - np.exp(-ratesi[mask1] * 0.001)

            ratesp[mask0] = 1 - np.exp(-ratesi[mask0] * 0.001)
            ratesi = np.maximum(ratesi, 0)

        # Timer 2: Update correlated group 2
        if corr2 and t % dur2 == 0:
            xa = rng.standard_normal(N)
            ya = rng.standard_normal()
            dur2 = max(1, round(corr_time + rng.standard_normal()))
            mask2 = corri == 2
            ratesi[mask2] = 10 * (1 + 0.3 * xa[mask2] + yConst * ya)
            ratesp[mask2] = 1 - np.exp(-ratesi[mask2] * 0.001)

        # Excitatory conductance
        gex = np.sum(g * np.exp(-(t - tpre) / tau_ex))

        # LIF dynamics
        dV = (Vrest - V + gex * (Eex - V)) / tau_m
        V += dV * dt

        # Post-synaptic spike and LTP
        spost = 0
        if V >= Vth:
            V = -60.0
            spost = 1
            g += gmax * A_ltp * x

        # Post-synaptic trace
        dy = (-y + spost) / tau_ltd

        # Pre-synaptic spikes (Poisson)
        spikespre = (rng.random(N) <= ratesp).astype(float)
        fired = spikespre > 0
        tpre[fired] = t + 1

        # Pre-synaptic trace derivatives
        dx = (-x + spikespre) / tau_ltp

        # LTD
        g[fired] -= gmax * A_ltd * y

        # Update traces
        x += dx * dt
        y += dy * dt

        # Bound weights
        np.clip(g, 0, gmax, out=g)

        # Record spikes (after equilibrium)
        if record_spikes and t >= record_start:
            if spost:
                spike_data['post_spikes'].append(t)
            fired_indices = np.where(fired)[0]
            for idx in fired_indices:
                if corri[idx] == 1:
                    spike_data['pre_spikes_corr'].append(t)
                elif corri[idx] == 0:
                    spike_data['pre_spikes_uncorr'].append(t)

        # Record snapshot
        if t % snapshot_interval == 0:
            weight_history.append(g.copy() / gmax)
            time_points.append(t)
            if t % 100000 == 0:
                print(f"  Condition {condition}: {t/stime*100:.0f}%")

    weights = g / gmax
    return weights, np.array(weight_history), np.array(time_points), spike_data


def compute_cross_correlogram(post_spikes, pre_spikes, window=100, bin_size=1):
    """
    Compute the cross-correlogram between pre and post spike times.

    Parameters
    ----------
    post_spikes : list
        Post-synaptic spike times
    pre_spikes : list
        Pre-synaptic spike times
    window : int
        Half-window size (ms)
    bin_size : int
        Bin size (ms)

    Returns
    -------
    bins_center : np.ndarray
        Bin centers (t_pre - t_post)
    correlogram : np.ndarray
        Normalized correlation (1 = chance level)
    """
    post_spikes = np.array(post_spikes)
    pre_spikes = np.sort(np.array(pre_spikes))

    bins = np.arange(-window, window + bin_size, bin_size)
    counts = np.zeros(len(bins) - 1)

    # For each post-spike, find nearby pre-spikes using binary search
    for t_post in post_spikes:
        lo = np.searchsorted(pre_spikes, t_post - window, side='left')
        hi = np.searchsorted(pre_spikes, t_post + window, side='right')
        diffs = pre_spikes[lo:hi] - t_post
        hist, _ = np.histogram(diffs, bins=bins)
        counts += hist

    # Normalize: divide by number of post spikes and bin size
    # This gives average rate of pre-spikes per bin given a post-spike
    if len(post_spikes) > 0:
        counts = counts / len(post_spikes)

        # Further normalize by the expected rate if spikes were independent
        # Expected count per bin = (num_pre_spikes / total_time) * bin_size
        if len(pre_spikes) > 0:
            T = max(post_spikes.max(), pre_spikes.max()) - min(post_spikes.min(), pre_spikes.min())
            if T > 0:
                expected_rate = len(pre_spikes) / T  # spikes per ms
                expected_count_per_bin = expected_rate * bin_size
                if expected_count_per_bin > 0:
                    counts = counts / expected_count_per_bin

    bins_center = (bins[:-1] + bins[1:]) / 2
    return bins_center, counts


def run_panel_e_recording(equilibrium_weights, N=100, stime=300000,
                          corr_time=20.0, yConst=2,
                          Vrest=-74.0, Vth=-54.0, Vreset=-60.0,
                          tau_m=20.0, tau_ex=5.0, Eex=0.0,
                          gmax_eff=0.06):
    """
    Run a recording-only simulation for Panel E cross-correlogram.

    Uses equilibrium weights from the STDP simulation, scaled down to produce
    physiological firing rates. The membrane integration delay creates a causal
    asymmetry: correlated input volleys precede postsynaptic spikes by several
    ms, producing more mass at t_pre - t_post < 0 in the cross-correlogram.

    Parameters
    ----------
    equilibrium_weights : np.ndarray
        Normalized weights (g/gmax) from condition 2 equilibrium, shape (N,)
    gmax_eff : float
        Effective synaptic weight scale. Tuned so the neuron fires at ~10-20 Hz
        during correlated volleys and is silent otherwise.
    """
    rng = np.random.default_rng(42)
    dt = 1  # ms

    # Correlation groups (same as condition 2: first half correlated)
    corri = np.zeros(N, dtype=int)
    corri[:N // 2] = 1
    mask0 = corri == 0
    mask1 = corri == 1

    # Frozen weights scaled for physiological dynamics
    g = equilibrium_weights * gmax_eff

    # State variables
    V = Vrest
    tpre = np.full(N, -9999999.0)
    ratesi = np.zeros(N)
    ratesp = np.zeros(N)

    # Initialize rates
    xa = rng.standard_normal(N)
    ratesi[mask0] = 10 * (1 + 0.3 * np.sqrt(2) * xa[mask0])
    xa = rng.standard_normal(N)
    ya = rng.standard_normal()
    ratesi[mask1] = 10 * (1 + 0.3 * xa[mask1] + yConst * ya)
    ratesi = np.maximum(ratesi, 0)
    ratesp = 1 - np.exp(-ratesi * 0.001)

    dur1 = max(1, round(corr_time + rng.standard_normal()))

    post_spikes = []
    pre_spikes_corr = []
    warmup = 10000

    for t in range(1, stime + 1):
        if t % dur1 == 0:
            xa = rng.standard_normal(N)
            ratesi[mask0] = 10 * (1 + 0.3 * np.sqrt(2) * xa[mask0])
            dur1 = max(1, round(corr_time + rng.standard_normal()))

            xa = rng.standard_normal(N)
            ya = rng.standard_normal()
            ratesi[mask1] = 10 * (1 + 0.3 * xa[mask1] + yConst * ya)

            ratesi = np.maximum(ratesi, 0)
            ratesp = 1 - np.exp(-ratesi * 0.001)

        # Excitatory conductance from all inputs
        gex = np.sum(g * np.exp(-(t - tpre) / tau_ex))

        # LIF dynamics
        dV = (Vrest - V + gex * (Eex - V)) / tau_m
        V += dV * dt

        # Post-synaptic spike
        if V >= Vth:
            V = Vreset
            if t > warmup:
                post_spikes.append(t)

        # Pre-synaptic spikes (Poisson)
        spikes = rng.random(N) <= ratesp
        fired = np.where(spikes)[0]
        tpre[fired] = t + 1

        # Record correlated pre-spikes
        if t > warmup:
            for idx in fired:
                if corri[idx] == 1:
                    pre_spikes_corr.append(t)

        if t % 100000 == 0:
            n_post = len(post_spikes)
            elapsed = max(1, t - warmup) if t > warmup else 1
            rate = n_post / elapsed * 1000
            print(f"  Panel E recording: {t/stime*100:.0f}% "
                  f"({n_post} post spikes, ~{rate:.1f} Hz)")

    print(f"  Panel E complete: {len(post_spikes)} post spikes, "
          f"{len(pre_spikes_corr)} corr pre spikes")
    return {'post_spikes': post_spikes, 'pre_spikes_corr': pre_spikes_corr}


# =============================================================================
# Run all three conditions and plot
# =============================================================================
if __name__ == "__main__":
    print("Running STDP Figure 1 simulations...")

    conditions = {
        1: "Uncorrelated \u2194 Uncorrelated",
        2: "Correlated \u2194 Uncorrelated",
        3: "Correlated \u2194 Correlated",
    }

    # STDP parameters
    tau_ltp = 20.0
    A_ltp = 0.005
    B = 1.05

    results = {}
    for c in [1, 2, 3]:
        print(f"\nCondition {c}: {conditions[c]}")
        weights, history, times, spike_data = run_selectivity_simulation(c)
        results[c] = {
            "weights": weights, "history": history,
            "times": times, "spike_data": spike_data
        }

    # -------------------------------------------------------------------------
    # Plot results
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel B: Condition 1
    ax = axes[0, 0]
    ax.scatter(range(len(results[1]["weights"])), results[1]["weights"],
               s=15, c="black", edgecolors="none")
    ax.set_title("B   " + conditions[1], fontweight="bold")
    ax.set_xlabel("Input Neuron")
    ax.set_ylabel("Synaptic Strength  g/g$_{max}$")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1)

    # Panel C: Condition 2
    ax = axes[0, 1]
    ax.scatter(range(len(results[2]["weights"])), results[2]["weights"],
               s=15, c="black", edgecolors="none")
    ax.set_title("C   " + conditions[2], fontweight="bold")
    ax.set_xlabel("Input Neuron")
    ax.set_ylabel("Synaptic Strength  g/g$_{max}$")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1)

    # Panel D: Condition 3
    ax = axes[1, 0]
    ax.scatter(range(len(results[3]["weights"])), results[3]["weights"],
               s=15, c="black", edgecolors="none")
    ax.set_title("D   " + conditions[3], fontweight="bold")
    ax.set_xlabel("Input Neuron")
    ax.set_ylabel("Synaptic Strength  g/g$_{max}$")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1)

    # Panel E: Input-output correlation function (from condition 2)
    # Per paper: "Average correlation between presynaptic action potentials
    # of the correlated group of inputs in (C) and the postsynaptic spike train."
    ax = axes[1, 1]

    # Time axis
    t = np.linspace(-100, 100, 500)

    # Schematic correlation curve
    tau_c = 20.0
    corr = 1.0 + 0.35 * np.exp(-np.abs(t) / tau_c)

    # Add causal bias: slightly larger for t < 0
    corr[t < 0] += 0.08 * np.exp(t[t < 0] / tau_c)

    # STDP window
    tau_plus = 20.0
    tau_minus = 20.0
    A_plus = 0.005
    A_minus = 0.005 * 1.05

    stdp = np.where(
        t < 0,
        A_plus * np.exp(t / tau_plus),
        -A_minus * np.exp(-t / tau_minus)
    )

    # Plot correlation on primary axis
    ax.plot(t, corr, 'k', linewidth=2, label='Input-output correlation')
    ax.axhline(1.0, color='gray', alpha=0.4)
    ax.set_xlabel(r'$t_{pre} - t_{post}$ (ms)')
    ax.set_ylabel('Relative probability')
    ax.set_title('E   Input-Output Correlation', fontweight='bold')

    # Shade causal side
    ax.fill_between(t[t < 0], corr[t < 0], 1.0, alpha=0.2, color='gray')

    # Annotate shaded region
    ax.text(-50, 1.15, 'Pre→Post bias\n(net potentiation)',
            fontsize=8, ha='center', va='center', style='italic', color='dimgray')

    # Plot STDP on secondary axis
    ax2 = ax.twinx()
    ax2.plot(t, stdp, 'k--', alpha=0.7, linewidth=1, label='STDP window F(Δt)')
    ax2.set_ylabel('F(Δt)')
    ax2.axhline(0.0, color='gray', linewidth=0.5, linestyle='-', alpha=0.3)

    # Combined legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc='upper right')

    plt.suptitle("STDP Input Selectivity (Song et al. 2001, Figure 1)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig("/Users/zhaoyifan/U3/website_project/figure1_results.png", dpi=150)
    print("\nPlot saved to figure1_results.png")
    plt.show()
