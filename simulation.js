/**
 * STDP (Spike-Timing-Dependent Plasticity) Simulation
 * Based on Song, Miller & Abbott (2000) Nature Neuroscience
 *
 * This module contains the core simulation algorithms.
 */

// =============================================================================
// GAUSSIAN RANDOM NUMBER GENERATOR (Box-Muller Transform)
// =============================================================================

export class GaussianRNG {
    constructor() {
        this.hasSpare = false;
        this.spare = 0;
    }

    reset() {
        this.hasSpare = false;
        this.spare = 0;
    }

    next() {
        if (this.hasSpare) {
            this.hasSpare = false;
            return this.spare;
        }
        this.hasSpare = true;
        const u = Math.random();
        const v = Math.random();
        const mag = Math.sqrt(-2.0 * Math.log(u));
        this.spare = mag * Math.cos(2.0 * Math.PI * v);
        return mag * Math.sin(2.0 * Math.PI * v);
    }
}

// =============================================================================
// FIGURE 7 SIMULATION: Input Selectivity Through Correlations
// =============================================================================

/**
 * Simulates how STDP leads to input selectivity based on correlation structure.
 *
 * Three conditions tested:
 *   1. All inputs uncorrelated → weights remain distributed
 *   2. Half correlated, half uncorrelated → correlated inputs win
 *   3. Two correlated groups competing → winner-take-all
 *
 * @param {number} conditionNumber - 1, 2, or 3
 * @param {object} params - simulation parameters
 * @param {function} onProgress - callback for progress updates (time, weights)
 * @returns {Array} final synaptic weights (normalized by gmax)
 */
export async function runSelectivitySimulation(conditionNumber, params, onProgress = null) {
    const rng = new GaussianRNG();

    // Adjust parameters for condition 3 (stronger competition needed)
    if (conditionNumber === 3) {
        params.gmax = 0.015 * 125;  // 1.875
        params.yConst = 15;
    } else {
        params.gmax = 0.015 * 50;   // 0.75
        params.yConst = 2;
    }

    // -------------------------------------------------------------------------
    // Set up correlation groups
    // corri[i] = 0: uncorrelated, 1: correlated group 1, 2: correlated group 2
    // -------------------------------------------------------------------------
    const corri = new Array(params.N).fill(0);
    if (conditionNumber === 2) {
        for (let i = 0; i < params.N / 2; i++) corri[i] = 1;
    }
    if (conditionNumber === 3) {
        for (let i = 0; i < params.N / 2; i++) corri[i] = 1;
        for (let i = Math.floor(params.N / 2); i < params.N; i++) corri[i] = 2;
    }

    // -------------------------------------------------------------------------
    // Initialize state variables
    // -------------------------------------------------------------------------
    const dt = 1;                                      // Time step (ms)
    let V = params.Vrest;                              // Membrane potential (mV)
    const x = new Array(params.N).fill(0);             // Pre-synaptic traces (for LTP)
    let y = 0;                                         // Post-synaptic trace (for LTD)
    const tpre = new Array(params.N).fill(-9999999);   // Last pre-spike times
    const ratesi = new Array(params.N).fill(0);        // Instantaneous firing rates
    const ratesp = new Array(params.N).fill(0);        // Spike probabilities

    // Initialize synaptic weights randomly in [0, gmax]
    const g = new Array(params.N);
    for (let i = 0; i < params.N; i++) {
        g[i] = Math.min(Math.max(Math.random() * params.gmax, 0), params.gmax);
    }

    // Correlation timing variables
    let dur1 = Math.round(Math.max(1, params.corr_time + rng.next()));
    let dur2 = Math.round(Math.max(1, params.corr_time + rng.next()));

    // Check which correlation groups exist
    const hasCorr1 = corri.some(c => c === 1);
    const hasCorr2 = corri.some(c => c === 2);

    // -------------------------------------------------------------------------
    // Initialize firing rates based on correlation structure
    // -------------------------------------------------------------------------
    let xa = new Array(params.N);

    // Uncorrelated inputs: independent noise
    for (let i = 0; i < params.N; i++) xa[i] = rng.next();
    for (let i = 0; i < params.N; i++) {
        if (corri[i] === 0) {
            ratesi[i] = 10 * (1 + 0.3 * Math.sqrt(2) * xa[i]);
        }
    }

    // Correlated group 1: shared + private noise
    if (hasCorr1) {
        for (let i = 0; i < params.N; i++) xa[i] = rng.next();
        const ya = rng.next();  // Shared noise
        for (let i = 0; i < params.N; i++) {
            if (corri[i] === 1) {
                ratesi[i] = 10 * (1 + 0.3 * Math.sqrt(2) * xa[i] + params.yConst * ya);
            }
        }
    }

    // Correlated group 2: shared + private noise
    if (hasCorr2) {
        for (let i = 0; i < params.N; i++) xa[i] = rng.next();
        const ya = rng.next();
        for (let i = 0; i < params.N; i++) {
            if (corri[i] === 2) {
                ratesi[i] = 10 * (1 + 0.3 * Math.sqrt(2) * xa[i] + params.yConst * ya);
            }
        }
    }

    // Convert rates to spike probabilities
    for (let i = 0; i < params.N; i++) {
        ratesi[i] = Math.max(ratesi[i], 0);
        ratesp[i] = 1 - Math.exp(-ratesi[i] * 0.001);
    }

    // -------------------------------------------------------------------------
    // Main simulation loop
    // -------------------------------------------------------------------------
    for (let t = 1; t <= params.stime; t++) {

        // Update firing rates when correlation window expires
        if (t % dur1 === 0) {
            // Regenerate uncorrelated inputs
            for (let i = 0; i < params.N; i++) xa[i] = rng.next();
            for (let i = 0; i < params.N; i++) {
                if (corri[i] === 0) {
                    ratesi[i] = 10 * (1 + 0.3 * Math.sqrt(2) * xa[i]);
                }
            }
            dur1 = Math.round(Math.max(1, params.corr_time + rng.next()));

            // Regenerate correlated group 1
            if (hasCorr1) {
                for (let i = 0; i < params.N; i++) xa[i] = rng.next();
                const ya = rng.next();
                for (let i = 0; i < params.N; i++) {
                    if (corri[i] === 1) {
                        ratesi[i] = 10 * (1 + 0.3 * xa[i] + params.yConst * ya);
                        ratesp[i] = 1 - Math.exp(-ratesi[i] * 0.001);
                    }
                }
            }
            for (let i = 0; i < params.N; i++) {
                if (corri[i] === 0) ratesp[i] = 1 - Math.exp(-ratesi[i] * 0.001);
                ratesi[i] = Math.max(ratesi[i], 0);
            }
        }

        // Update correlated group 2 on separate timer
        if (hasCorr2 && t % dur2 === 0) {
            for (let i = 0; i < params.N; i++) xa[i] = rng.next();
            const ya = rng.next();
            dur2 = Math.round(Math.max(1, params.corr_time + rng.next()));
            for (let i = 0; i < params.N; i++) {
                if (corri[i] === 2) {
                    ratesi[i] = 10 * (1 + 0.3 * xa[i] + params.yConst * ya);
                    ratesp[i] = 1 - Math.exp(-ratesi[i] * 0.001);
                }
            }
        }

        // ---------------------------------------------------------------------
        // Calculate total excitatory conductance (sum of exponential kernels)
        // ---------------------------------------------------------------------
        let gex = 0;
        for (let i = 0; i < params.N; i++) {
            gex += g[i] * Math.exp(-(t - tpre[i]) / params.tau_ex);
        }

        // ---------------------------------------------------------------------
        // Leaky Integrate-and-Fire neuron dynamics
        // dV/dt = (Vrest - V + gex*(Eex - V)) / tau_m
        // ---------------------------------------------------------------------
        const dV = (params.Vrest - V + gex * (params.Eex - V)) / params.tau_m;
        V = V + dV * dt;

        // ---------------------------------------------------------------------
        // Post-synaptic spike and LTP
        // ---------------------------------------------------------------------
        let spost = 0;
        if (V >= params.Vth) {
            V = -60;  // Reset potential
            spost = 1;

            // LTP: if post fires after pre, strengthen synapse
            // Weight change proportional to pre-synaptic trace x[i]
            for (let i = 0; i < params.N; i++) {
                g[i] = g[i] + params.gmax * params.A_ltp * x[i];
            }
        }

        // Update post-synaptic trace (exponential decay)
        const dy = (-y + spost) / params.tau_ltd;

        // ---------------------------------------------------------------------
        // Pre-synaptic spikes (Poisson process)
        // ---------------------------------------------------------------------
        const spikespre = new Array(params.N);
        for (let i = 0; i < params.N; i++) {
            spikespre[i] = Math.random() <= ratesp[i] ? 1 : 0;
            if (spikespre[i]) tpre[i] = t + 1;
        }

        // Calculate dx (but don't update x yet, matching MATLAB order)
        const dx = new Array(params.N);
        for (let i = 0; i < params.N; i++) {
            dx[i] = (-x[i] + spikespre[i]) / params.tau_ltp;
        }

        // ---------------------------------------------------------------------
        // LTD: if pre fires after post, weaken synapse
        // Weight change proportional to post-synaptic trace y
        // ---------------------------------------------------------------------
        for (let i = 0; i < params.N; i++) {
            if (spikespre[i]) {
                g[i] = g[i] - params.gmax * params.A_ltd * y;
            }
        }

        // Update traces (after LTD, matching MATLAB order)
        for (let i = 0; i < params.N; i++) {
            x[i] = x[i] + dx[i] * dt;
        }
        y = y + dy * dt;

        // Bound weights to [0, gmax]
        for (let i = 0; i < params.N; i++) {
            g[i] = Math.max(Math.min(g[i], params.gmax), 0);
        }

        // Progress callback (every 50000 ms for UI updates)
        if (onProgress && t % 50000 === 0) {
            const normalizedWeights = g.map(w => w / params.gmax);
            await onProgress(t, params.stime, normalizedWeights);
        }
    }

    return g.map(w => w / params.gmax);
}

// =============================================================================
// FIGURE 4 SIMULATION: Temporal Coding Through Latency
// =============================================================================

/**
 * Simulates how STDP implements temporal coding based on input latency.
 *
 * Each input has a different latency (time of burst onset relative to t=0).
 * STDP strengthens early inputs (that fire before post-spike) and weakens
 * late inputs (that fire after post-spike), sharpening temporal selectivity.
 *
 * @param {object} params - simulation parameters
 * @param {function} onProgress - callback (iteration, total, weights, voltage)
 * @returns {object} { latencies, weights }
 */
export async function runLatencySimulation(params, onProgress = null) {
    const rng = new GaussianRNG();
    const dt = 1;
    const int_start = -50;   // Integration window start (ms)
    const int_end = 50;      // Integration window end (ms)
    const int_length = int_end - int_start + 1;

    // Generate input latencies (Gaussian distribution)
    const latencies = new Array(params.N);
    for (let i = 0; i < params.N; i++) {
        latencies[i] = rng.next() * params.latency_std;
    }

    // Initialize all weights equally
    const g = new Array(params.N).fill(0.003);

    // Spike probability during burst
    const firingp = 1 - Math.exp(-params.burstrate * 0.001);

    let initialVoltage = null;

    // -------------------------------------------------------------------------
    // Main loop: repeated stimulus presentations
    // -------------------------------------------------------------------------
    for (let trial = 1; trial <= params.ttimes; trial++) {
        const tpre = new Array(params.N).fill(-9999999);
        let V = params.Vrest;
        const x = new Array(params.N).fill(0);
        let y = 0;
        const Vs = new Array(int_length).fill(0);

        let idx = 0;
        for (let t = int_start; t <= int_end; t++) {
            // Calculate excitatory conductance
            let gex = 0;
            for (let i = 0; i < params.N; i++) {
                gex += g[i] * Math.exp(-(t - tpre[i]) / params.tau_ex);
            }

            // Integrate-and-fire dynamics
            const dV = (params.Vrest - V + gex * (params.Eex - V)) / params.tau_m;
            V = V + dV * dt;
            Vs[idx] = V;

            // Post-synaptic spike
            let spost = 0;
            if (V >= params.Vth) {
                V = -60;
                Vs[idx] = 0;
                spost = 1;

                // LTP for recently active synapses
                for (let i = 0; i < params.N; i++) {
                    g[i] = g[i] + params.gmax * params.A_ltp * x[i];
                }
            }

            const dy = (-y + spost) / params.tau_ltd;

            // Pre-synaptic spikes (only during each input's burst window)
            const spikespre = new Array(params.N);
            for (let i = 0; i < params.N; i++) {
                const inBurstWindow = t >= latencies[i] && t < (latencies[i] + params.burstdur);
                spikespre[i] = (Math.random() <= firingp && inBurstWindow) ? 1 : 0;
                if (spikespre[i]) tpre[i] = t + 1;
            }

            // Calculate dx (but don't update x yet, matching MATLAB order)
            const dx = new Array(params.N);
            for (let i = 0; i < params.N; i++) {
                dx[i] = (-x[i] + spikespre[i]) / params.tau_ltp;
            }

            // LTD for late pre-synaptic spikes
            for (let i = 0; i < params.N; i++) {
                if (spikespre[i]) {
                    g[i] = g[i] - params.gmax * params.A_ltd * y;
                }
            }

            // Update traces (after LTD, matching MATLAB order)
            for (let i = 0; i < params.N; i++) {
                x[i] = x[i] + dx[i] * dt;
            }
            y = y + dy * dt;

            // Bound weights
            for (let i = 0; i < params.N; i++) {
                g[i] = Math.max(Math.min(g[i], params.gmax), 0);
            }

            idx++;
        }

        // Store initial voltage trace
        if (trial === 1) {
            initialVoltage = [...Vs];
        }

        // Progress callback
        if (onProgress && trial % 100 === 0) {
            const normalizedWeights = g.map(w => w / params.gmax);
            await onProgress(trial, params.ttimes, latencies, normalizedWeights, Vs, initialVoltage);
        }
    }

    return {
        latencies,
        weights: g.map(w => w / params.gmax)
    };
}
