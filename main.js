import { Chart, registerables } from 'chart.js';
import './style.css';
import { runSelectivitySimulation, runLatencySimulation } from './simulation.js';

Chart.register(...registerables);

// =============================================================================
// FIGURE 7: Input Selectivity Simulation UI
// =============================================================================

class STDPSimulation {
    constructor() {
        this.isRunning = false;
        this.charts = {};
        this.initializeCharts();
        this.bindEvents();
    }

    initializeCharts() {
        const commonOptions = {
            responsive: true,
            maintainAspectRatio: true,
            aspectRatio: 1.5,
            animation: false,
            devicePixelRatio: window.devicePixelRatio || 1,
            plugins: { legend: { display: false } },
            scales: {
                x: { title: { display: true, text: 'Input Neuron' }, min: 0, max: 100 },
                y: { title: { display: true, text: 'g/g_max' }, min: -0.05, max: 1.05 }
            }
        };

        // Condition charts
        ['condition1', 'condition2', 'condition3'].forEach((chartId) => {
            const ctx = document.getElementById(`${chartId}Chart`).getContext('2d');
            this.charts[chartId] = new Chart(ctx, {
                type: 'scatter',
                data: {
                    datasets: [{
                        data: [],
                        backgroundColor: 'rgba(0, 0, 0, 0.8)',
                        borderColor: 'black',
                        borderWidth: 1,
                        pointRadius: 3
                    }]
                },
                options: commonOptions
            });
        });

        // STDP Learning Rule curve
        const stdpCtx = document.getElementById('stdpCurveChart').getContext('2d');
        this.charts.stdpCurve = new Chart(stdpCtx, {
            type: 'scatter',
            data: {
                datasets: [{
                    data: [],
                    borderColor: 'black',
                    borderWidth: 2,
                    pointRadius: 0,
                    showLine: true,
                    fill: false
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: true,
                aspectRatio: 1.5,
                animation: false,
                devicePixelRatio: window.devicePixelRatio || 1,
                plugins: { legend: { display: false } },
                scales: {
                    x: {
                        title: { display: true, text: 'Δt (ms)' },
                        min: -100,
                        max: 100
                    },
                    y: {
                        title: { display: true, text: 'Δw' },
                        min: -0.006,
                        max: 0.006
                    }
                }
            }
        });

        // Draw initial STDP curve
        this.updateSTDPCurve();
    }

    updateSTDPCurve() {
        const tau_ltp = parseFloat(document.getElementById('tau_ltp').value);
        const tau_ltd = parseFloat(document.getElementById('tau_ltp').value);
        const A_ltp = parseFloat(document.getElementById('A_ltp').value);
        const A_ltd = A_ltp * parseFloat(document.getElementById('B').value);

        const data = [];

        // LTD side (Δt < 0: post before pre)
        for (let dt = -100; dt < 0; dt += 1) {
            const dw = -A_ltd * Math.exp(dt / tau_ltd);
            data.push({ x: dt, y: dw });
        }

        // Zero point
        data.push({ x: 0, y: 0 });

        // LTP side (Δt > 0: pre before post)
        for (let dt = 1; dt <= 100; dt += 1) {
            const dw = A_ltp * Math.exp(-dt / tau_ltp);
            data.push({ x: dt, y: dw });
        }

        this.charts.stdpCurve.data.datasets[0].data = data;
        this.charts.stdpCurve.update();
    }

    bindEvents() {
        document.getElementById('runCondition1').addEventListener('click', () => this.runCondition(1));
        document.getElementById('runCondition2').addEventListener('click', () => this.runCondition(2));
        document.getElementById('runCondition3').addEventListener('click', () => this.runCondition(3));
        document.getElementById('stopSimulation').addEventListener('click', () => this.stopSimulation());

        // Update STDP curve when parameters change
        ['tau_ltp', 'A_ltp', 'B'].forEach(id => {
            document.getElementById(id).addEventListener('change', () => this.updateSTDPCurve());
        });
    }

    getParameters() {
        return {
            N: parseInt(document.getElementById('N').value),
            tau_ltp: parseFloat(document.getElementById('tau_ltp').value),
            tau_ltd: parseFloat(document.getElementById('tau_ltp').value),
            A_ltp: parseFloat(document.getElementById('A_ltp').value),
            A_ltd: parseFloat(document.getElementById('A_ltp').value) * parseFloat(document.getElementById('B').value),
            gmax: parseFloat(document.getElementById('gmax').value),
            Vrest: parseFloat(document.getElementById('Vrest').value),
            Vth: parseFloat(document.getElementById('Vth').value),
            tau_m: parseFloat(document.getElementById('tau_m').value),
            stime: parseInt(document.getElementById('stime').value),
            tau_ex: parseFloat(document.getElementById('tau_ex').value),
            Eex: parseFloat(document.getElementById('Eex').value),
            corr_time: parseFloat(document.getElementById('corr_time').value),
            yConst: parseFloat(document.getElementById('yConst').value)
        };
    }

    async runCondition(conditionNumber) {
        if (this.isRunning) return;

        this.isRunning = true;
        const buttons = ['runCondition1', 'runCondition2', 'runCondition3'];
        buttons.forEach(id => document.getElementById(id).disabled = true);

        const conditionNames = {
            1: 'Uncorrelated ↔ Uncorrelated',
            2: 'Correlated ↔ Uncorrelated',
            3: 'Correlated ↔ Correlated'
        };

        document.getElementById('status').textContent = `Running ${conditionNames[conditionNumber]}...`;

        const params = this.getParameters();

        try {
            // Progress callback for UI updates - update the condition chart directly
            const onProgress = async (t, total, weights) => {
                if (!this.isRunning) throw new Error('Stopped');

                const progress = (t / total) * 100;
                document.getElementById('progressBar').style.width = `${progress}%`;
                document.getElementById('status').textContent = `${conditionNames[conditionNumber]}: ${Math.round(progress)}%`;

                const weightsData = weights.map((w, i) => ({ x: i, y: w }));
                this.charts[`condition${conditionNumber}`].data.datasets[0].data = weightsData;
                this.charts[`condition${conditionNumber}`].update('none');

                await new Promise(resolve => setTimeout(resolve, 5)); // Small delay like MATLAB's pause(0.0005)
            };

            const finalWeights = await runSelectivitySimulation(conditionNumber, params, onProgress);

            // Store final results
            const finalWeightsData = finalWeights.map((w, i) => ({ x: i, y: w }));
            this.charts[`condition${conditionNumber}`].data.datasets[0].data = finalWeightsData;
            this.charts[`condition${conditionNumber}`].update();

            document.getElementById('progressBar').style.width = '100%';
            document.getElementById('status').textContent = `${conditionNames[conditionNumber]} completed`;

        } catch (error) {
            if (error.message !== 'Stopped') {
                console.error('Simulation error:', error);
                document.getElementById('status').textContent = 'Error occurred';
            }
        }

        this.isRunning = false;
        buttons.forEach(id => document.getElementById(id).disabled = false);
    }

    stopSimulation() {
        this.isRunning = false;
        document.getElementById('status').textContent = 'Simulation stopped';
        const buttons = ['runCondition1', 'runCondition2', 'runCondition3'];
        buttons.forEach(id => document.getElementById(id).disabled = false);
    }
}

// =============================================================================
// FIGURE 4: Latency Simulation UI
// =============================================================================

class Figure4Simulation {
    constructor() {
        this.isRunning = false;
        this.charts = {};
        this.chartsInitialized = false;
        this.bindEvents();
    }

    initializeCharts() {
        if (this.chartsInitialized) return;
        this.chartsInitialized = true;

        const weightsOptions = {
            responsive: true,
            maintainAspectRatio: false,
            animation: false,
            plugins: { legend: { display: false } },
            scales: {
                x: { title: { display: true, text: 'Latency (ms)' }, min: -55, max: 55 },
                y: { title: { display: true, text: 'g/g_max' }, min: 0, max: 1 }
            }
        };

        const voltageOptions = {
            responsive: true,
            maintainAspectRatio: false,
            animation: false,
            plugins: { legend: { display: false } },
            scales: {
                x: { title: { display: true, text: 'Time (ms)' }, min: -50, max: 50 },
                y: { title: { display: true, text: 'Voltage (mV)' } }
            }
        };

        // Initial weights chart
        const initWeightsCtx = document.getElementById('f4_initialWeights').getContext('2d');
        this.charts.initialWeights = new Chart(initWeightsCtx, {
            type: 'scatter',
            data: { datasets: [{ data: [], backgroundColor: 'rgba(0, 0, 0, 0.8)', borderColor: 'black', borderWidth: 1, pointRadius: 2 }] },
            options: weightsOptions
        });

        // Initial voltage chart
        const initVoltageCtx = document.getElementById('f4_initialVoltage').getContext('2d');
        this.charts.initialVoltage = new Chart(initVoltageCtx, {
            type: 'scatter',
            data: { datasets: [{ data: [], borderColor: 'black', borderWidth: 1, fill: false, pointRadius: 0, showLine: true }] },
            options: voltageOptions
        });

        // Current weights chart
        const currWeightsCtx = document.getElementById('f4_currentWeights').getContext('2d');
        this.charts.currentWeights = new Chart(currWeightsCtx, {
            type: 'scatter',
            data: { datasets: [{ data: [], backgroundColor: 'rgba(0, 0, 0, 0.8)', borderColor: 'black', borderWidth: 1, pointRadius: 2 }] },
            options: weightsOptions
        });

        // Current voltage chart
        const currVoltageCtx = document.getElementById('f4_currentVoltage').getContext('2d');
        this.charts.currentVoltage = new Chart(currVoltageCtx, {
            type: 'scatter',
            data: { datasets: [{ data: [], borderColor: 'black', borderWidth: 1, fill: false, pointRadius: 0, showLine: true }] },
            options: voltageOptions
        });
    }

    bindEvents() {
        document.getElementById('runFigure4').addEventListener('click', () => this.runSimulation());
        document.getElementById('stopFigure4').addEventListener('click', () => this.stopSimulation());
    }

    getParameters() {
        return {
            N: parseInt(document.getElementById('f4_N').value),
            tau_ltp: parseFloat(document.getElementById('f4_tau_ltp').value),
            tau_ltd: parseFloat(document.getElementById('f4_tau_ltp').value),
            A_ltp: parseFloat(document.getElementById('f4_A_ltp').value),
            A_ltd: parseFloat(document.getElementById('f4_A_ltp').value) * parseFloat(document.getElementById('f4_B').value),
            gmax: parseFloat(document.getElementById('f4_gmax').value),
            Vrest: parseFloat(document.getElementById('f4_Vrest').value),
            Vth: parseFloat(document.getElementById('f4_Vth').value),
            tau_m: parseFloat(document.getElementById('f4_tau_m').value),
            tau_ex: parseFloat(document.getElementById('f4_tau_ex').value),
            Eex: parseFloat(document.getElementById('f4_Eex').value),
            ttimes: parseInt(document.getElementById('f4_ttimes').value),
            latency_std: parseFloat(document.getElementById('f4_latency_std').value),
            burstdur: parseFloat(document.getElementById('f4_burstdur').value),
            burstrate: parseFloat(document.getElementById('f4_burstrate').value)
        };
    }

    async runSimulation() {
        if (this.isRunning) return;

        this.initializeCharts();
        this.isRunning = true;
        document.getElementById('runFigure4').disabled = true;
        document.getElementById('f4_status').textContent = 'Initializing...';

        const params = this.getParameters();

        try {
            // Progress callback
            const onProgress = async (trial, total, latencies, weights, voltage, initialVoltage) => {
                if (!this.isRunning) throw new Error('Stopped');

                const progress = (trial / total) * 100;
                document.getElementById('f4_progressBar').style.width = `${progress}%`;
                document.getElementById('f4_status').textContent = `Iteration ${trial}/${total}`;

                // Update initial weights/voltage on first callback
                if (trial === 100) {
                    const initialWeightsData = latencies.map((lat, i) => ({ x: lat, y: 0.003 / params.gmax }));
                    this.charts.initialWeights.data.datasets[0].data = initialWeightsData;
                    this.charts.initialWeights.update('none');

                    if (initialVoltage) {
                        const voltageData = initialVoltage.map((v, i) => ({ x: -50 + i, y: v }));
                        this.charts.initialVoltage.data.datasets[0].data = voltageData;
                        this.charts.initialVoltage.update('none');
                    }
                }

                // Update current weights/voltage
                const weightsData = latencies.map((lat, i) => ({ x: lat, y: weights[i] }));
                this.charts.currentWeights.data.datasets[0].data = weightsData;
                this.charts.currentWeights.update('none');

                const voltageData = voltage.map((v, i) => ({ x: -50 + i, y: v }));
                this.charts.currentVoltage.data.datasets[0].data = voltageData;
                this.charts.currentVoltage.update('none');

                await new Promise(resolve => setTimeout(resolve, 1));
            };

            const { latencies, weights } = await runLatencySimulation(params, onProgress);

            // Final update
            const finalWeightsData = latencies.map((lat, i) => ({ x: lat, y: weights[i] }));
            this.charts.currentWeights.data.datasets[0].data = finalWeightsData;
            this.charts.currentWeights.update();

            document.getElementById('f4_progressBar').style.width = '100%';
            document.getElementById('f4_status').textContent = 'Simulation completed';

        } catch (error) {
            if (error.message !== 'Stopped') {
                console.error('Simulation error:', error);
                document.getElementById('f4_status').textContent = 'Error occurred';
            }
        }

        this.isRunning = false;
        document.getElementById('runFigure4').disabled = false;
    }

    stopSimulation() {
        this.isRunning = false;
        document.getElementById('f4_status').textContent = 'Simulation stopped';
        document.getElementById('runFigure4').disabled = false;
    }
}

// =============================================================================
// Tab Switching
// =============================================================================

let figure4Sim = null;

function setupTabs() {
    const tabButtons = document.querySelectorAll('.tab-button');
    const tabContents = document.querySelectorAll('.tab-content');

    tabButtons.forEach(button => {
        button.addEventListener('click', () => {
            const tabId = button.getAttribute('data-tab');

            tabButtons.forEach(btn => btn.classList.remove('active'));
            button.classList.add('active');

            tabContents.forEach(content => content.classList.remove('active'));
            document.getElementById(tabId).classList.add('active');

            if (tabId === 'figure4' && figure4Sim) {
                figure4Sim.initializeCharts();
            }
        });
    });
}

// Initialize on page load
document.addEventListener('DOMContentLoaded', () => {
    new STDPSimulation();
    figure4Sim = new Figure4Simulation();
    setupTabs();
});
