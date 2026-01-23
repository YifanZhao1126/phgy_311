import { Chart, registerables } from 'chart.js';
import './style.css';
Chart.register(...registerables);

class STDPSimulation {
    constructor() {
        this.isRunning = false;
        this.hasSpare = false;
        this.spare = 0;
        this.charts = {};
        this.initializeCharts();
        this.bindEvents();
    }

    initializeCharts() {
        const commonOptions = {
            responsive: true,
            maintainAspectRatio: false,
            animation: false,
            plugins: {
                legend: { display: false }
            },
            scales: {
                x: { title: { display: true, text: 'Input Neuron' }, min: 0, max: 100 },
                y: { title: { display: true, text: 'g/g_max' }, min: -0.05, max: 1.05 }
            }
        };

        // Current simulation chart
        const currentCtx = document.getElementById('currentChart').getContext('2d');
        this.charts.current = new Chart(currentCtx, {
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

        // Condition charts
        ['condition1', 'condition2', 'condition3'].forEach((chartId, index) => {
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
    }

    bindEvents() {
        document.getElementById('runCondition1').addEventListener('click', () => {
            this.runCondition(1);
        });
        document.getElementById('runCondition2').addEventListener('click', () => {
            this.runCondition(2);
        });
        document.getElementById('runCondition3').addEventListener('click', () => {
            this.runCondition(3);
        });
        document.getElementById('stopSimulation').addEventListener('click', () => {
            this.stopSimulation();
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
        document.getElementById('status').textContent = 'Initializing...';

        const params = this.getParameters();
        
        try {
            await this.runSelectivitySimulation(conditionNumber, params);
        } catch (error) {
            console.error('Simulation error:', error);
            document.getElementById('status').textContent = 'Error occurred';
        }

        this.isRunning = false;
        buttons.forEach(id => document.getElementById(id).disabled = false);
    }

    // Use native JavaScript random - Math.random() is fast enough for simulations
    matlabRand() {
        return Math.random();
    }

    generateGaussianRandom() {
        if (this.hasSpare) {
            this.hasSpare = false;
            return this.spare;
        }
        this.hasSpare = true;
        const u = this.matlabRand();
        const v = this.matlabRand();
        const mag = Math.sqrt(-2.0 * Math.log(u));
        this.spare = mag * Math.cos(2.0 * Math.PI * v);
        return mag * Math.sin(2.0 * Math.PI * v);
    }

    // Direct translation of MATLAB simSTDP2 function
    async runSelectivitySimulation(conditionNumber, params) {
        // Reset Gaussian spare for each run
        this.hasSpare = false;

        // Adjust parameters for condition 3 like in MATLAB
        if (conditionNumber === 3) {
            params.gmax = 0.015 * 125; // 1.875
            params.yConst = 15;
        } else {
            params.gmax = 0.015 * 50; // 0.75  
            params.yConst = 2;
        }
        
        const conditionNames = {
            1: 'Uncorrelated ↔ Uncorrelated',
            2: 'Correlated ↔ Uncorrelated', 
            3: 'Correlated ↔ Correlated'
        };
        
        document.getElementById('status').textContent = `Running ${conditionNames[conditionNumber]}...`;
        
        // Set up correlation identifiers based on condition
        const corri = new Array(params.N).fill(0);
        if (conditionNumber === 2) {
            for (let i = 0; i < params.N / 2; i++) {
                corri[i] = 1;
            }
        }
        if (conditionNumber === 3) {
            for (let i = 0; i < params.N / 2; i++) {
                corri[i] = 1;
            }
            for (let i = Math.floor(params.N / 2); i < params.N; i++) {
                corri[i] = 2;
            }
        }
        
        // Direct MATLAB translation starts here
        let V = params.Vrest;
        const dt = 1;
        
        const ratesi = new Array(params.N).fill(0);
        const ratesp = new Array(params.N).fill(0);
        const x = new Array(params.N).fill(0);
        let y = 0;
        const tpre = new Array(params.N).fill(-9999999);
        
        // MATLAB: g = min(max(rand(N,1).*gmax, 0),gmax)
        const g = new Array(params.N);
        for (let i = 0; i < params.N; i++) {
            const randVal = this.matlabRand() * params.gmax;
            g[i] = Math.min(Math.max(randVal, 0), params.gmax);
        }
        
        let dur1 = Math.round(Math.max(1, params.corr_time + this.generateGaussianRandom()));
        let dur2 = Math.round(Math.max(1, params.corr_time + this.generateGaussianRandom()));
        
        // Uncorrelated input
        let xa = new Array(params.N);
        for (let i = 0; i < params.N; i++) {
            xa[i] = this.generateGaussianRandom();
        }
        for (let i = 0; i < params.N; i++) {
            if (corri[i] === 0) {
                ratesi[i] = 10 * (1 + 0.3 * Math.sqrt(2) * xa[i]);
            }
        }
        
        let corr1 = 0;
        if (corri.some(c => c === 1)) {
            corr1 = 1;
            for (let i = 0; i < params.N; i++) {
                xa[i] = this.generateGaussianRandom();
            }
            const ya = this.generateGaussianRandom();
            for (let i = 0; i < params.N; i++) {
                if (corri[i] === 1) {
                    ratesi[i] = 10 * (1 + 0.3 * Math.sqrt(2) * xa[i] + params.yConst * ya);
                }
            }
        }
        
        let corr2 = 0;
        if (corri.some(c => c === 2)) {
            corr2 = 1;
            for (let i = 0; i < params.N; i++) {
                xa[i] = this.generateGaussianRandom();
            }
            const ya = this.generateGaussianRandom();
            for (let i = 0; i < params.N; i++) {
                if (corri[i] === 2) {
                    ratesi[i] = 10 * (1 + 0.3 * Math.sqrt(2) * xa[i] + params.yConst * ya);
                }
            }
        }
        
        for (let i = 0; i < params.N; i++) {
            ratesi[i] = Math.max(ratesi[i], 0);
            if (corri[i] === 0) ratesp[i] = 1 - Math.exp(-ratesi[i] * 0.001);
            if (corri[i] === 1) ratesp[i] = 1 - Math.exp(-ratesi[i] * 0.001);
            if (corri[i] === 2) ratesp[i] = 1 - Math.exp(-ratesi[i] * 0.001);
        }
        
        // Main simulation loop - exact MATLAB translation
        for (let t = 1; t <= params.stime; t++) {
            
            if (t % dur1 === 0) {
                // Generate new interval
                for (let i = 0; i < params.N; i++) {
                    xa[i] = this.generateGaussianRandom();
                }
                for (let i = 0; i < params.N; i++) {
                    if (corri[i] === 0) {
                        ratesi[i] = 10 * (1 + 0.3 * Math.sqrt(2) * xa[i]);
                    }
                }
                dur1 = Math.round(Math.max(1, params.corr_time + this.generateGaussianRandom()));
                
                for (let i = 0; i < params.N; i++) {
                    xa[i] = this.generateGaussianRandom();
                }
                const ya = this.generateGaussianRandom();
                if (corr1) {
                    for (let i = 0; i < params.N; i++) {
                        if (corri[i] === 1) {
                            ratesi[i] = 10 * (1 + 0.3 * xa[i] + params.yConst * ya);
                        }
                    }
                    for (let i = 0; i < params.N; i++) {
                        if (corri[i] === 1) {
                            ratesp[i] = 1 - Math.exp(-ratesi[i] * 0.001);
                        }
                    }
                }
                for (let i = 0; i < params.N; i++) {
                    if (corri[i] === 0) {
                        ratesp[i] = 1 - Math.exp(-ratesi[i] * 0.001);
                    }
                }
                for (let i = 0; i < params.N; i++) {
                    ratesi[i] = Math.max(ratesi[i], 0);
                }
            }
            
            if (corr2 && t % dur2 === 0) {
                for (let i = 0; i < params.N; i++) {
                    xa[i] = this.generateGaussianRandom();
                }
                const ya = this.generateGaussianRandom();
                dur2 = Math.round(Math.max(1, params.corr_time + this.generateGaussianRandom()));
                
                for (let i = 0; i < params.N; i++) {
                    if (corri[i] === 2) {
                        ratesi[i] = 10 * (1 + 0.3 * xa[i] + params.yConst * ya);
                        ratesp[i] = 1 - Math.exp(-ratesi[i] * 0.001);
                    }
                }
            }
            
            // Excitatory input kernel
            let gex = 0;
            for (let i = 0; i < params.N; i++) {
                gex += g[i] * Math.exp(-(t - tpre[i]) / params.tau_ex);
            }
            
            // Integrate-and-fire neuron model
            const dV = (params.Vrest - V + gex * (params.Eex - V)) / params.tau_m;
            V = V + dV * dt;
            
            let spost = 0;
            if (V >= params.Vth) {
                V = -60;
                spost = 1;
                for (let i = 0; i < params.N; i++) {
                    g[i] = g[i] + params.gmax * params.A_ltp * x[i];
                }
            }
            const dy = (-y + spost) / params.tau_ltd;
            
            const spikespre = new Array(params.N);
            for (let i = 0; i < params.N; i++) {
                spikespre[i] = this.matlabRand() <= ratesp[i] ? 1 : 0;
                if (spikespre[i]) tpre[i] = t + 1;
            }
            
            for (let i = 0; i < params.N; i++) {
                const dx = (-x[i] + spikespre[i]) / params.tau_ltp;
                x[i] = x[i] + dx * dt;
            }
            
            // Spike-timing-dependent plasticity
            for (let i = 0; i < params.N; i++) {
                if (spikespre[i]) {
                    g[i] = g[i] - params.gmax * params.A_ltd * y;
                }
            }
            
            y = y + dy * dt;
            
            // Bounds conductances
            for (let i = 0; i < params.N; i++) {
                g[i] = Math.max(Math.min(g[i], params.gmax), 0);
            }
            
            if (t % 50000 === 0) {
                const progress = (t / params.stime) * 100;
                document.getElementById('progressBar').style.width = `${progress}%`;
                document.getElementById('status').textContent = `${conditionNames[conditionNumber]}: ${Math.round(progress)}%`;

                const weightsData = g.map((w, i) => ({
                    x: i,
                    y: w / params.gmax
                }));
                this.charts.current.data.datasets[0].data = weightsData;
                this.charts.current.update('none');

                await new Promise(resolve => setTimeout(resolve, 0));
            }
        }
        
        // Store final results
        const finalWeightsData = g.map((w, i) => ({
            x: i,
            y: w / params.gmax
        }));
        
        this.charts[`condition${conditionNumber}`].data.datasets[0].data = finalWeightsData;
        this.charts[`condition${conditionNumber}`].update();
        
        document.getElementById('progressBar').style.width = '100%';
        document.getElementById('status').textContent = `${conditionNames[conditionNumber]} completed`;
    }

    stopSimulation() {
        this.isRunning = false;
        document.getElementById('status').textContent = 'Simulation stopped';
        const buttons = ['runCondition1', 'runCondition2', 'runCondition3'];
        buttons.forEach(id => document.getElementById(id).disabled = false);
    }
}

class Figure4Simulation {
    constructor() {
        this.isRunning = false;
        this.hasSpare = false;
        this.spare = 0;
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
            plugins: {
                legend: { display: false }
            },
            scales: {
                x: { title: { display: true, text: 'Latency (ms)' }, min: -55, max: 55 },
                y: { title: { display: true, text: 'g/g_max' }, min: 0, max: 1 }
            }
        };

        const voltageOptions = {
            responsive: true,
            maintainAspectRatio: false,
            animation: false,
            plugins: {
                legend: { display: false }
            },
            scales: {
                x: { title: { display: true, text: 'Time (ms)' }, min: -50, max: 50 },
                y: { title: { display: true, text: 'Voltage (mV)' } }
            }
        };

        // Initial weights chart
        const initWeightsCtx = document.getElementById('f4_initialWeights').getContext('2d');
        this.charts.initialWeights = new Chart(initWeightsCtx, {
            type: 'scatter',
            data: {
                datasets: [{
                    data: [],
                    backgroundColor: 'rgba(0, 0, 0, 0.8)',
                    borderColor: 'black',
                    borderWidth: 1,
                    pointRadius: 2
                }]
            },
            options: weightsOptions
        });

        // Initial voltage chart
        const initVoltageCtx = document.getElementById('f4_initialVoltage').getContext('2d');
        this.charts.initialVoltage = new Chart(initVoltageCtx, {
            type: 'scatter',
            data: {
                datasets: [{
                    data: [],
                    borderColor: 'black',
                    borderWidth: 1,
                    fill: false,
                    pointRadius: 0,
                    showLine: true
                }]
            },
            options: voltageOptions
        });

        // Current weights chart
        const currWeightsCtx = document.getElementById('f4_currentWeights').getContext('2d');
        this.charts.currentWeights = new Chart(currWeightsCtx, {
            type: 'scatter',
            data: {
                datasets: [{
                    data: [],
                    backgroundColor: 'rgba(0, 0, 0, 0.8)',
                    borderColor: 'black',
                    borderWidth: 1,
                    pointRadius: 2
                }]
            },
            options: weightsOptions
        });

        // Current voltage chart
        const currVoltageCtx = document.getElementById('f4_currentVoltage').getContext('2d');
        this.charts.currentVoltage = new Chart(currVoltageCtx, {
            type: 'scatter',
            data: {
                datasets: [{
                    data: [],
                    borderColor: 'black',
                    borderWidth: 1,
                    fill: false,
                    pointRadius: 0,
                    showLine: true
                }]
            },
            options: voltageOptions
        });
    }

    bindEvents() {
        document.getElementById('runFigure4').addEventListener('click', () => {
            this.runSimulation();
        });
        document.getElementById('stopFigure4').addEventListener('click', () => {
            this.stopSimulation();
        });
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

    matlabRand() {
        return Math.random();
    }

    generateGaussianRandom() {
        if (this.hasSpare) {
            this.hasSpare = false;
            return this.spare;
        }
        this.hasSpare = true;
        const u = this.matlabRand();
        const v = this.matlabRand();
        const mag = Math.sqrt(-2.0 * Math.log(u));
        this.spare = mag * Math.cos(2.0 * Math.PI * v);
        return mag * Math.sin(2.0 * Math.PI * v);
    }

    async runSimulation() {
        if (this.isRunning) return;

        // Initialize charts if not done yet
        this.initializeCharts();

        this.isRunning = true;
        document.getElementById('runFigure4').disabled = true;
        document.getElementById('f4_status').textContent = 'Initializing...';

        const params = this.getParameters();

        try {
            await this.runLatencySimulation(params);
        } catch (error) {
            console.error('Simulation error:', error);
            document.getElementById('f4_status').textContent = 'Error occurred';
        }

        this.isRunning = false;
        document.getElementById('runFigure4').disabled = false;
    }

    async runLatencySimulation(params) {
        this.hasSpare = false;

        const dt = 1;
        const int_start = -50;
        const int_end = 50;
        const int_length = int_end - int_start + 1;

        // Generate latencies with Gaussian distribution (15ms std dev)
        const latencies = new Array(params.N);
        for (let i = 0; i < params.N; i++) {
            latencies[i] = this.generateGaussianRandom() * params.latency_std;
        }

        // Initialize conductances
        const g = new Array(params.N).fill(0.003);
        const firingp = 1 - Math.exp(-params.burstrate * 0.001);

        // Plot initial weights
        const initialWeightsData = latencies.map((lat, i) => ({
            x: lat,
            y: g[i] / params.gmax
        }));
        this.charts.initialWeights.data.datasets[0].data = initialWeightsData;
        this.charts.initialWeights.update('none');

        // Main simulation loop
        for (let xx = 1; xx <= params.ttimes; xx++) {
            if (!this.isRunning) break;

            const tpre = new Array(params.N).fill(-9999999);
            let V = params.Vrest;
            const x = new Array(params.N).fill(0);
            let y = 0;
            const Vs = new Array(int_length).fill(0);

            let idx = 0;
            for (let t = int_start; t <= int_end; t++) {
                // Excitatory input kernel
                let gex = 0;
                for (let i = 0; i < params.N; i++) {
                    gex += g[i] * Math.exp(-(t - tpre[i]) / params.tau_ex);
                }

                // Integrate-and-fire neuron model
                const dV = (params.Vrest - V + gex * (params.Eex - V)) / params.tau_m;
                V = V + dV * dt;
                Vs[idx] = V;

                let spost = 0;
                if (V >= params.Vth) {
                    V = -60;
                    Vs[idx] = 0;
                    spost = 1;
                    // LTP
                    for (let i = 0; i < params.N; i++) {
                        g[i] = g[i] + params.gmax * params.A_ltp * x[i];
                    }
                }
                const dy = (-y + spost) / params.tau_ltd;

                // Pre-synaptic spikes: only during burst window
                const spikespre = new Array(params.N);
                for (let i = 0; i < params.N; i++) {
                    const inBurstWindow = t >= latencies[i] && t < (latencies[i] + params.burstdur);
                    spikespre[i] = (this.matlabRand() <= firingp && inBurstWindow) ? 1 : 0;
                    if (spikespre[i]) tpre[i] = t + 1;
                }

                // Update pre-synaptic traces
                for (let i = 0; i < params.N; i++) {
                    const dx = (-x[i] + spikespre[i]) / params.tau_ltp;
                    x[i] = x[i] + dx * dt;
                }

                // LTD
                for (let i = 0; i < params.N; i++) {
                    if (spikespre[i]) {
                        g[i] = g[i] - params.gmax * params.A_ltd * y;
                    }
                }

                y = y + dy * dt;

                // Bound conductances
                for (let i = 0; i < params.N; i++) {
                    g[i] = Math.max(Math.min(g[i], params.gmax), 0);
                }

                idx++;
            }

            // Plot initial voltage trace after first iteration
            if (xx === 1) {
                const voltageData = Vs.map((v, i) => ({
                    x: int_start + i,
                    y: v
                }));
                this.charts.initialVoltage.data.datasets[0].data = voltageData;
                this.charts.initialVoltage.update('none');
            }

            // Update progress and plots (every 100 iterations for performance)
            if (xx % 100 === 0) {
                const progress = (xx / params.ttimes) * 100;
                document.getElementById('f4_progressBar').style.width = `${progress}%`;
                document.getElementById('f4_status').textContent = `Iteration ${xx}/${params.ttimes}`;

                const weightsData = latencies.map((lat, i) => ({
                    x: lat,
                    y: g[i] / params.gmax
                }));
                this.charts.currentWeights.data.datasets[0].data = weightsData;
                this.charts.currentWeights.update('none');

                const voltageData = Vs.map((v, i) => ({
                    x: int_start + i,
                    y: v
                }));
                this.charts.currentVoltage.data.datasets[0].data = voltageData;
                this.charts.currentVoltage.update('none');

                await new Promise(resolve => setTimeout(resolve, 1));
            }
        }

        // Final update
        const finalWeightsData = latencies.map((lat, i) => ({
            x: lat,
            y: g[i] / params.gmax
        }));
        this.charts.currentWeights.data.datasets[0].data = finalWeightsData;
        this.charts.currentWeights.update();

        document.getElementById('f4_progressBar').style.width = '100%';
        document.getElementById('f4_status').textContent = 'Simulation completed';
    }

    stopSimulation() {
        this.isRunning = false;
        document.getElementById('f4_status').textContent = 'Simulation stopped';
        document.getElementById('runFigure4').disabled = false;
    }
}

// Tab switching functionality
let figure4Sim = null;

function setupTabs() {
    const tabButtons = document.querySelectorAll('.tab-button');
    const tabContents = document.querySelectorAll('.tab-content');

    tabButtons.forEach(button => {
        button.addEventListener('click', () => {
            const tabId = button.getAttribute('data-tab');

            // Update active button
            tabButtons.forEach(btn => btn.classList.remove('active'));
            button.classList.add('active');

            // Update active content
            tabContents.forEach(content => content.classList.remove('active'));
            document.getElementById(tabId).classList.add('active');

            // Initialize Figure 4 charts when tab is first shown
            if (tabId === 'figure4' && figure4Sim) {
                figure4Sim.initializeCharts();
            }
        });
    });
}

// Initialize the simulations when the page loads
document.addEventListener('DOMContentLoaded', () => {
    new STDPSimulation();
    figure4Sim = new Figure4Simulation();
    setupTabs();
});