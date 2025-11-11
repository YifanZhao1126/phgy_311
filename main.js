import { Chart, registerables } from 'chart.js';
Chart.register(...registerables);

class STDPSimulation {
    constructor() {
        this.isRunning = false;
        this.hasSpare = false;
        this.spare = 0;
        this.charts = {};
        // MATLAB's default random seed behavior
        this.matlabSeed = 0; // MATLAB default 
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

    // MATLAB-compatible random number generator
    // Using same algorithm as MATLAB's default rand() function
    matlabRand() {
        // Simple linear congruential generator matching MATLAB's behavior
        this.matlabSeed = (this.matlabSeed * 134775813 + 1) % Math.pow(2, 32);
        return this.matlabSeed / Math.pow(2, 32);
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
        // Reset to MATLAB's default seed behavior
        this.matlabSeed = 0;
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
        // DEBUG: Let's check what MATLAB actually does
        console.log('Before setting corri, conditionNumber:', conditionNumber);
        const corri = new Array(params.N).fill(0);
        if (conditionNumber === 2) {
            // MATLAB: corri(1:end/2) = 1; (first half)
            for (let i = 0; i < params.N / 2; i++) {
                corri[i] = 1;
            }
            console.log('Condition 2: First 50 neurons set to corri=1');
        }
        if (conditionNumber === 3) {
            // MATLAB: corri(1:end/2) = 1; corri(end/2+1:end) = 2;
            for (let i = 0; i < params.N / 2; i++) {
                corri[i] = 1;
            }
            for (let i = Math.floor(params.N / 2); i < params.N; i++) {
                corri[i] = 2;
            }
            console.log('Condition 3: First 50 = corri=1, Last 50 = corri=2');
        }
        console.log('Corri setup:', corri.slice(0, 10), '...', corri.slice(-10));
        
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
            
            if (t % 10000 === 0) {
                const progress = (t / params.stime) * 100;
                document.getElementById('progressBar').style.width = `${progress}%`;
                document.getElementById('status').textContent = `${conditionNames[conditionNumber]}: ${Math.round(progress)}%`;
                
                const weightsData = g.map((w, i) => ({
                    x: i,
                    y: w / params.gmax
                }));
                this.charts.current.data.datasets[0].data = weightsData;
                this.charts.current.update('none');
                
                await new Promise(resolve => setTimeout(resolve, 1));
            }
        }
        
        // Store final results
        const finalWeightsData = g.map((w, i) => ({
            x: i,
            y: w / params.gmax
        }));
        
        // DEBUG: Check final weights pattern
        console.log(`Final weights for condition ${conditionNumber}:`);
        console.log('First 10 weights:', g.slice(0, 10).map(w => (w/params.gmax).toFixed(3)));
        console.log('Last 10 weights:', g.slice(-10).map(w => (w/params.gmax).toFixed(3)));
        console.log('Average first half:', (g.slice(0, 50).reduce((a,b) => a+b, 0)/50/params.gmax).toFixed(3));
        console.log('Average last half:', (g.slice(50, 100).reduce((a,b) => a+b, 0)/50/params.gmax).toFixed(3));
        
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

// Initialize the simulation when the page loads
document.addEventListener('DOMContentLoaded', () => {
    new STDPSimulation();
});