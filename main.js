// Plotly is loaded from CDN in index.html (same as example)
import { runSelectivitySimulation, runLatencySimulation } from './simulation.js';
import './style.css';

// =============================================================================
// Plotly Configuration (same as example)
// =============================================================================

var config = {responsive: true};

// =============================================================================
// Figure 1 Layouts
// =============================================================================

var layout1 = {
    title: '<b>Uncorrelated \u2194 Uncorrelated</b>',
    xaxis: {
        title: '<b>Input Neuron</b>',
        range: [0, 100],
        zeroline: false
    },
    yaxis: {
        title: '<b>Synaptic Strength g/g<sub>max</sub></b>',
        range: [0, 1],
        zeroline: false
    },
    autosize: true,
    paper_bgcolor: '#c7c7c7',
    showlegend: false
};

var layout2 = {
    title: '<b>Correlated \u2194 Uncorrelated</b>',
    xaxis: {
        title: '<b>Input Neuron</b>',
        range: [0, 100],
        zeroline: false
    },
    yaxis: {
        title: '<b>Synaptic Strength g/g<sub>max</sub></b>',
        range: [0, 1],
        zeroline: false
    },
    autosize: true,
    paper_bgcolor: '#c7c7c7',
    showlegend: false
};

var layout3 = {
    title: '<b>Correlated \u2194 Correlated</b>',
    xaxis: {
        title: '<b>Input Neuron</b>',
        range: [0, 100],
        zeroline: false
    },
    yaxis: {
        title: '<b>Synaptic Strength g/g<sub>max</sub></b>',
        range: [0, 1],
        zeroline: false
    },
    autosize: true,
    paper_bgcolor: '#c7c7c7',
    showlegend: false
};

var layoutHist = {
    title: '<b>Weight Distribution</b>',
    xaxis: {
        title: '<b>g/g<sub>max</sub></b>',
        zeroline: false
    },
    yaxis: {
        title: '<b>Count</b>',
        zeroline: false
    },
    autosize: true,
    paper_bgcolor: '#c7c7c7',
    showlegend: false
};

// =============================================================================
// Figure 4 Layouts
// =============================================================================

var layoutF4InitWeights = {
    title: '<b>Initial Synaptic Weights</b>',
    xaxis: {
        title: '<b>Latency (ms)</b>',
        range: [-55, 55],
        zeroline: false
    },
    yaxis: {
        title: '<b>g/g<sub>max</sub></b>',
        range: [0, 1],
        zeroline: false
    },
    autosize: true,
    paper_bgcolor: '#c7c7c7',
    showlegend: false
};

var layoutF4InitVoltage = {
    title: '<b>Initial Voltage Trace</b>',
    xaxis: {
        title: '<b>Time (ms)</b>',
        range: [-50, 50],
        zeroline: false
    },
    yaxis: {
        title: '<b>Voltage (mV)</b>',
        range: [-80, 0],
        zeroline: false
    },
    autosize: true,
    paper_bgcolor: '#c7c7c7',
    showlegend: false
};

var layoutF4CurrWeights = {
    title: '<b>Synaptic Weights</b>',
    xaxis: {
        title: '<b>Latency (ms)</b>',
        range: [-55, 55],
        zeroline: false
    },
    yaxis: {
        title: '<b>g/g<sub>max</sub></b>',
        range: [0, 1],
        zeroline: false
    },
    autosize: true,
    paper_bgcolor: '#c7c7c7',
    showlegend: false
};

var layoutF4CurrVoltage = {
    title: '<b>Voltage Trace</b>',
    xaxis: {
        title: '<b>Time (ms)</b>',
        range: [-50, 50],
        zeroline: false
    },
    yaxis: {
        title: '<b>Voltage (mV)</b>',
        range: [-80, 0],
        zeroline: false
    },
    autosize: true,
    paper_bgcolor: '#c7c7c7',
    showlegend: false
};

// =============================================================================
// FIGURE 1: Input Selectivity Simulation
// =============================================================================

var isRunning = false;

function getParameters() {
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

function initializePlots() {
    // Initialize empty scatter plots for Figure 1
    var emptyTrace = {
        x: [],
        y: [],
        mode: 'markers',
        type: 'scatter',
        marker: {
            color: 'rgb(0, 0, 0)',
            size: 6
        },
        cliponaxis: false
    };

    Plotly.newPlot('condition1Chart', [emptyTrace], layout1, config);
    Plotly.newPlot('condition2Chart', [emptyTrace], layout2, config);
    Plotly.newPlot('condition3Chart', [emptyTrace], layout3, config);

    // Initialize empty histogram
    var emptyHist = {
        x: [],
        type: 'histogram',
        marker: {
            color: 'rgba(0, 0, 0, 0.7)'
        },
        xbins: {
            start: 0,
            end: 1,
            size: 0.1
        }
    };
    Plotly.newPlot('histogramChart', [emptyHist], layoutHist, config);
}

function updateConditionPlot(conditionNumber, weights) {
    var chartId = 'condition' + conditionNumber + 'Chart';
    var x = [];
    var y = [];
    for (var i = 0; i < weights.length; i++) {
        x.push(i);
        y.push(weights[i]);
    }

    // Use restyle for smoother updates without flashing
    Plotly.restyle(chartId, {
        x: [x],
        y: [y]
    }, [0]);
}

function updateHistogram(weights) {
    // Use restyle for smoother updates
    Plotly.restyle('histogramChart', {
        x: [weights]
    }, [0]);
}

async function runCondition(conditionNumber) {
    if (isRunning) return;

    isRunning = true;
    var buttons = ['runCondition1', 'runCondition2', 'runCondition3'];
    buttons.forEach(function(id) { document.getElementById(id).disabled = true; });

    var conditionNames = {
        1: 'Uncorrelated \u2194 Uncorrelated',
        2: 'Correlated \u2194 Uncorrelated',
        3: 'Correlated \u2194 Correlated'
    };

    document.getElementById('status').textContent = 'Running ' + conditionNames[conditionNumber] + '...';

    var params = getParameters();

    try {
        var onProgress = async function(t, total, weights) {
            if (!isRunning) throw new Error('Stopped');

            var progress = (t / total) * 100;
            document.getElementById('progressBar').style.width = progress + '%';
            document.getElementById('status').textContent = conditionNames[conditionNumber] + ': ' + Math.round(progress) + '%';

            updateConditionPlot(conditionNumber, weights);

            await new Promise(function(resolve) { setTimeout(resolve, 5); });
        };

        var finalWeights = await runSelectivitySimulation(conditionNumber, params, onProgress);

        updateConditionPlot(conditionNumber, finalWeights);
        updateHistogram(finalWeights);

        document.getElementById('progressBar').style.width = '100%';
        document.getElementById('status').textContent = conditionNames[conditionNumber] + ' completed';

    } catch (error) {
        if (error.message !== 'Stopped') {
            console.error('Simulation error:', error);
            document.getElementById('status').textContent = 'Error occurred';
        }
    }

    isRunning = false;
    buttons.forEach(function(id) { document.getElementById(id).disabled = false; });
}

function stopSimulation() {
    isRunning = false;
    document.getElementById('status').textContent = 'Simulation stopped';
    var buttons = ['runCondition1', 'runCondition2', 'runCondition3'];
    buttons.forEach(function(id) { document.getElementById(id).disabled = false; });
}

// =============================================================================
// FIGURE 4: Latency Simulation
// =============================================================================

var f4_isRunning = false;
var f4_chartsInitialized = false;

function getF4Parameters() {
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

function initializeF4Plots() {
    if (f4_chartsInitialized) return;
    f4_chartsInitialized = true;

    var emptyScatter = {
        x: [],
        y: [],
        mode: 'markers',
        type: 'scatter',
        marker: {
            color: 'rgb(0, 0, 0)',
            size: 4
        },
        cliponaxis: false
    };

    var emptyLine = {
        x: [],
        y: [],
        mode: 'lines',
        type: 'scatter',
        line: {
            color: 'rgb(0, 0, 0)',
            width: 2
        }
    };

    Plotly.newPlot('f4_initialWeights', [emptyScatter], layoutF4InitWeights, config);
    Plotly.newPlot('f4_initialVoltage', [emptyLine], layoutF4InitVoltage, config);
    Plotly.newPlot('f4_currentWeights', [emptyScatter], layoutF4CurrWeights, config);
    Plotly.newPlot('f4_currentVoltage', [emptyLine], layoutF4CurrVoltage, config);
}

async function runFigure4() {
    if (f4_isRunning) return;

    initializeF4Plots();
    f4_isRunning = true;
    document.getElementById('runFigure4').disabled = true;
    document.getElementById('f4_status').textContent = 'Initializing...';

    var params = getF4Parameters();

    try {
        var onProgress = async function(trial, total, latencies, weights, voltage, initialVoltage) {
            if (!f4_isRunning) throw new Error('Stopped');

            var progress = (trial / total) * 100;
            document.getElementById('f4_progressBar').style.width = progress + '%';
            document.getElementById('f4_status').textContent = 'Iteration ' + trial + '/' + total;

            // Update initial weights/voltage on first callback
            if (trial === 100) {
                Plotly.restyle('f4_initialWeights', {
                    x: [latencies],
                    y: [new Array(latencies.length).fill(0.003 / params.gmax)]
                }, [0]);

                if (initialVoltage) {
                    var initTimeAxis = [];
                    for (var i = 0; i < initialVoltage.length; i++) {
                        initTimeAxis.push(-50 + i);
                    }
                    Plotly.restyle('f4_initialVoltage', {
                        x: [initTimeAxis],
                        y: [initialVoltage]
                    }, [0]);
                }
            }

            // Update current weights using restyle (no flashing)
            Plotly.restyle('f4_currentWeights', {
                x: [latencies],
                y: [weights]
            }, [0]);

            // Update current voltage using restyle (no flashing)
            var timeAxis = [];
            for (var i = 0; i < voltage.length; i++) {
                timeAxis.push(-50 + i);
            }
            Plotly.restyle('f4_currentVoltage', {
                x: [timeAxis],
                y: [voltage]
            }, [0]);

            await new Promise(function(resolve) { setTimeout(resolve, 1); });
        };

        var result = await runLatencySimulation(params, onProgress);

        // Final update using restyle
        Plotly.restyle('f4_currentWeights', {
            x: [result.latencies],
            y: [result.weights]
        }, [0]);

        document.getElementById('f4_progressBar').style.width = '100%';
        document.getElementById('f4_status').textContent = 'Simulation completed';

    } catch (error) {
        if (error.message !== 'Stopped') {
            console.error('Simulation error:', error);
            document.getElementById('f4_status').textContent = 'Error occurred';
        }
    }

    f4_isRunning = false;
    document.getElementById('runFigure4').disabled = false;
}

function stopFigure4() {
    f4_isRunning = false;
    document.getElementById('f4_status').textContent = 'Simulation stopped';
    document.getElementById('runFigure4').disabled = false;
}

// =============================================================================
// Tab Switching
// =============================================================================

function setupTabs() {
    var tabButtons = document.querySelectorAll('.tab-button');
    var tabContents = document.querySelectorAll('.tab-content');

    tabButtons.forEach(function(button) {
        button.addEventListener('click', function() {
            var tabId = button.getAttribute('data-tab');

            tabButtons.forEach(function(btn) { btn.classList.remove('active'); });
            button.classList.add('active');

            tabContents.forEach(function(content) { content.classList.remove('active'); });
            document.getElementById(tabId).classList.add('active');

            if (tabId === 'figure4') {
                initializeF4Plots();
            }
        });
    });
}

// =============================================================================
// Initialize on page load
// =============================================================================

document.addEventListener('DOMContentLoaded', function() {
    initializePlots();
    setupTabs();

    // Bind event listeners
    document.getElementById('runCondition1').addEventListener('click', function() { runCondition(1); });
    document.getElementById('runCondition2').addEventListener('click', function() { runCondition(2); });
    document.getElementById('runCondition3').addEventListener('click', function() { runCondition(3); });
    document.getElementById('stopSimulation').addEventListener('click', stopSimulation);

    document.getElementById('runFigure4').addEventListener('click', runFigure4);
    document.getElementById('stopFigure4').addEventListener('click', stopFigure4);
});
