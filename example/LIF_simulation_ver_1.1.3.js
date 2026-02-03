var config = {responsive: true};

var layout = {
  title: '<b>Membrane Potential</b>',
  xaxis: {
    title: '<b>Time (s)</b>'
  },
  yaxis: {
    title: '<b>V<sub>m</sub> (mV)</b>'
  }, 
  autosize: true, 
  paper_bgcolor: '#c7c7c7', 
  showlegend: false
}; 

var layout2 = {
  title: '<b>Membrane Current</b>',
  xaxis: {
    title: '<b>Time (s)</b>'
  },
  yaxis: {
    title: '<b>I<sub>m</sub> (pA)</b>'
  }, 
  autosize: true, 
  paper_bgcolor: '#c7c7c7', 
  showlegend: false
}; 

var layout3 = {
  title: '<b>Short-Term Depression</b>',
  xaxis: {
    title: '<b>Time (s)</b>'
  },
  yaxis: {
    title: '<b>D, fraction of RRP filled</b>'
  }, 
  autosize: true, 
  paper_bgcolor: '#c7c7c7', 
  showlegend: false
};

var layout4 = {
  title: '<b>Short-Term Facilitation</b>',
  xaxis: {
    title: '<b>Time (s)</b>'
  },
  yaxis: {
    title: '<b>F, calcium concentration</b>'
  }, 
  autosize: true, 
  paper_bgcolor: '#c7c7c7', 
  showlegend: false
};

var makePulseTrain = function(tStart,freq,nPulse,wDur) { 
	var pulseTrain1 = []; 
	for (var j = 0; j < nPulse; j++) {
		if ((1/freq*j+tStart*1e-3) < (wDur*1e-3)) {
			pulseTrain1.push(1/freq*j+tStart*1e-3); 
		}
	}
	return pulseTrain1; 
}

var PulseNow = function(t,deltaT,PulseTrain) { 
	var n = PulseTrain.length; 
	if (n > 0) { 
		for (var i = 0; i < n; i++) { 
			if (Math.abs(t-PulseTrain[i]) < (deltaT/2)) { 
				return 1; 
			}
		} 
	}
	return 0; 
}	

var currFlag = 0; 
var synFlag = 0; 
var extTrainFlag = 0; 

var deltaT = 0.0001;	//s 
var tStart = 30;	//ms
var tau1 = 0.1; 	//ms
var tau2 = 5.26; 	//ms
var tauM = 25; 	//ms
var E_leak = -70; 	//mV
var freq = 50; 	//Hz
var R_input = 200; 	//MOhm
var g_syn = 1; 	//nS
var E_syn = 0; 	//mV 
var C_M = tauM*1e-3/(R_input*1e6)*1e12; 	//pF

var tau_rise = tau1*tau2/(tau2-tau1); 
var nFactor = 1/(Math.pow((tau1/tau2),(tau_rise/tau2))-Math.pow((tau1/tau2),(tau_rise/tau1))); 
var peakLatency = tau_rise*Math.log(tau2/tau1); 

var wDur = 200; 	//ms 
var n = wDur*1e-3/deltaT; 	//number of points in the graph 
var t = 0; 

var x = []; 
var vMembrane = [];
var y1 = []; 
var y2 = []; 
var g_syn_wave = []; 
var iMembrane = []; 
var synRecov = []; 
var synFRecov = []; 

vMembrane.push(E_leak*1e-3); 
synRecov.push(1); 
synFRecov.push(1); 
y1.push(0); 
y2.push(0); 

var nPulse = 1; 
if (extTrainFlag) {
	var pulseTrain = makePulseTrain(tStart,freq,nPulse,wDur);
} else {
	var pulseTrain = []; 
}

var synDep = 0.8; 	//depression amount 
var tauDepr = 0.7; 	//s 
var synFacil = 0; 	//facilitation amount 
var tauFacil = 0.1; 	//s 

var currDur = 120; 	//ms 
var currMag = 0.2; 	//nA 

var currInj = 0; 
var currSyn = 0; 
var totalCurr1 = 0; 
var totalCurr2 = 0; 

var spikingNow = 0; 
var AP_threshold = -40; 	//mV 
var AP_max = 20; 	//mV 
var AP_reset = -45; 	//mV 

var shiftCurrX = 0; 	//ms 

var old1 = []; 
var old2 = []; 
var old3 = []; 
var old4 = []; 


var usercurrFlag = document.getElementById("currInj"); 
var usersynFlag = document.getElementById("syn"); 
var userextTrainFlag = document.getElementById("extTrain"); 

var userlastFlag = document.getElementById("lastTrace"); 

var usercurrMag = document.getElementById("amplitude"); 
var userE_leak = document.getElementById("restingMemPo"); 
var usertauM = document.getElementById("membraneTimCo"); 
var usercurrDur = document.getElementById("duration"); 
var userR_input = document.getElementById("inputRes"); 
var userC_M = document.getElementById("membraneCap"); 
var userAP_threshold = document.getElementById("APthres"); 
var userAP_max = document.getElementById("APPeak"); 
var userAP_reset = document.getElementById("APReset"); 
var userwDur = document.getElementById("simulationDur"); 
var usershiftCurrX = document.getElementById("currentOffset"); 
var userfreq = document.getElementById("pulseFreq"); 
var usernPulse = document.getElementById("numberPulse"); 
var usertStart = document.getElementById("timeStart"); 
var userg_syn = document.getElementById("synapticConductance"); 
var userE_syn = document.getElementById("reversalPo"); 
var usertau1 = document.getElementById("tauOne"); 
var usertau2 = document.getElementById("tauTwo"); 
var usersynDep = document.getElementById("STD"); 
var usersynFacil = document.getElementById("STF"); 
var usertauDepr = document.getElementById("tauD"); 
var usertauFacil = document.getElementById("tauF"); 

window.addEventListener("load", (event) => {
  usercurrMag.value = currMag;
  userE_leak.value = E_leak;
  usertauM.value = tauM; 
  usercurrDur.value = currDur; 
  userR_input.value = R_input; 
  userC_M.value = C_M.toFixed(3); 
  userAP_threshold.value = AP_threshold; 
  userAP_max.value = AP_max; 
  userAP_reset.value = AP_reset; 
  userwDur.value = wDur; 
  usershiftCurrX.value = shiftCurrX; 
  userfreq.value = freq; 
  usernPulse.value = nPulse; 
  usertStart.value = tStart; 
  userg_syn.value = g_syn; 
  userE_syn.value = E_syn; 
  usertau1.value = tau1; 
  usertau2.value = tau2; 
  usersynDep.value = synDep; 
  usersynFacil.value = synFacil; 
  usertauDepr.value = tauDepr; 
  usertauFacil.value = tauFacil; 

  document.getElementById("presets").value = "L5_PC";	// Defines startup preset. JSj 7 Sep 2020
  submitpreset(event);
  
  graph(); 
  draw(); 
  
  old1 = vMembrane.slice(); 
  old2 = iMembrane.slice(); 
  old3 = synRecov.slice(); 
  old4 = synFRecov.slice(); 
});


var submit = function(event) {
  event.preventDefault(); 
  
  old1 = vMembrane.slice(); 
  old2 = iMembrane.slice(); 
  old3 = synRecov.slice(); 
  old4 = synFRecov.slice(); 
  
  wDur = parseFloat(userwDur.value); 
  E_leak = parseFloat(userE_leak.value);
  tauM = parseFloat(usertauM.value); 
  R_input = parseFloat(userR_input.value); //Updating R_input as trauncated form from input field (calculation precision concern)
  AP_threshold = parseFloat(userAP_threshold.value); 
  AP_reset = parseFloat(userAP_reset.value); 
  AP_max = parseFloat(userAP_max.value); 
  tau1 = parseFloat(usertau1.value); 
  tau2 = parseFloat(usertau2.value); 
  synDep = parseFloat(usersynDep.value); 
  synFacil = parseFloat(usersynFacil.value); 
  tauDepr = parseFloat(usertauDepr.value); 
  tauFacil = parseFloat(usertauFacil.value); 
  currDur = parseFloat(usercurrDur.value); 
  currMag = parseFloat(usercurrMag.value); 
  shiftCurrX = parseFloat(usershiftCurrX.value); 
  freq = parseFloat(userfreq.value); 
  nPulse = parseFloat(usernPulse.value); 
  tStart = parseFloat(usertStart.value); 
  g_syn = parseFloat(userg_syn.value); 
  E_syn = parseFloat(userE_syn.value); 
  spikingNow = 0; 
  if (usercurrFlag.checked) { 
	currFlag = 1; 
  } else { 
	currFlag = 0; 
  } 
  if (userextTrainFlag.checked) { 
	extTrainFlag = 1; 
	pulseTrain = makePulseTrain(tStart,freq,nPulse,wDur); 
  } else { 
	extTrainFlag = 0; 
	nPulse = 1; 
	pulseTrain = []; 
  } 
  if (usersynFlag.checked) { 
	synFlag = 1; 
  } else { 
	synFlag = 0; 
  } 

  C_M = tauM*1e-3/(R_input*1e6)*1e12; 
  userC_M.value = C_M.toFixed(3); 
  tau_rise = tau1*tau2/(tau2-tau1); 
  nFactor = 1/(Math.pow((tau1/tau2),(tau_rise/tau2))-Math.pow((tau1/tau2),(tau_rise/tau1))); 
  peakLatency = tau_rise*Math.log(tau2/tau1); 
  n = wDur*1e-3/deltaT; 
  
  graph();
  draw(); 
}


var graph = function() { 
  t = 0; 
  x = []; 
  vMembrane = [];
  y1 = []; 
  y2 = []; 
  g_syn_wave = []; 
  iMembrane = []; 
  synRecov = []; 
  synFRecov = []; 

  x.push(t); 
  vMembrane.push(E_leak*1e-3); 
  synRecov.push(1); 
  synFRecov.push(1); 
  y1.push(0); 
  y2.push(0); 
  
  if (vMembrane[0] >= AP_threshold*1e-3) {
	  spikingNow = 1; 
  } 
  
  for (var i = 0; i < n-1; i++) {
	t += deltaT; 
	x.push(t); 
    if (synFlag) { 
		if (PulseNow(t,deltaT,pulseTrain)) { 
			y1.push(synRecov[i]*synFRecov[i]); 
			y2.push(synRecov[i]*synFRecov[i]+y2[i]); 
			synRecov.push(synRecov[i]*synDep); 
			synFRecov.push(synFRecov[i]+synFacil); 
		} else { 
			y1.push(y1[i] + deltaT/(tau1*1e-3)*(-y1[i])); 
			y2.push(y2[i] + deltaT/(tau2*1e-3)*(-y2[i])); 
			synRecov.push(synRecov[i] + deltaT/tauDepr*(1-synRecov[i])); 
			synFRecov.push(synFRecov[i]+ deltaT/tauFacil*(1-synFRecov[i])); 
		} 
	}
	currInj = 0; 
	for (var j = 0; j < nPulse; j++) { 
		if ((t >= 1/freq*j + tStart*1e-3 + shiftCurrX*1e-3) && (t < 1/freq*j + tStart*1e-3 + currDur*1e-3 + shiftCurrX*1e-3)) { 
			currInj += currMag*1e-9; 
		} 
	} 
	g_syn_wave.push(g_syn*1e-9*(nFactor*(y2[i]-y1[i]))); 
	currSyn = (vMembrane[i]-E_syn*1e-3)*g_syn_wave[i]; 
	totalCurr1 = 0; 
	totalCurr2 = 0; 
	if (currFlag) { 
		totalCurr1 += currInj; 
		totalCurr2 -= currInj; 
	}
	if (synFlag) { 
		totalCurr1 += currSyn; 
		totalCurr2 += currSyn; 
	} 
	iMembrane.push(totalCurr1*1e12); 
	if (spikingNow) { 
		vMembrane.push(AP_reset*1e-3); 
		spikingNow = 0; 
	} else { 
		vMembrane.push(vMembrane[i]+deltaT/(tauM*1e-3)*(R_input*1e6*(-totalCurr2)+E_leak*1e-3-vMembrane[i])); 
	} 
	if (vMembrane[i+1] > AP_threshold*1e-3) { 
		vMembrane[i+1] = AP_max*1e-3; 
		spikingNow = 1; 
	} 
  } 
  
  for (var i = 0; i < vMembrane.length; i++) { 
	vMembrane[i] = vMembrane[i]*1e3; 
  } 
  
}

var draw = function() {
  var V_membrane = {
  x: x,
  y: vMembrane,
  type: 'scatter', 
  name: 'V_membrane', 
  line: {
    color: 'rgb(219, 64, 82)',
    width: 2
    }
  };
  
  var oldV_membrane = {
  x: x,
  y: old1,
  type: 'scatter', 
  mode: 'lines',
  name: 'last trace', 
  line: {
    color: 'rgb(255, 165, 0)',
    width: 1.5,
	dash: 'dot'
    }
  };

  var i_membrane = {
    x: x,
    y: iMembrane,
    type: 'scatter', 
    name: 'i_membrane', 
    line: {
      color: 'rgb(55, 128, 191)',
      width: 2
    }
  }; 
  
  var oldi_membrane = {
    x: x,
    y: old2,
    type: 'scatter', 
	mode: 'lines', 
    name: 'last trace', 
    line: {
      color: 'rgb(73, 246, 252)',
      width: 1.5, 
	  dash: 'dot'
    }
  }; 
  
  var shortTermD = {
    x: x,
    y: synRecov,
    type: 'scatter', 
    name: 'shortTermD', 
    line: {
      color: 'rgb(124, 252, 0)',
      width: 2
    }
  };
  
  var oldshortTermD = {
    x: x,
    y: old3,
    type: 'scatter', 
	mode: 'lines', 
    name: 'last trace', 
    line: {
      color: 'rgb(124, 252, 0)',
      width: 2, 
	  dash: 'dot' 
    }
  };
  
  var shortTermF = {
    x: x,
    y: synFRecov,
    type: 'scatter', 
    name: 'shortTermF', 
    line: {
      color: 'rgb(255, 165, 0)',
      width: 2
    }
  };
  
  var oldshortTermF = {
    x: x,
    y: old4,
    type: 'scatter', 
	mode: 'lines', 
    name: 'last trace', 
    line: {
      color: 'rgb(255, 165, 0)',
      width: 2, 
	  dash: 'dot'
    }
  };
  
  if (userlastFlag.checked) {
	  Plotly.newPlot('myDiv', [V_membrane, oldV_membrane], layout, config); 
	  Plotly.newPlot('myDiv2', [i_membrane, oldi_membrane], layout2, config);
	  Plotly.newPlot('myDiv3', [shortTermD, oldshortTermD], layout3, config);  
	  Plotly.newPlot('myDiv4', [shortTermF, oldshortTermF], layout4, config);  
  } else {
	  Plotly.newPlot('myDiv', [V_membrane], layout, config); 
	  Plotly.newPlot('myDiv2', [i_membrane], layout2, config);
	  Plotly.newPlot('myDiv3', [shortTermD], layout3, config);  
	  Plotly.newPlot('myDiv4', [shortTermF], layout4, config);  
  }
    
}

var submitRin = function(event) {
  event.preventDefault(); 
  
  old1 = vMembrane.slice(); 
  old2 = iMembrane.slice(); 
  old3 = synRecov.slice(); 
  old4 = synFRecov.slice(); 
  
  wDur = parseFloat(userwDur.value); 
  E_leak = parseFloat(userE_leak.value);
  C_M = parseFloat(userC_M.value); //Updating C_M as trauncated form from input field (calculation precision concern)
  R_input = parseFloat(userR_input.value); 
  AP_threshold = parseFloat(userAP_threshold.value); 
  AP_reset = parseFloat(userAP_reset.value); 
  AP_max = parseFloat(userAP_max.value); 
  tau1 = parseFloat(usertau1.value); 
  tau2 = parseFloat(usertau2.value); 
  synDep = parseFloat(usersynDep.value); 
  synFacil = parseFloat(usersynFacil.value); 
  tauDepr = parseFloat(usertauDepr.value); 
  tauFacil = parseFloat(usertauFacil.value); 
  currDur = parseFloat(usercurrDur.value); 
  currMag = parseFloat(usercurrMag.value); 
  shiftCurrX = parseFloat(usershiftCurrX.value); 
  freq = parseFloat(userfreq.value); 
  nPulse = parseFloat(usernPulse.value); 
  tStart = parseFloat(usertStart.value); 
  g_syn = parseFloat(userg_syn.value); 
  E_syn = parseFloat(userE_syn.value); 
  spikingNow = 0; 
  if (usercurrFlag.checked) { 
	currFlag = 1; 
  } else { 
	currFlag = 0; 
  } 
  if (userextTrainFlag.checked) { 
	extTrainFlag = 1; 
	pulseTrain = makePulseTrain(tStart,freq,nPulse,wDur); 
  } else { 
	extTrainFlag = 0; 
	nPulse = 1; 
	pulseTrain = []; 
  } 
  if (usersynFlag.checked) { 
	synFlag = 1; 
  } else { 
	synFlag = 0; 
  } 

  tauM = (C_M*1e-12)*(R_input*1e6)*1e3; 
  usertauM.value = tauM.toFixed(3); 
  tau_rise = tau1*tau2/(tau2-tau1); 
  nFactor = 1/(Math.pow((tau1/tau2),(tau_rise/tau2))-Math.pow((tau1/tau2),(tau_rise/tau1))); 
  peakLatency = tau_rise*Math.log(tau2/tau1); 
  n = wDur*1e-3/deltaT; 
  
  graph();
  draw(); 
}

var submitCM = function(event) {
  event.preventDefault(); 
  
  old1 = vMembrane.slice(); 
  old2 = iMembrane.slice(); 
  old3 = synRecov.slice(); 
  old4 = synFRecov.slice(); 
  
  wDur = parseFloat(userwDur.value); 
  E_leak = parseFloat(userE_leak.value);
  tauM = parseFloat(usertauM.value); //Updating tauM as trauncated form from input field (calculation precision concern)
  C_M = parseFloat(userC_M.value); 
  AP_threshold = parseFloat(userAP_threshold.value); 
  AP_reset = parseFloat(userAP_reset.value); 
  AP_max = parseFloat(userAP_max.value); 
  tau1 = parseFloat(usertau1.value); 
  tau2 = parseFloat(usertau2.value); 
  synDep = parseFloat(usersynDep.value); 
  synFacil = parseFloat(usersynFacil.value); 
  tauDepr = parseFloat(usertauDepr.value); 
  tauFacil = parseFloat(usertauFacil.value); 
  currDur = parseFloat(usercurrDur.value); 
  currMag = parseFloat(usercurrMag.value); 
  shiftCurrX = parseFloat(usershiftCurrX.value); 
  freq = parseFloat(userfreq.value); 
  nPulse = parseFloat(usernPulse.value); 
  tStart = parseFloat(usertStart.value); 
  g_syn = parseFloat(userg_syn.value); 
  E_syn = parseFloat(userE_syn.value); 
  spikingNow = 0; 
  if (usercurrFlag.checked) { 
	currFlag = 1; 
  } else { 
	currFlag = 0; 
  } 
  if (userextTrainFlag.checked) { 
	extTrainFlag = 1; 
	pulseTrain = makePulseTrain(tStart,freq,nPulse,wDur); 
  } else { 
	extTrainFlag = 0; 
	nPulse = 1; 
	pulseTrain = []; 
  } 
  if (usersynFlag.checked) { 
	synFlag = 1; 
  } else { 
	synFlag = 0; 
  } 

  R_input = tauM*1e-3/(C_M*1e-12)*1e-6; 
  userR_input.value = R_input.toFixed(3); 
  tau_rise = tau1*tau2/(tau2-tau1); 
  nFactor = 1/(Math.pow((tau1/tau2),(tau_rise/tau2))-Math.pow((tau1/tau2),(tau_rise/tau1))); 
  peakLatency = tau_rise*Math.log(tau2/tau1); 
  n = wDur*1e-3/deltaT; 
  
  graph();
  draw(); 
}


document.getElementById("currInj").addEventListener("change", submit); 
document.getElementById("syn").addEventListener("change", submit); 
document.getElementById("extTrain").addEventListener("change", submit); 

document.getElementById("amplitude").addEventListener("change", submit); 
document.getElementById("restingMemPo").addEventListener("change", submit); 
document.getElementById("membraneTimCo").addEventListener("change", submit); 
document.getElementById("duration").addEventListener("change", submit); 
document.getElementById("inputRes").addEventListener("change", submitRin); 
document.getElementById("membraneCap").addEventListener("change", submitCM); 
document.getElementById("APthres").addEventListener("change", submit); 
document.getElementById("APPeak").addEventListener("change", submit); 
document.getElementById("APReset").addEventListener("change", submit); 
document.getElementById("simulationDur").addEventListener("change", submit); 
document.getElementById("currentOffset").addEventListener("change", submit); 
document.getElementById("pulseFreq").addEventListener("change", submit); 
document.getElementById("numberPulse").addEventListener("change", submit); 
document.getElementById("timeStart").addEventListener("change", submit); 
document.getElementById("synapticConductance").addEventListener("change", submit); 
document.getElementById("reversalPo").addEventListener("change", submit); 
document.getElementById("tauOne").addEventListener("change", submit); 
document.getElementById("tauTwo").addEventListener("change", submit); 
document.getElementById("STD").addEventListener("change", submit); 
document.getElementById("STF").addEventListener("change", submit); 
document.getElementById("tauD").addEventListener("change", submit); 
document.getElementById("tauF").addEventListener("change", submit); 

document.getElementById("presets").addEventListener("change", submitpreset); 
document.getElementById("lastTrace").addEventListener("change", draw); 


function JSdict() {
    this.Keys = [];
    this.Values = [];
}

// Check if dictionary extensions aren't implemented yet.
// Returns value of a key
if (!JSdict.prototype.getVal) {
    JSdict.prototype.getVal = function (key) {
        if (key == null) { 
			window.alert("Preset name cannot be null");
            return;
        }
        for (var i = 0; i < this.Keys.length; i++) {
            if (this.Keys[i] == key) {
                return this.Values[i];
            }
        } 
		window.alert("Preset not found!");
        return;
    }
}


// Check if dictionary extensions aren't implemented yet.
// Adds a unique key value pair
if (!JSdict.prototype.add) {
    JSdict.prototype.add = function (key, val) {
        // Allow only strings or numbers as keys
        if (typeof (key) == "number" || typeof (key) == "string") {
            if (key == null || val == null) {
				window.alert("Preset name or parameter settings cannot be null"); 
                return;
            }
            if (keysLength != valsLength) {
				window.alert("Dictionary inconsistent. Keys length don't match values!"); 
                return;
            }
            var keysLength = this.Keys.length;
            var valsLength = this.Values.length; 
            for (var i = 0; i < keysLength; i++) {
                if (this.Keys[i] == key) {
					this.Values[i] = val;
					window.alert("Preset updated");
                    return;
                }
            }
            this.Keys.push(key);
            this.Values.push(val);
			return; 
        }
        else {
			window.alert("Only number or string can be Preset names!"); 
            return;
        }
    }
}

// Check if dictionary extensions aren't implemented yet.
// Removes a key value pair
if (!JSdict.prototype.remove) {
    JSdict.prototype.remove = function (key) {
        if (key == null) {
			window.alert("Preset name cannot be null");
            return;
        }
        if (keysLength != valsLength) {
            return "Dictionary inconsistent. Keys length don't match values!";
        }
        var keysLength = this.Keys.length;
        var valsLength = this.Values.length;
        var flag = false;
        for (var i = 0; i < keysLength; i++) {
            if (this.Keys[i] == key) {
                this.Keys.shift(key);
                this.Values.shift(this.Values[i]);
                flag = true;
                break;
            }
        }
        if (!flag) { 
			window.alert("Preset does not exist");
            return;
        }
    }
}

function submitpreset(event) {
  var p = document.getElementById("presets").value; 
  if (p == "NMJ") { 
	usercurrFlag.checked = false; 
	userextTrainFlag.checked = true; 
	usersynFlag.checked = true; 
	usercurrMag.value = 0.1;
	userE_leak.value = -90;
	usertauM.value = 1; 
	usercurrDur.value = 120; 
	userR_input.value = 1; 
	userAP_threshold.value = 140; 
	userAP_max.value = 20; 
	userAP_reset.value = -45; 
	userwDur.value = 80; 
	usershiftCurrX.value = 0; 
	userfreq.value = 50; 
	usernPulse.value = 1; 
	usertStart.value = 30; 
	userg_syn.value = 5000; 
	userE_syn.value = 0; 
	usertau1.value = 0.1; 
	usertau2.value = 5.26; 
	usersynDep.value = 0.8; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.7; 
	usertauFacil.value = 0.1; 
  } 	  
  if (p == "L5_PC") { 
	usercurrFlag.checked = true; 
	userextTrainFlag.checked = false; 
	usersynFlag.checked = true; 
	usercurrMag.value = 0.1;
	userE_leak.value = -70;
	usertauM.value = 25; 
	usercurrDur.value = 120; 
	userR_input.value = 200; 
	userAP_threshold.value = -40; 
	userAP_max.value = 20; 
	userAP_reset.value = -45; 
	userwDur.value = 200; 
	usershiftCurrX.value = 0; 
	userfreq.value = 50; 
	usernPulse.value = 1; 
	usertStart.value = 30; 
	userg_syn.value = 1; 
	userE_syn.value = 0; 
	usertau1.value = 0.1; 
	usertau2.value = 5.26; 
	usersynDep.value = 0.8; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.7; 
	usertauFacil.value = 0.1; 
  } 	  
  if (p == "L5_PC_synapse") { 
	usercurrFlag.checked = false; 
	userextTrainFlag.checked = true; 
	usersynFlag.checked = true; 
	usercurrMag.value = 0.1;
	userE_leak.value = -70;
	usertauM.value = 25; 
	usercurrDur.value = 120; 
	userR_input.value = 200; 
	userAP_threshold.value = -40; 
	userAP_max.value = 20; 
	userAP_reset.value = -45; 
	userwDur.value = 600; 
	usershiftCurrX.value = 0; 
	userfreq.value = 10; 
	usernPulse.value = 5; 
	usertStart.value = 30; 
	userg_syn.value = 1; 
	userE_syn.value = 0; 
	usertau1.value = 0.1; 
	usertau2.value = 5.26; 
	usersynDep.value = 0.6; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.221; 
	usertauFacil.value = 0.169; 
  } 	  
  if (p == "L5_PC_lnF") { 
	usercurrFlag.checked = false; 
	userextTrainFlag.checked = true; 
	usersynFlag.checked = true; 
	usercurrMag.value = 0.1;
	userE_leak.value = -70;
	usertauM.value = 25; 
	usercurrDur.value = 120; 
	userR_input.value = 200; 
	userAP_threshold.value = -40; 
	userAP_max.value = 20; 
	userAP_reset.value = -45; 
	userwDur.value = 200; 
	usershiftCurrX.value = 0; 
	userfreq.value = 50; 
	usernPulse.value = 5; 
	usertStart.value = 30; 
	userg_syn.value = 12.5; 
	userE_syn.value = 0; 
	usertau1.value = 0.1; 
	usertau2.value = 5.26; 
	usersynDep.value = 1; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.7; 
	usertauFacil.value = 0.1; 
  } 	  
  if (p == "L5_PC_synapseNoDep") { 
	usercurrFlag.checked = false; 
	userextTrainFlag.checked = true; 
	usersynFlag.checked = true; 
	usercurrMag.value = 0.1;
	userE_leak.value = -70;
	usertauM.value = 25; 
	usercurrDur.value = 120; 
	userR_input.value = 200; 
	userAP_threshold.value = -40; 
	userAP_max.value = 20; 
	userAP_reset.value = -45; 
	userwDur.value = 200; 
	usershiftCurrX.value = 0; 
	userfreq.value = 10; 
	usernPulse.value = 5; 
	usertStart.value = 30; 
	userg_syn.value = 1; 
	userE_syn.value = 0; 
	usertau1.value = 0.1; 
	usertau2.value = 5.26; 
	usersynDep.value = 1; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.221; 
	usertauFacil.value = 0.169; 
  } 	  
  if (p == "L5_MC_synapse") { 
	usercurrFlag.checked = false; 
	userextTrainFlag.checked = true; 
	usersynFlag.checked = true; 
	usercurrMag.value = 0.1;
	userE_leak.value = -70;
	usertauM.value = 25; 
	usercurrDur.value = 120; 
	userR_input.value = 200; 
	userAP_threshold.value = -40; 
	userAP_max.value = 20; 
	userAP_reset.value = -45; 
	userwDur.value = 600; 
	usershiftCurrX.value = 0; 
	userfreq.value = 10; 
	usernPulse.value = 5; 
	usertStart.value = 30; 
	userg_syn.value = 1; 
	userE_syn.value = 0; 
	usertau1.value = 1.4; 
	usertau2.value = 5.26; 
	usersynDep.value = 1; 
	usersynFacil.value = 0.5; 
	usertauDepr.value = 0.221; 
	usertauFacil.value = 0.169; 
  } 	  
  if (p == "L5_BC") { 
	usercurrFlag.checked = false; 
	userextTrainFlag.checked = true; 
	usersynFlag.checked = true; 
	usercurrMag.value = 0.1;
	userE_leak.value = -70;
	usertauM.value = 10; 
	usercurrDur.value = 120; 
	userR_input.value = 160; 
	userAP_threshold.value = -37; 
	userAP_max.value = 20; 
	userAP_reset.value = -80; 
	userwDur.value = 200; 
	usershiftCurrX.value = 0; 
	userfreq.value = 10; 
	usernPulse.value = 1; 
	usertStart.value = 30; 
	userg_syn.value = 13.0258; 
	userE_syn.value = 0; 
	usertau1.value = 1.4; 
	usertau2.value = 3; 
	usersynDep.value = 1; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.221; 
	usertauFacil.value = 0.169; 
  } 	  
  if (p == "L5_BC_synapse") { 
	usercurrFlag.checked = false; 
	userextTrainFlag.checked = true; 
	usersynFlag.checked = true; 
	usercurrMag.value = 0.1;
	userE_leak.value = -70;
	usertauM.value = 10; 
	usercurrDur.value = 120; 
	userR_input.value = 160; 
	userAP_threshold.value = -37; 
	userAP_max.value = 20; 
	userAP_reset.value = -80; 
	userwDur.value = 200; 
	usershiftCurrX.value = 0; 
	userfreq.value = 30; 
	usernPulse.value = 5; 
	usertStart.value = 30; 
	userg_syn.value = 1; 
	userE_syn.value = 0; 
	usertau1.value = 1.4; 
	usertau2.value = 3; 
	usersynDep.value = 0.25; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.44; 
	usertauFacil.value = 0.169; 
  } 	  
  if (p == "L5_BC_slow") { 
	usercurrFlag.checked = false; 
	userextTrainFlag.checked = true; 
	usersynFlag.checked = true; 
	usercurrMag.value = 0.1;
	userE_leak.value = -70;
	usertauM.value = 10; 
	usercurrDur.value = 120; 
	userR_input.value = 160; 
	userAP_threshold.value = -37; 
	userAP_max.value = 20; 
	userAP_reset.value = -80; 
	userwDur.value = 200; 
	usershiftCurrX.value = 0; 
	userfreq.value = 10; 
	usernPulse.value = 1; 
	usertStart.value = 30; 
	userg_syn.value = 11.152; 
	userE_syn.value = 0; 
	usertau1.value = 1.4; 
	usertau2.value = 5; 
	usersynDep.value = 1; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.221; 
	usertauFacil.value = 0.169; 
  } 	  
  if (p == "L5_BC_train") { 
	usercurrFlag.checked = false; 
	userextTrainFlag.checked = true; 
	usersynFlag.checked = true; 
	usercurrMag.value = 0.1;
	userE_leak.value = -70;
	usertauM.value = 10; 
	usercurrDur.value = 120; 
	userR_input.value = 160; 
	userAP_threshold.value = -37.005; 
	userAP_max.value = 20; 
	userAP_reset.value = -80; 
	userwDur.value = 200; 
	usershiftCurrX.value = 0; 
	userfreq.value = 10; 
	usernPulse.value = 1; 
	usertStart.value = 30; 
	userg_syn.value = 22.0934; 
	userE_syn.value = 0; 
	usertau1.value = 1.4; 
	usertau2.value = 3; 
	usersynDep.value = 1; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.221; 
	usertauFacil.value = 0.169; 
  } 	  
  if (p == "L5_BC_slow_train") { 
	usercurrFlag.checked = false; 
	userextTrainFlag.checked = true; 
	usersynFlag.checked = true; 
	usercurrMag.value = 0.1;
	userE_leak.value = -70;
	usertauM.value = 10; 
	usercurrDur.value = 120; 
	userR_input.value = 160; 
	userAP_threshold.value = -37.005; 
	userAP_max.value = 20; 
	userAP_reset.value = -80; 
	userwDur.value = 200; 
	usershiftCurrX.value = 0; 
	userfreq.value = 10; 
	usernPulse.value = 1; 
	usertStart.value = 30; 
	userg_syn.value = 17.7621; 
	userE_syn.value = 0; 
	usertau1.value = 1.4; 
	usertau2.value = 5; 
	usersynDep.value = 1; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.221; 
	usertauFacil.value = 0.169; 
  } 	  
  if (p == "pulse_summation") { 
	usercurrFlag.checked = true; 
	userextTrainFlag.checked = true; 
	usersynFlag.checked = false; 
	usercurrMag.value = 0.32;
	userE_leak.value = -70;
	usertauM.value = 25; 
	usercurrDur.value = 10; 
	userR_input.value = 200; 
	userAP_threshold.value = -40; 
	userAP_max.value = 20; 
	userAP_reset.value = -45; 
	userwDur.value = 300; 
	usershiftCurrX.value = 0; 
	userfreq.value = 10; 
	usernPulse.value = 3; 
	usertStart.value = 30; 
	userg_syn.value = 1; 
	userE_syn.value = 0; 
	usertau1.value = 0.1; 
	usertau2.value = 5.26; 
	usersynDep.value = 0.8; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.7; 
	usertauFacil.value = 0.1; 
  } 	  
  if (p == "exponential") { 
	usercurrFlag.checked = true; 
	userextTrainFlag.checked = false; 
	usersynFlag.checked = false; 
	usercurrMag.value = 130;
	userE_leak.value = 0;
	usertauM.value = 25; 
	usercurrDur.value = 1; 
	userR_input.value = 200; 
	userAP_threshold.value = 2e6; 
	userAP_max.value = 20; 
	userAP_reset.value = -45; 
	userwDur.value = 200; 
	usershiftCurrX.value = 0; 
	userfreq.value = 50; 
	usernPulse.value = 1; 
	usertStart.value = 30; 
	userg_syn.value = 1; 
	userE_syn.value = 0; 
	usertau1.value = 0.1; 
	usertau2.value = 5.26; 
	usersynDep.value = 0.8; 
	usersynFacil.value = 0; 
	usertauDepr.value = 0.7; 
	usertauFacil.value = 0.1; 
  } 	  
  submit(event); 
} 
