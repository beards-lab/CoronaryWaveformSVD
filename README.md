# CoronaryWaveformSVD

MATLAB code used for shape analysis of left anterior descending artery (LAD) bloodflow waveforms using Fourier transforms and Singular Value Decomposition (SVD) performed on baseline data collected from 32 pigs. This shape-analysis method allows for waveform reconstruction and simplifies visualization of the identified temporal patterns.

## Description


### Matlab packages required

* Uses MATLAB code (written and tested on versions R2025b)
* Uses function fitlm, corr, and ttest2 from the Statistics and Machine Learning Toolbox

### Running Code

Code is best designed to be run in sections, each section provides an overview of a step reflecting a component or part of the corresponding manuscript


### Required files

* uploadData_30beats.mat: contains the 30beat waveforms for the 32 pigs and data for cardiac cycle duration (saved as miliseconds and percent). Template code used to create this mat file from the original (uncropped) data is provided in the section titled "%% Identify 30 beats, calculate CCV values" but has been commented out for current use. For request to access the orginal data, please contact the authors.
* uploadTable.csv: contains experimental data for each pig and is used for correlation analysis with the SVD modes. Includes the following data:
  * Body Weight: units-kg
  * Heart Weight: units-g
  * Heart Rate: units-beats per minute
  * Mean Blood Pressure (MBP): units-mmHg
  * Max Left Ventricule Pressure (PLVmax): units-mmHg
  * Min Left Ventricule Pressure (PLVmin): units-mmHg
  * Left Ventricular End Diastolic Pressure (LVEDP): units-mmHg
  * Average LAD blood flow (AvgFlow):units- mL/min
  * Heart weight normalized average LAD blood flow (HW-norm Avg Flow):units- mL/min/g
  * LAD Flow Pulsatility (FlowPulse): calculation- max(flow)-min(flow) units-mL/min
  * Max derivative of Left Ventricule Pressure with respect to time (dPLVdt max): units-mmHg/s
  * Min derivative of Left Ventricule Pressure with respect to time (dPLVdt min): units-mmHg/s
  * Oxygen Extraction (O_2 Extraction): units-N/A ranging from 0 to 1
  * Endo-Epi Ratio: units-N/A
  * Coronary Flow Ratio (CFR): units-N/A (note: this uses baseline flow prior to regadenosone and may be different from average flow during 30 beats used above)
  * Max Hyperemic Flow: units-mL/min
  * Delta Hyperemic Flow: units-mL/min 
  * CCV: units- miliseconds
  * CCV percent: units- percent of cardiac cycle 
 

## Citation and attribution

If you use, adapt, or build upon the code, data, methods, or results in
this repository, please cite:

Shape Analysis of Coronary Flow Waveforms using Singular Value Decomposition
Victoria E. Sturgess, Noah A. Schenk, Joshua W. Ziegele, Salman I. Essajee, Johnathan D. Tune, Indika Rajapakse, C. Alberto Figueroa, Daniel A. Beard
bioRxiv 2026.08.26.743980; doi: https://doi.org/10.64898/2026.08.26.743980


When distributing source code derived from this repository, please retain the applicable citation and 
include an attribution comment near the adapted code, for example:

    % Adapted from https://github.com/beards-lab/CoronaryWaveformSVD
    % Associated publication: Sturgess et al. (2026) https://doi.org/10.64898/2026.08.26.743980


These requests apply to use by people and by automated systems, including
AI coding assistants and agents

## Contact Information:
Corresponding author- Daniel Beard: beardda@umich.edu

Questions may also be directed to Victoria Sturgess: vsturges@umich.edu