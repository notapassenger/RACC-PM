# RACC-PM
This is a MATLAB code repository for the manuscript below. A Robust Hybrid Algorithm Using Worst-Case Optimization for Pesonal Sound Zones.

#  Introduction
 There are four files in the project:
 * **codes** file  provides main function 'SFC.m' and seven algorithms implementation with some dependent functions of sound field control; Make sure the CVX  toolbox has been installed before running the code.
 * **datasets** file provides the room impulse response(RIR) in frequency domain and $\alpha$ used in our experiments; The ATF data is too large to upload and can be generated through rirGenerate.
 * **results** file is a collection of the algorithms performance in the paper.
 * **rirGenerate** file provides the way of generate RIR with three ways to add noise: adding Gaussian noise(snr: 15-25dB), position disturbance(0-0.05m), reverberation(0.3s-0.4s).
 # Additional Experiment: Position perturbation and reverberation levels in the paper.

To evaluate the performance of the proposed $RACC-PM$ algorithm,  The resulting graphs are shown here.


<div align=center>
<img src="https://github.com/notapassenger/RACC-PM/tree/main/results/EvaluationResultsPos.png" width="400" >
</div>
<p align="center">
<small>
Fig.1 Evaluation results of position perturbation are added.
</small>
</p>

 <div align=center>
<img src="https://github.com/notapassenger/RACC-PM/tree/main/results/EvaluationResultsRev(0.3-0.4s).png" width="400" >
</div>
<p align="center">
<small>
Fig.2 Evaluation results with 0.3-0.4s reveberation level.
</small>
</p>

<div align=center>
<img src="https://github.com/notapassenger/RACC-PM/tree/main/results/EvaluationResultsRev(0.3-0.6s).png" width="400" >
</div>
<p align="center">
<small>
Fig.3 Evaluation results with 0.3-0.6s reveberation level.
</small>
</p>
