# RACC-PM
This is a MATLAB code repository for the manuscript below. 
'A Robust Hybrid Algorithm for Personal Sound Zones: A Worst-Case Optimization Approach'.

#  Introduction
 There are four files in the project:
 * **codes** file  provides main function 'SFC.m' and seven algorithms implementation with some dependent functions of sound field control; Make sure the CVX  toolbox has been installed before running the code.
 * **datasets** file provides the room impulse response(RIR) in frequency domain and $\alpha$ used in our experiments; The ATF data is too large to upload and can be generated through rirGenerate.
 * **results** file is a collection of the algorithms performance in the paper.
 * **rirGenerate** file provides the way of generate RIR with three ways to add noise: adding Gaussian noise(snr: 15-25dB), position perturbation(0-0.05m), reverberation(0.3s-0.4s).
 # Additional Experiment: Adding position perturbation and reverberation levels in the paper.

To evaluate the performance of the proposed RACC-PM and other algorithms,  The resulting graphs are shown here.


<div align=center>
<img src=/results/EvaluationResultsPos.png width="400" >
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

The corresponding parameter settings are shown in TABLE I.

<p align="center">
<small>
TABLE I: Parameters settings on different datasets
</small>
</p>
<table border="1" width="500px" cellspacing="10" align="center">
<tr>
  <th align="center"> Dataset </th>
  <th align="center"> ACC-PM </th>
  <th align="center"> VAST-NF </th>
  <th align="center"> wc-RACC </th>	
  <th align="center"> POTDC-RACC </th>
  <th align="center"> RACC-PM </th>
</tr>
<tr>
  <td rowspan="1" align="center">Position perturbation</td>
  <td rowspan="3" align="center"> $\kappa=0.7$</th>
  <td rowspan="3" align="center"> $\mu = 1$, $V=L/2=8$</th>
  <td rowspan="3" align="center"> $\gamma_{\rm{B}}\approx\epsilon_B^2$, $\gamma_{\rm{D}}\approx\epsilon_D^2$</th>
  <td rowspan="3" align="center"> $\alpha_{l}, \alpha_{u}$ [34], $\eta = \epsilon_B$, $\gamma_{\rm{D}} = \gamma_{\rm{D}}$
     in wc-RACC</th>
  <td rowspan="3" align="center"> $\sqrt{e_w} = \Vert \mathbf{w}_{\rm{ACC-PM}} \Vert$, $\rho = 0.1$, $\mu = 1$, $\alpha = {\rm{AC}}_{\rm{ACC(True)}}$, 
    $\epsilon_{B} = 0.0001\sqrt{{\rm{tr}}(\mathbf{H}_{\rm{B}}^{\rm{H}}\mathbf{H}_{\rm{B}})}$, 
    $\epsilon_D = 0.0001\sqrt{{\rm{tr}}(\mathbf{H}_{\rm{D}}^{\rm{H}}\mathbf{H}_{\rm{D}})}$
</th>
</tr>


<tr>
  <td rowspan="1" align="center">Reverberation(0.3-0.4s)</td>
  
</tr>
<tr>
  <td rowspan="1" align="center">Reverberation(0.3-0.6s)</td>
</tr>
</table>
