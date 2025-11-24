# Dynamic Entrainment of Neuronal Spiking — Code and Data

This repository accompanies the journal article **_Dynamic entrainment of neuronal spiking: A phase-plane analysis and a data-driven approach_**, which is currently under review. It contains all code and datasets used for data generation, machine learning model development, and figure production.

## Repository Structure

### DataGeneration
MATLAB scripts that were used to generate the synthetic dataset for training the machine learning model.

### MachineLearning
The generated dataset (`.csv`) and the Jupyter notebooks used to train, evaluate, and select the machine learning model.

### Plotting
MATLAB scripts for generating phase-plane plots, time-series plots, and the XPP code used to plot the bifurcation diagram.

---

## Workflow for Dynamic Entrainment

1. **Generate the training dataset** (in DataGeneration).

2. **Train and save the machine learning model** (in MachineLearning).

3. **Run** `MorrisLecar_After_First_Pulse.m` (in Plotting) for given **Amp** and **L**.  
   - Run with **item = 1 uncommented** to start from the steady state.

4. **Run** `NeuralNetwork3_savedModel_I0_80.ipynb` (in MachineLearning).  
   - Inputs: **Amp**, **L**, and the output from Step 3.  
   - The model predicts **tnext**, the required time gap to the second pulse for 1:1 entrainment.

5. **Update initial conditions** for the next pulse.  
   - Run `MorrisLecar_After_First_Pulse.m` again with **item = 1 commented**.  
   - Use `p2` from the previous run as the new initial voltage.  
   - Use the Morris–Lecar phase plane to find the corresponding **n** value.

6. **Iterate** between the MATLAB script and the machine learning model to compute a sequence of **t_next** values. From this, you can get the pulse onsets to achieve 1:1 entrainment.

7. **Visualize results** by inserting the computed onset times into the `pulse_starts` array in `DynamicEntrainment_Pulse_MorrisLecar_IC.m` (in Plotting), then run the script to generate the voltage time series.

---

If you use this code for research, please consider citing the associated article once published.
