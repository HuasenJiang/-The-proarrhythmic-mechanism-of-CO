# -The-proarrhythmic-mechanism-of-CO
**This repository contains the code and experimental data for this study: *In silico* assessment of pharmacotherapy for carbon monoxide induced arrhythmias in healthy and failing human hearts, Huasen Jiang et al., 2022**
## Graphical Abstract
![abstract](https://user-images.githubusercontent.com/71884708/178404056-891523fd-c9ae-46ec-ada4-91767eba7be9.jpg)
## Abstract
**Background:** Carbon monoxide (CO) is gaining increased attention in the air pollution-induced arrhythmias. The severe cardiotoxic consequences of CO urgently require an effective pharmacotherapy to treat it. However, existing evidence demonstrates that CO can induce arrhythmias by directly affecting ion channels rather than the indirect pathway of tissue hypoxia, making it empirically difficult to find suitable drugs. 

**Objective:** To find an effective pharmacotherapy for the treatment of CO-induced arrhythmias through a virtual pathological tissue model that acts as an in silico drug-screening platform. 

**Methods:** Two pathological models describing the CO effects on healthy and failing hearts, respectively, were constructed as control baseline models. After this, we first assessed the efficacy of some common antiarrhythmic drugs like ranolazine, amiodarone, nifedipine, etc, by incorporating their ion channel-level effects into the cell model. Cellular biomarkers like action potential duration and tissue-level biomarkers such as the QT interval from pseudo-ECGs were obtained to assess the drug efficacy. In addition, we also evaluated multiple specific IKr activators in a similar way to multichannel blocking drugs, as the IKr activator showed great potency in dealing with the CO-induced pathological changes.

**Results:** Simulation results showed that the tested seven antiarrhythmic drugs failed to rescue the heart from CO-induced arrhythmias in terms of the action potential and the ECG manifestation. Some of them even worsened the condition of arrhythmogenesis. In contrast, IKr activators like HW-0168 effectively alleviated the proarrhythmic effects of CO.

**Conclusions:** Current antiarrhythmic drugs including the ranolazine suggested in previous studies did not achieve therapeutic effects for the cardiotoxicity of CO. The specific IKr activator, which is less concerned in current antiarrhythmic strategy, is a promising pharmacotherapy for the treatment of CO-induced arrhythmias. 
## Structure of the repository
***Baseline:** It contains code for single-cell and one-dimensional tissue models of healthy and failing hearts affected by CO.

***Data:** It contains data on CO pathogenicity, data on changes in ionic currents in failing hearts, and data on four specific IKr activators and six multichannel blockers.

***Drug Simulation:** It contains codes for single-cell and one-dimensional tissue models for two of the drugs studied in detail in experiments, ranolazine and HW-0168, and for one-dimensional tissue for six multichannel blockers.
