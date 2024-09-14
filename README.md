# Gene Expression Analysis of a Glioblastoma dataset 

Authors (@slack): Manal Agdada (@Manal), Pascal Onaho (@PascalOnaho), Rahma Nabil Sallam (@rahmanabil2002), Hagar Haitham Elazab (@HBONH33), Salma ismail (@Salmaismail28), Idahosa Clinton (@doc_Idahosa), YaraHaitham (@YaraHaitham), Izuchukwu Obilom (@Zucchini)

Github Repo: https://github.com/manal-agdada/Stage-2---HackBio-internship-

## Software used for analysis

R studio 4.4.0

## Website used for enrichment analysis

ShinyGO http://bioinformatics.sdstate.edu/go/


## Introduction

This study analyzes and visualizes a gene expression dataset through heatmaps and functional enrichment analysis to explore patterns in differentially expressed genes and their biological significance. The dataset chosen comprises the top 500+ differentially expressed genes in glioblastoma under different conditions.

## Heatmap generation

We generated heatmaps for the dataset using the `heatmap.2()` function from the gplots package in R, employing two color palettes: a sequential palette to show gradient changes (Figure 1A), and a diverging palette to highlight deviations from a central value (Figure 1B). Sequential palettes work best for data with a progression from low to high values, whereas diverging palettes suit data with both positive and negative values. Three heatmap variants were then created: one clustering only genes (Figure 2A), one clustering only samples (Figure 2B), and one clustering both genes and samples (Figure 2C). This helps to assess relationships and patterns in the data, whether considering genes, samples, or both.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfVT2SGfXarjeq8pGCWX2JWuVxCgsn71fG4iWUug-hKwyPZdQNtcjtNyJfRcfyqvuYmuz3YV1MN8q_ADHEkuoRSi4jn8QIGdpNvuBClXkoNu4IT1HQDE6ghJA4tzcTC0XusrvTHaby85f15PU2eN79Otrkd?key=zicOZvtXSkbx22_fpayXbg)

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdqpg4rn-es493VonB-x7VyrKFHyWEjbe-atT3GRJ_djsFuiFRaCZ296DKEzIphAP_2LVpIWFBbxVZrDbimgFaNoPQGg_7PpRLgBdlJlL0v5dtCsGtknSZPdcfpN-RNeCUZoWQw3RR8icIjQay1HPE_8_rB?key=zicOZvtXSkbx22_fpayXbg)

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXficvm8Y_QzTA1Dh4Ke82dP6wmtPVdh_WLjdWmkfeGLnxYMg7ZC3CdNFKc2BAejs858V6f6xzHe9kD0o4ft51KYF943SwbIRaniGlTcpe8G6XLukInroGckG80dB0JH-1Nyc7j_3as46gVKb_sj3JW2gj5j?key=zicOZvtXSkbx22_fpayXbg)

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeQ2rH88MJz63hxmrLrqqnOiwwZZ-4Jt66xnwoTMwRKJyQQAZ400LjE8T_ZTTVpak5SBEH6Sfisu6J6xEeyO_jnb5NQhAUscI7CC_VrTeJFWw7UXWJqqSingIJQ1FU6YNFpuwRlaLnD0EQBgvNKMqaXRpo?key=zicOZvtXSkbx22_fpayXbg)

## Identification of up- and down-regulated genes

Careful assessment of the clustered heatmaps revealed similarities among samples categorized as “02A” and those labeled “01A” or “01B”, identifying two distinct groups of five similar samples each. Then, we calculated the fold change and relative p-values between them. Visualizing the p-values against the fold change helped justify the selection of cutoffs for subsequent gene subsetting, distinguishing up- and down-regulated genes in the dataset (Figure 3). The chosen cutoffs led to the identification of 57 up-regulated genes and 100 down-regulated genes.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeMBUmRHen7DNhRQgGtpvTjT5eQSaqWA1dLscp99p6odWZPBUrBqnZdOzbBbtWk95FXo5K5OjCG7YbCao0c46NLh_iLGdO6sgrUurOH0YkDIT7_snzIqaaaeNAxN92sT_dDd0zANwVc1gj7CVqshshwraY5?key=zicOZvtXSkbx22_fpayXbg)

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXc_UQ1CfDR5sI9QU_KxuIjoAwwVIirRFAju6ZZdz5oOTwqHyZX-zDjQyNWpxr--5-oRMfAfQMXWA9s3CcSkaHl7-elWRuBJ9eSMgqCQTgHPasCfasDrMIeRflkNQhdD6R1TRSaaMLcIHeLB3yQMOKRR4Jc3?key=zicOZvtXSkbx22_fpayXbg)

## Functional enrichment analysis

We performed a functional enrichment analysis of these deregulated genes to identify associated biological processes and pathways. ShinyGO was used to visualize enriched pathways using the GO biological process as a pathway database. We identified the top five deregulated pathways (Figure 4), with the top three being:

1\.  Regionalization refers to the process by which specific regions within an organism or tissue are defined and organized during development \[1].

2\.  Pattern specification process involves the establishment of spatial and temporal patterns during development guiding the formation of distinct structures and ensuring proper organ and tissue organization \[2].

3\.  Sensory organ development entails the processes involved in the formation and maturation of sensory organs, including their structure and function, essential for detecting and responding to environmental stimuli \[3].

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcGg9aQaXSFCCHMVxC4mdKFqJeTdBq4dLZr-_4rz4c7CBEdBRHLRrvWSVPLPpxyBS__EXFgo5zdwnY0_IB88FLfpkL0Dh_6S7gYdZpKk8Apvjz1lKW-N0s2L-Q4wB01KVq_mEp8_qMdhyYhBCJ1AogqM7f6?key=zicOZvtXSkbx22_fpayXbg)

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeNHts-O6VeVzKF2Taveiiv5YsQrVoa14M0aU-eo14k8aKCrAyS4X4NGfKSupln10kTuHtxT_3sZfiTmFzXGHI1bIDV8JpZsEKNCiYSQ6QrbrQniILZVeUlnJGyWRcs-pgtKR_yUIvr-T-tJFM57gn4DV5e?key=zicOZvtXSkbx22_fpayXbg)

## Conclusion

This analysis revealed key pathways linked to gene expression in glioblastoma. These findings enhance our understanding of the biological processes involved in this disease.

## References 

1\. Szczesniak R, Rice JL, Brokamp C, Ryan P, Pestian T, Ni Y, Andrinopoulou ER, Keogh RH, Gecili E, Huang R, Clancy JP, Collaco JM. Influences of environmental exposures on individuals living with cystic fibrosis. Expert Rev Respir Med. 2020 Jul;14(7):737-748.

2\. Ramos R, Swedlund B, Ganesan AK, Morsut L, Maini PK, Monuki ES, Lander AD, Chuong CM, Plikus MV. Parsing patterns: Emerging roles of tissue self-organization in health and disease. Cell. 2024 Jun 20;187(13):3165-3186.

3\. Ungefroren, H., Witte, D. and Lehnert, H. (2018), The role of small GTPases of the Rho/Rac family in TGF-β-induced EMT and cell motility in cancer. Dev. Dyn., 247: 451-461.


