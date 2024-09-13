<!--StartFragment-->


### **Gene Expression Analysis of a Glioblastoma dataset**

Authors (@slack): Manal Agdada (@Manal), Pascal Onaho (@PascalOnaho), Rahma Nabil Sallam (@rahmanabil2002), Hagar Haitham Elazab (@HBONH33), Salma ismail (@Salmaismail28), Idahosa Clinton (@doc\_Idahosa), YaraHaitham (@YaraHaitham), Izuchukwu Obilom (@Zucchini)

Github Repo: <https://github.com/manal-agdada/Stage-2---HackBio-internship-/blob/main/stage2.md>

Github code link: <https://github.com/manal-agdada/Stage-2---HackBio-internship-/blob/main/stage2.R>

This study analyzes and visualizes a gene expression dataset through heatmaps and functional enrichment analysis to explore patterns in differentially expressed genes and their biological significance. The dataset chosen comprises the top 500+ differentially expressed genes in glioblastoma under different conditions.


### **Heatmap generation**

We generated heatmaps for the dataset using the \`heatmap.2()\` function from the _gplots_ package in R, employing two color palettes: a sequential palette to show gradient changes (Figure 1A), and a diverging palette to highlight deviations from a central value (Figure 1B). Sequential palettes work best for data with a progression from low to high values, whereas diverging palettes suit data with both positive and negative deviations. Three heatmap variants were then created: one clustering only genes (Figure 2A), one clustering only samples (Figure 2B), and one clustering both genes and samples (Figure 2C). This helps to assess relationships and patterns in the data, whether considering genes, samples, or both.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfDhqB3AZ_CJ8f6D9WWE0xXft1el4UeTPEckYJ6zqiGFFxwBHewCUYzu0iJSdvrnbDxuD9_igfD3SqEtv3YAB2p3r6k-exa8QHrQ5jbE0x_NuhKg7M2XPD70R65ILXJI9PJNCrzlsDHX60DMqUbYE-Zs74T?key=zicOZvtXSkbx22_fpayXbg)

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdQjy5jGlIKqLAz7B8X-evuiHU6UEqNusnOA55H85OcAumAjr-ujzFy45orwQ9xYk5BzB9g3DpCWyguvhmvo2reKGHAUy4-tqc_hN6hsS_w5d06WEYzTnGeSb8GetDewg_LcWpFD3eI279heJxW0wyDdOBw?key=zicOZvtXSkbx22_fpayXbg)

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXei26jQcEC3m5Y3nu-Z9buoUQBhtNtFNEfTBNjvMvccVBfqB8_PVPnz4gQXTb5CBbmkVFChFCLfC55jsCfxuvzxrCFS8jfnvDMKp4kJweeb4CIxZAuWtwiiIYYpcEHSGFiMKWMxY8I90KOQWraW6h9CG2HV?key=zicOZvtXSkbx22_fpayXbg)

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfiMyNdTD0JYKXeqBEdOGk1OPU7g1Ek3REkcarb5gKEogTe1dnY6MWefTU1hWUGMCgEKonuGyVA5sYs5ktIwJpuvZw8z9f14jC1Dnlb2LdSG-AnqbqYNhe5u59ZPjMQvHoF_evn2Tqoe8tTrfhfgyOLNEQ?key=zicOZvtXSkbx22_fpayXbg)


### **Identification of up- and down-regulated genes**

Careful assessment of the clustered heatmaps revealed similarities among samples with the “02A” and those with “01A” or “01B”. Consulting the TCGA website (<https://www.cancer.gov/ccg/research/genome-sequencing/tcga>), we found that “02” refers to recurrent tumors, and “01” refers to primary tumors. Principal Component Analysis (PCA) confirmed that “01” samples are clustered tightly, while “02” samples clustered more loosely (Figure 3). This led to the identification of 19 upregulated and 70 downregulated genes through differential expression analysis using standard cut-offs for fold change and adjusted p-values.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeQVSFiNFUpgrR1v6Q44d3Da9MJyUP2oUW50WrJ-tfgk3gUHAkFp3DnxCpGJ9F5RDhbp9zwJld346mAXqoOyw7uNtGPFr-S4r2vYQR2uTT5LDwfXieONjWI_QcTCy-Huq8pMEihaCJyiyPkuNT_7aDQS89K?key=zicOZvtXSkbx22_fpayXbg)

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXc4RnETK6bS2DiGjci0vKJlHbqLNCz6ivEwtzopK0VFNTliOn7eY2-iIBEurYQFVdpN_nccuYNhdC6361DNnw_RReBy_H8rGmgD1dhMctaZ-3bKeZIr0mmOb5LFnN2PupBTAtzxW9wwM6UxUoUa8NT4WxwX?key=zicOZvtXSkbx22_fpayXbg)


### **Functional enrichment analysis** 

We performed a functional enrichment analysis of these deregulated genes to identify associated biological processes and pathways. ShinyGO was used to visualize enriched pathways using all available gene sets as a pathway database. We identified the top ten deregulated pathways (Figure 4), with the top three being:

1\.  The CXC Chemokine domain pathway involves proteins that recruit immune cells via CXCR proteins \[1].

2\.  The CXCR Chemokine Receptor Binding pathway involves G-protein coupled receptors that respond to specific molecules such as chemokines to regulate immune cell function \[1].

3\.  The MEISSER BRAIN HCP WITH H3K27ME3 pathway involves genes with H3K27 trimethylation marks in the brain \[3].

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXecPAthq7Ur3knd27i_41tSaT16vAyBgR9shZJJYg5R88pE9Oa4TuSEJxgAVEFkOYYI5rgfEyDuQ1v7vIE5Pj6YakqW_DRyIxaavnC5ppUOrHuhwUrtB_NwkxoBjzxNwt7EboXHq5PqbtTmS_IYDcjK81fh?key=zicOZvtXSkbx22_fpayXbg)

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXczVWQF4PFgvX4EFv2uaf5vUY-IWRC8y3nwpZwQWZmSiZcbT35JDvMaQE1xn7LWCp3kqYo3kLUWaEuiYWijk3W8j2nY-0gSmXGHxIrY2rWibd3IyABE-9DNlTCKjzhtfa76m4dl6Ek84eJmuE201CMv9grG?key=zicOZvtXSkbx22_fpayXbg)


### **Conclusion**

This analysis revealed key pathways linked to gene expression in glioblastoma. These findings enhance our understanding of the biological processes involved in this disease.


### **References** 

1\. García-Cuesta EM, Santiago CA, Vallejo-Díaz J, Juarranz Y, Rodríguez-Frade JM, Mellado M. **The Role of the CXCL12/CXCR4/ACKR3 Axis in Autoimmune Diseases**. Front Endocrinol (Lausanne). 2019 Aug 27;10:585.

2\. Meissner A, Mikkelsen TS, Gu H, Wernig M, Hanna J, Sivachenko A, Zhang X, Bernstein BE, Nusbaum C, Jaffe DB, Gnirke A, Jaenisch R, Lander ES. **Genome-scale DNA methylation maps of pluripotent and differentiated cells**. Nature. 2008 Aug 7;454(7205):766-70. 


<!--EndFragment-->
