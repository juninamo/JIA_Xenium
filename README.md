# Spatial transcriptomics reveals immune–stromal crosstalk within the synovium of patients with juvenile idiopathic arthritis

Analytic codes for JIA-synovium Xenium data analysis

📄 **Published in *JCI Insight* (2026)**: [https://insight.jci.org/articles/view/198074](https://insight.jci.org/articles/view/198074)

R package for spatial neighborhood analysis (spatialCooccur) can be downloaded from [here](https://github.com/juninamo/spatialCooccur).

## 💾 Data availability
- Raw and processed Xenium data are available at GEO under accession [GSE326837](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326837).
- The RA CITE-seq dataset used for comparison is available from the AMP RA/SLE Network under controlled access (SynID: syn52297840) and is not redistributed here.

## 📁 Folder structure

```
├── codes/
│   ├── Overview.ipynb                               # (1) Overview of Xenium data, (2) Association test of clinical variables, and (3) Integrative analysis with GWAS data
│   ├── Spatial_Niche.ipynb                          # Spatial niche analysis
│   ├── Spatial_Neighborhood_Analysis_FibEndo.ipynb  # (1) Spatial Neighborhood analysis, and (2) Spatial association between fibroblasts and endothelial cells
│   ├── Colocalization_Score_MacroT.ipynb            # (1) Spatial association between macrophages and endothelial cells, and (2) Colocalization score analysis of macrophages and T cells
│   ├── TLS.ipynb                                    # Analysis of tertiary lymphoid structures (TLS)
│   ├── vsRA.ipynb                                   # Comparison with RA scRNA-seq data
│   └── JIA[1-9].html                                # Interactive plot of JIA Xenium data
│
└── README.md                                        # This file
└── LICENSE                                          # License file
```


## 📝 Citation 
Inamo J, Fierkens R, Clay MR, Jonsson AH, Lin C, Hayes K, Rogers N, Leach H, Yomogida K. Spatial transcriptomics reveals immune–stromal crosstalk within the synovium of patients with juvenile idiopathic arthritis. *JCI Insight*. 2026;11(1):e198074. doi:[10.1172/jci.insight.198074](https://doi.org/10.1172/jci.insight.198074)

Preprint: [*medRxiv* 2025](https://www.medrxiv.org/content/10.1101/2025.08.05.25332835v1), doi:[10.1101/2025.08.05.25332835](https://doi.org/10.1101/2025.08.05.25332835)

## 📬 Contact
For questions or issues related to this tutorial, please contact;

**Name:** Jun Inamo  
**Email:** juninamo@keio.jp  
**Affiliation:** Division of Rheumatology, Department of Internal Medicine, Keio University School of Medicine

The data presented here comes from the [Yomogida lab](https://www.yomogidalab.com/).

## ✅ Acknowledgments
This work was supported by the Uehara Memorial Foundation Postdoctoral Fellowship, a Grant-in-Aid for Japan Society for the Promotion of Science Overseas Research Fellows, the Mochida Memorial Foundation for Medical and Pharmaceutical Research (to J.I.), K08DK128544 (K.Y.). 
We acknowledge the Accelerating Medicines Partnership (AMP): Rheumatoid Arthritis and Systemic Lupus Erythematosus (AMP RA/SLE) Network for providing the data used in this study, specifically the RA CITE-seq dataset (SynID: syn52297840, Release 1.5/V7). 

## License
This repository is provided under the MIT License.

&nbsp;&nbsp;