# CCpol23+(2025) MBMD

This repository contains the reference implementation of the **CCpol23+ (2025)** potential energy surface (PES) for water.  
The model improves the description of short- and long-range three-body interactions and extends the framework introduced in:

**Ommair et al., “How Long-Range Are Three-Body “Exchange” Interactions in Liquid Water?,” J. Chem. Theory Comput. (2026) 22 (6)**  
https://doi.org/10.1021/acs.jctc.5c02107

**Current version written by:** *Ommair Ishaque (Nov 8, 2025)*

---

## 📘 Overview

CCpol23+ (2025) computes three body interaction energies of water trimers using an improved $F_{S^3}$ potential energy surface will remaining part of the potential is same as given in Góra et al using the same parameters. Another fit is made with improved optimized $F_{S^2}+F_{S^3}$ form. For both fit parameters are given in the input directories.  
The trimer configuration space is separated using the average monomer separation:

### PES Regions

| PES Component | Region Definition | Purpose |
|--------------|-------------------|---------|
| **Short-range** | Ravg < 5.5 Å | Exchange-repulsion and short-range induction |
| **Long-range** | Ravg > 4.5 Å | Smooth decay of long-range 3-body terms |
| **Overlap** | 4.5 Å – 5.5 Å | Blended SR/LR energy using switching function s(x) |
