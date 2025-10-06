## 16S rRNA Gene Marker Analysis

### Seawater DNA Extraction
- Environmental DNA (eDNA) was extracted from Sterivex cartridges following, using a protocol adapted from the *QIAGEN DNeasy Blood & Tissue Kit* (Qiagen, North Rhine–Westphalia, Germany).  
- Cartridges were thawed on ice, and the outlet was sealed with a Luer-lock cap to prevent leakage.  
- An extraction buffer (440 µL total) was prepared and added through the inlet, consisting of:  
  - 40 µL **Proteinase K** (provided by the kit)  
  - 200 µL **Buffer AL** (provided by the kit)   
  - 220 µL **PBS** (Gibco, NY, USA)  
- The inlet was then sealed with another Luer-lock cap (TrueCare, FL, USA), and both ends were wrapped in Parafilm.  
- Cartridges were placed horizontally on a rotator inside a hybridization incubator (Robbins Scientific Model 400, Sunnyvale, CA, USA) and incubated at **56 °C for 4 h** while rotating at 20 rpm.  
- To recover the lysate, the inlet was connected to a 2 mL microcentrifuge tube, secured inside a sterile 50 mL conical tube, and centrifuged at **5,000 × g for 2 min**.  
- The recovered lysate was mixed with **200 µL of 100% ethanol**, vortexed, and transferred to a *DNeasy spin column* (Qiagen).  
- The column was washed sequentially with **500 µL Buffer AW1** and **500 µL Buffer AW2**, each followed by centrifugation at **5,000 × g for 1 min**.  
- DNA was eluted in **100 µL Buffer AE** into a new 1.5 mL tube and stored at **−80 °C** until library preparation.

---

### Library Preparation and Sequencing
- The V4 hypervariable region of the 16S rRNA gene was amplified using the following primer pair:  
  - Forward primer 515F (5′-GTGYCAGCMGCCGCGGTAA-3′) (*Parada et al., 2016*)  
  - Reverse primer 806R (5′-GGACTACNVGGGTWTCTAAT-3′) (*Apprill et al., 2015*)  
- Primers were synthesized with Illumina adapter overhangs following the *Earth Microbiome Project* protocol (*Caporaso et al., 2012*).  

- **PCR reactions** (25 µL total volume) contained:  
  - 12.5 µL AccuStart II PCR ToughMix (2×) (Quanta BioSciences, Gaithersburg, MD, USA)  
  - 7.5 µL ultra-pure water 
  - 2.5 µL forward primer (10 µM)  
  - 2.5 µL **reverse primer (10 µM)
  - 1 µL template DNA
     
- **Thermal cycling conditions:**  
  - Initial denaturation: 94 °C for 3 min  
  - 25 cycles of:  
    - Denaturation at 94 °C for 45 s  
    - Annealing at 50 °C for 60 s  
    - Extension at 72 °C for 90 s  
  - Final extension: 72 °C for 10 min  
  - Hold: 4 °C  

- PCR products were visualized by 1.5% agarose gel electrophoresis (85 V, 120 min).  
- Amplicons were purified using AMPure XP Magnetic Beads (Beckman Coulter Life Sciences, IN, USA), following the manufacturer’s protocol.  
- DNA concentration was quantified using the Quant-iT 1X dsDNA High Sensitivity (HS) Assay Kit (Thermo Fisher Scientific, MA, USA) on a BioTek Synergy H1 multi-mode plate reader.  
- Libraries were pooled equimolarly and sequenced using Illumina MiSeq 2×250 bp paired-end runs at the *Oregon State University Center for Quantitative Life Sciences (CQLS)*.

---

### Notes 
- All reagents were handled in a laminar flow hood, and negative controls were included throughout the extraction and amplification processes.
