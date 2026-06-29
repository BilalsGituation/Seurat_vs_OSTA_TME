Date released: 2026-06-29

Author: Bilal Asser

(Hardware needs minimum 16GB memory to run this)

## Purpose of this work:
This workflow provides a practical intro to Tumour Microenvironment (TME) analysis using Spatial Transcriptomics data (10X Visium), for learners who have some experience with single cell and spatial transcriptomics analysis in Seurat through the vignettes.

## Quick note on breast cancer TME:
DCIS is Ductal Carcinoma in situ. These cells are preinvasive cancer cells still confined to the milk ducts.

## Quick concept tutorial for later:
- Normalised Mutual Information (NMI) is a measure of the similarity of two images.
- NiftyReg is a visual pattern detection algorithm which finds and registers statistical transformations on a source image to make it more similar to a query image (increasing their NMI score).
(statistical transformation examples: rotation, scale, shear, translation)
    - If two images lack shared structural patterns, NiftyReg cannot artificially increase their similarity.

## Aims: 
1) Spatial transcriptomics workflow interoperation by translating the OSTA chapter 12 workflow to the Seurat framework.
2) Visually compare differences in spatial maps produced by our deconvolution using `Seurat::FindTransferAnchors()` and OSTA's deconvolution using RCTD.
3) Use NMI to quantify the change in spatial alignment between spatial maps from `Seurat::FindTransferAnchors()` and RCTD after affine registration with the NiftyReg algorithm.
4) Calculate match rates between Seurat predictions and ground truth supplied by 10X.

(Bonus) Visually compare Seurat's UMAP clustering vs the ground truth cell type IDs provided by 10X genomics.

## Questions: 
- Do the spatial maps from this TME dataset visually differ when produced with Seurat’s integration-based deconvolution instead of RCTD’s probabilistic deconvolution? 
- What is the increase in NMI achieved by optimally aligning the Seurat map to the RCTD map using affine transformations? 
- How close are Seurat's cell type predictions to the annotation provided by 10X?

## Links to project inspiration:

Use this link to see the original workflow the OSTA ebook authors wrote, that I adapted this script from

https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html

To skip to the images
(Note: Deconvolution result images are flipped vertically in this link)

https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#visualization

## Limitations:
- This version only calculates a confusion matrix on the Seurat results.
- The sample is one spatial transcriptomics slide/slice.
- No analysis is made into why one prediction algorithm captures certain oncological info better than the other.
- Did not find markers for clusters.

## Conclusion ####

In a sentence:
The breast cancer cells were difficult to identify, even using reference-based deconvolution.

### What evidence went into this?
#### 1) NMI Similarity (Seurat vs RCTD)

- Both methods are making similar cell type assignments
- Neither method is producing radically different spatial patterns
 
#### 2) Low Match Rates in tumour or preinvasive DCIS cells (Predictions vs Ground Truth)
DCIS1: 54.9% correct

DCIS2: 60.0% correct

Tumor: 60.9% correct

Stromal: 99.3% correct (not breast cancer, nor hard to identify)

### In a paragraph:
The moderate match rates between deconvolution predictions and ground truth annotations, combined with the high similarity between Seurat and RCTD predictions, suggest that the spatial transcriptomics deconvolution of the breast TME is fundamentally challenging due to the intermingled nature of cell populations. Both methods produced similar spatial patterns of predicted cell types, indicating that the observed limitations reflect biological complexity rather than methodological inadequacy.

### TME observations:
- DCIS1 and DCIS2 spots do not overlap - distinct populations.
- myoepithelial cells line the ducts and their spots overlap with both DCIS subtypes - DUCTAL carcinoma in situ (more in Seurat than RCTD image).
- DCIS2 and tumour regions do some bordering of each other.
- Stromal cells in the breast cancer microenvironment are non-cancerous cells that provide structural and biochemical support to the tumour, affecting its progression, metastasis, and therapy resistance. Stromal cells are concentrated around the tumour which makes sense because of the well-characterised desmoplastic response.
- Endothelial spots are scattered througout but concentrated around tumour because of angiogenesis. (clearer with Seurat?)
- Macrophages are recruited to stromal-DCIS interfaces - immune infiltration. Similar to Tumour-associated macrophages (TAMs) but at site of preinvasive population. Some TAMs seen too. (more in RCTD)

### Final evaluation:
Biologically meaningful spatial patterns that make oncological sense are still captured by the prediction algorithms.

And as always, remember the limitations.

Please create an issue if this workflow could be improved or made more robust.

## References:

Written at end of `TME.R` script, most called using `citation()`