Head Model Comparison (Initial vs Final)

Overview
This project compares two versions of head models (initial vs final) and the
resulting FEM electric fields from electrode simulations. The workflow is
designed to quantify voxel-level segmentation changes and relate them to
field magnitude/direction changes. Data is sourced from both Beth Israel (BI)
and TIME, but the comparison is always within a subject (initial vs final),
not BI vs TIME.

Core concepts
- "Initial" vs "final" refers to two different head model runs for the same
  subject (e.g., run1 vs run2 or dilated vs original masks).
- Voxel masks are compared per tissue (wm, gm, csf, bone, skin, eyes, air).
- Electric field comparisons are computed at tetrahedral element centers.
- ROI shells/spheres around electrodes are used to link tissue change and
  field change by distance.

Key entry points
- analyze_heads.m: batch runner for one or more subjects.
- analyze_head.m: per-subject analysis of masks and (optionally) electric
  fields; builds ROI and segmentation summaries.

Script details

analyze_heads.m
- Purpose: orchestrates batch analysis across subjects for BI or TIME (each
  subject compares initial vs final within that dataset).
- Inputs (hard-coded):
  - TIME flag switches between TIME and BI data roots.
  - just_masks=1 runs voxel-only analysis (no electric field inputs).
  - is_DE_run=1 uses dilated/eroded masks for "final" when enabled.
- Workflow:
  - Builds paths for initial/final mask directories and FEM outputs.
  - Calls analyze_head() for each subject.
  - When just_masks=1, generates axial/coronal figures and saves per-subject
    Segmentation.mat in a *_VoxCount folder.
  - When just_masks=0, calls brain_blocking() for additional field analysis.
- Outputs:
  - *_VoxCount folders with segmentation summaries and slice images.
  - If full run, BB_* directories from brain_blocking.

analyze_head.m
- Purpose: top-level per-head analysis.
- Inputs:
  - electrodes: array of electrode positions; if empty, skips E-field logic.
  - ROIcenter: used for voxel ROI scaling/centering.
  - geometrySTR: geometry .mat/.geo path.
  - pre/post: .mat paths containing E_vector (pre and post fields).
  - mask_prepI/mask_prepF: paths to mask_prep folders for initial/final.
  - ROI_radii: vector of ROI radii (voxels).
  - name: subject identifier.
- Behavior:
  - If electrodes is empty, only voxel diffs are computed.
  - Otherwise:
    - Loads geometry, computes element centers.
    - Loads E_vector for pre/post, computes field magnitude and angle change.
    - Loads voxel segmentations with get_voxel_data().
    - For each ROI radius:
      - Converts voxel masks into xyz points (scale_and_center).
      - Builds ROI shells/spheres around electrode centers.
      - Averages field change in each ROI.
    - Calls Visualize(...) for plotting (see note below).
- Outputs:
  - EAngle, EMagnitude arrays for field differences.
  - ROI struct containing voxel and field summaries by radius.
  - Segmentation struct with initial/final/diff masks and voxel counts.

scale_and_center.m
- Purpose: converts voxel masks to xyz coordinates, scales to geometry space,
  and extracts ROI spheres/shells.
- Inputs: initial/final/diff masks, geometry, ROI center and radii.
- Outputs:
  - segmentationI/F/D: xyz point clouds per tissue.
  - segmentationDDeletion/Add: separate deletion/addition point clouds.
  - .shell fields for ROI shells (radOld..radNew).

visualize_head.m
- Purpose: generates plots linking voxel changes, conductivity changes, and
  E-field changes across ROI shells/spheres.
- Outputs:
  - Visualize_* folders with PNGs and figures for:
    - 3D centering proof: voxel change points and field vectors in geometry.
    - Axial slices by tissue index for initial/final/diff masks.
    - Axial slices by conductivity-weighted masks and diff.
    - Voxel change volume vs ROI radius (sphere and shell variants).
    - Voxel additions vs radius and voxel deletions vs radius.
    - Field magnitude change vs ROI radius (overall and by tissue).
    - Conductivity change vs field change scatter plots for ROI spheres.
    - Conductivity change summaries over distance (additions, deletions,
      total, and per-voxel averages).
    - Dice coefficient analyses:
      - Unweighted and weighted Dice vs ROI radius (spheres and shells).
      - Dice vs field change scatter plots (per tissue).
      - Weighted vs unweighted Dice correlation comparisons.
  - Plot filename map (saved under Visualize_*):
    - 3D centering proof: `3D_centering_proof.png`
    - Axial tissue slices: `Axial_Slices<name>.png`
    - Distance vs voxel change: `DistVsVoxChange<name>.png`
    - Distance vs voxel additions: `DistVsVoxAdd<name>.png`
    - Distance vs voxel deletions: `DistVsVoxDel<name>.png`
    - Conductivity change vs field change (spheres): `Cond_vs_Field_Change<name>.png`
    - ROI shell vox/field/radius summary: `ROI_Shell_Vox_Field_Radii_Plots<name>.png`
    - ROI shell additions summary: `ROI_Shell_Vox_Add_Field_Radii_Plots<name>.png`
    - ROI shell deletions summary: `ROI_Shell_Vox_Del_Field_Radii_Plots<name>.png`
    - Avg conductivity additions vs field: `Vox_Tot_Add_Cond_vs_Field<name>.png`
    - Avg conductivity additions vs radius: `Radius_vs_Vox_Tot_Add_Cond<name>.png`
    - Avg conductivity deletions vs field: `Vox_Tot_Del_Cond_vs_Field<name>.png`
    - Avg conductivity deletions vs radius: `Radius_vs_Vox_Tot_Del_Cond<name>.png`
    - Avg conductivity changes vs field: `Vox_Tot_Change_Cond_vs_Field<name>.png`
    - Avg conductivity changes vs radius: `Radius_vs_Vox_Change_Cond_Tot<name>.png`
    - Per-voxel conductivity additions vs field: `Vox_Add_Cond_vs_Field<name>.png`
    - Per-voxel conductivity additions vs radius: `Radius_vs_Vox_Add_Cond<name>.png`
    - Per-voxel conductivity deletions vs field: `Vox_Del_Cond_vs_Field<name>.png`
    - Per-voxel conductivity deletions vs radius: `Radius_vs_Vox_Del_Cond<name>.png`
    - Per-voxel conductivity changes vs field: `Vox_Change_Cond_vs_Field<name>.png`
    - Per-voxel conductivity changes vs radius: `Radius_vs_Vox_Change_Cond<name>.png`
    - Dice vs radius (spheres): `Radius_vs_Dice_Sphere<name>.png`
    - Dice vs field (spheres): `Dice_Spheres_vs_Field_Mag_Change<name>.png`
    - Weighted Dice vs field (spheres): `Weighted_Dice_Spheres_vs_Field_Mag_Change<name>.png`
    - Unweighted vs weighted Dice correlation: `UnWeighted_Dice_vs_Weighted_Sphere<name>.png`
    - Dice vs radius (shells): `Radius_vs_Dice_Shell<name>.png`

brain_blocking.m
- Purpose: ROI-based analysis in brain tissues (wm/gm) to compare field and
  conductivity changes along "cylinders" between electrodes and ROIs.
- Outputs:
  - BB_* folders with plots comparing conductivity vs field changes by ROI.

combine_plots.m
- Purpose: aggregates per-subject BB_* plots into combined scatter plots.
- Inputs: expects BB_* directories and subject list.

vox_aggregate_counts.m
- Purpose: aggregates voxel change counts across *_VoxCount folders.
- Outputs:
  - VoxCountAgg folder with BI and TIME bar plots, plus VoxCounts.mat.

dilate_final.m
- Purpose: generates dilated/eroded mask_prep variants for BI and TIME runs.
- Uses nifti_load/nifti_save and thickenbinvol/thinbinvol to modify masks.

Data expectations
- mask_prep folders contain tissue NIfTI masks (one per tissue).
- geometry files contain:
  - Geometry.node, Geometry.cell, and tissue indices (CI).
  - E_vector fields are stored in .mat files (pre/post).

External dependencies and missing functions
- get_voxel_data: used in analyze_head.m to load voxel masks; not defined in
  this directory.
- Visualize: analyze_head.m calls Visualize(...), but the local file is
  visualize_head.m. If needed, either rename or update the call.
- nifti_load/nifti_save, thickenbinvol/thinbinvol, elemvolume: expected to be
  available on the MATLAB path (commonly from SimNIBS or related toolboxes).

Typical usage
1) For mask-only comparison:
   - Set just_masks=1 in analyze_heads.m.
   - Run analyze_heads.m; review *_VoxCount outputs.

2) For full field + mask comparison:
   - Set just_masks=0 and ensure pre/post E_vector and geometry files exist.
   - Run analyze_heads.m; review Visualize_* and BB_* outputs.

Notes
- ROI_radii is currently logspace(0,250,20) in analyze_heads.m; confirm this
  is intended (values are in voxels).
- Paths are hard-coded to /Volumes/spiral/...; adjust for your environment.
