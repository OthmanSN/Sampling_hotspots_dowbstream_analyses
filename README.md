# Sampling_hotspots_downstream_analyses
Prediction of Dicroglossidae sampling blind spots across Asian ecozones, analysed on best-supported speices delimitation Model A1

Overview
This directory contains downstream analyses, summary tables, and figure outputs used to examine geographic patterns of taxonomic mismatch and sampling coverage across Dicroglossidae datasets. The analyses were based on the best-supported species delimitation framework, Model A1, and were designed to identify regional inconsistencies between GenBank-based labels and revised assignments, as well as areas of low sampling density that may require further genomic investigation. The workflow integrates cleaned assignment tables, regional summaries, geographic coordinates, and graphical outputs including alluvial plots, heatmaps, bar charts, and regional bubble maps.

Main input files
dicroglossid_modelA1_vs_genbank_full_georef.csv

Main cleaned dataset used for downstream analyses. This file includes:
sample or sequence identifier
GenBank-based label
Model A1 assignment
match or mismatch status
locality and country information
region grouping
latitude and longitude coordinates

This file is the primary input for regional summaries, mapping, and mismatch analyses.

dicroglossid_alluvial_input.csv 
Input table used for generating alluvial plots comparing GenBank-based assignments and revised Model A1 assignments.

dicroglossid_pair_summary.csv
Summarised pairwise counts of GenBank labels versus Model A1 assignments. Used for pair heatmaps and contingency testing.

dicroglossid_region_summary.csv
Regional summary table used in earlier plotting steps. Final region-based outputs were regenerated directly from the updated georeferenced main dataset.

Main R workflow
Sampling_hotspots_downstream_analyses.R

R script used to:
-Clean and standardise taxon labels
-harmonise country, locality, and region fields
-generate species- and clade-level alluvial plots
-generate country and region bar plots
-generate region and pair heatmaps
-test whether reassignment patterns were structured using contingency tables and chi-square analysis
-build regional bubble maps of mismatch proportion
-export publication-ready figures and supplementary tables
-Main output figures

-Alluvial plots
Fig_main_species_alluvial_revise_label.pdf/png
Species-level alluvial plot showing reassignment from GenBank labels to Model A1 assignments.

-FigS_clade_alluvial_revised_nolabel.pdf/png
Clade-level alluvial plot used as supplementary material.

-Regional summaries
Fig_country_barplot_updated.pdf/png
Country-level bar plot summarising match and mismatch counts.
Fig_region_barplot_updated.pdf/png

-Region-level bar plot summarising match and mismatch counts.
Fig_region_heatmap_updated.pdf/png

-Heatmap showing counts of Model A1 assignments across region groups.
Pairwise reassignment summaries
Fig_pair_heatmap.pdf/png

-Heatmap of GenBank labels versus Model A1 assignments.
Fig_pair_heatmap_mismatch.pdf/png
Mismatch-only heatmap highlighting reassignment structure.

-Geographic outputs
Fig_region_bubble_*.pdf/png
Bubble maps summarising regional mismatch proportion.

-In the final version:
Black points represent all georeferenced sequence entries;
bubble size represents total records per region;
bubble fill represents mismatch proportion;
bubble positions are based on median sampling coordinates within each region, with minor manual adjustment only where needed for readability.

Region definitions

Final region groupings used in the analyses include:

Hainan (insular)
Taiwan (insular)
Japan (archipelago)
Lower Yangtze region
Eastern China mainland
Mainland Southeast Asia
Sundaic (Maritime Southeast Asia)
Other / unclear

These groupings were used for descriptive comparison of mismatch patterns and for identifying sampling blind spots across ecozones.

Interpretation of key outputs

The downstream analyses were designed to address three main questions:

How often do GenBank labels agree with the best-supported Model A1 assignment?
Are mismatches distributed evenly across regions, or concentrated in specific parts of Asia?
Which regions combine low sampling density with relatively high mismatch, suggesting priorities for future genomic sampling?

The combined figure panels and summary tables indicated that mismatch was highest in the Lower Yangtze region and several insular regions, while lower mismatch was observed in parts of Mainland Southeast Asia. Regions with limited sampling but substantial mismatch were interpreted as likely blind spots for future data collection.

Notes
Final analyses should be run from the georeferenced dataset rather than earlier cleaned versions without complete coordinates.
Region summaries and bubble maps should always be regenerated if region labels, country assignments, or coordinates are updated.
Supplementary taxonomic interpretation, including type material and detailed lineage notes, should be kept separate from these downstream descriptive outputs unless explicitly required.
Suggested citation within manuscript text

These outputs correspond to the regional mismatch and sampling coverage analyses associated with Fig. 5 and Table S16 of the manuscript.
