# lso_EM_manuscript

synapsEM_analysis_macro_sm_v5.txt
-   Plasma membrane is annotated using freehand selection instead of freehand line to obtain area measurements
-   Synaptic vesicle is annotated using freehand selection to accommodate differences in vesicle shape. An ellipse is fit to each freehand selection.
-   “Endosome” is used for mitochondria.
-   “Active zone” is used for synaptic contact and “Dense Projection” is used for active zones
o   Note that not every synaptic profile analyzed had an active zone.
-   List of additional measurements added to the macro export text file:
o   Plasma membrane
§  Area (in pixels)
o   Synaptic vesicle
§  Radius – long radius of the ellipse
§  Major axis
§  Minor axis
§  Area (in pixels)
§  Aspect ratio (AR)
-   Added an roi.zip output in addition to each text file as a fail-safe measure
import_analysis_data.m
-   edits to process the additional measurements
start_analysis.m
-   commented out dp distance analysis
start_analysis_dp.m
-   original start_analysis file
min_2D_dist.m
-   edits to accommodate multiple “Dense Projections” measurements (some synaptic profiles analyzed had multiple active zones).
aggregate_syn_data.m
-   new file that sorts the raw data output from start_analysis.m 
- see comments in code for a list of variables
-   includes additional analyses to calculate vesicle density, synapse type based on aspect ratio, and nearest neighbors.
find_nn3.m
-   new file that finds the 3 nearest neighboring vesicles for each vesicle
 
Parameters:
-   pixel size = 0.8155 (nm/pixel)
-   bin size: 50
 
Procedure: raw_and_dist_data = start_analysis; agg_data = aggregate_syn_data(raw_and_dist_data);
 
Procedure (for distance from active zone): dist_data_dp = start_analysis_dp;
 
WARNINGS: 
-   MATLAB codes will only work with text files exported from synapsEM_analysis_macro_sm_edited_v5 due to the additional measurements that are outputted to the text file.
-   Text files with “Dense Projections” must be analyzed separately to obtain distance data to active zones (run start_analysis_dp instead of start_analysis)
 
