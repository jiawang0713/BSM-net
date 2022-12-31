# BSM-net
Code for the project: *Variable Selection for High-dimensional Nodal Attributes in Social Networks with Degree Heterogeneity*. Please see the folder structure below. 

---------------------------------------------------------------------------------------
## simulation

### proposed_method

* `net_nonsp_worker.jl`: the proposed method (BSM-net) written in Julia 1.3.1.

* `net_nonsp_worker.PBS`: A PBS script used to submit jobs to the scheduler

* `net_nonsp_mgr.bashrc`: A bash shell script file used to submit jobs to the scheduler

* `net_sp_worker.jl`: the proposed method (BSM-net.sp) written in Julia 1.3.1.

* `net_sp_worker.PBS`: A PBS script used to submit jobs to the scheduler

* `net_sp_mgr.bashrc`: A bash shell script file used to submit jobs to the scheduler

### existing_method

* `existing_method_pall.R`: existing methods (LASSO/SCAD/MCP) written in R 3.6.1
    

**NOTE**: The data sets are firstly generated in `net_sp_worker.jl`, then called in `net_nonsp_worker.jl` and `existing_method_LLA.R`.

-------------------------------------------------------------------------------------------

## real_data

### facebook_net

#### ----- codes

* `rd_sp_fb.jl`: the proposed method (BSM-net.sp) on Facebook friendship data set written in Julia 1.3.1.

* `fb_net_plot.R`: R codes for rendering results and plots.

#### ----- data and results

* `Y.csv`: the adjacency matrix Y (constructed based on https://snap.stanford.edu/data/ego-Facebook.html).

* `X.csv`: the node-level features X (constructed based on https://snap.stanford.edu/data/ego-Facebook.html).

* `21featnames`: anonymous feature names for those 21 survey questions about users profiles (downloaded from https://snap.stanford.edu/data/ego-Facebook.html).

 ---------------------------------------------------------------------------------------------------------------------------------
 
### citation_net

#### ----- codes

* `rd_sp_citation.jl`: the proposed method (BSM-net.sp) on paper citation data set written in Julia 1.3.1.

* `create_X.R`: R codes for LDA on paper abstract (downloaded from https://www.aminer.cn/citation) to construct node-features X.

* `create_Y.R`: R codes to construct the adjacency matrix Y based on citation data set https://www.aminer.cn/citation.

#### ----- data and results

* `citation_Y_1207_n=476.csv`: the adjacency matrix Y (constructed based on https://www.aminer.cn/citation).

* `citation_X_1207_n=476.csv`: the node-level features X (constructed based on LDA on paper abstract https://www.aminer.cn/citation).

* `Terms_1207_n=476.csv`: the 5 most frequent terms for each topics. (created in `create_X.R`)

* `selected_indx_1207_n=476.csv`: indices of selected papers. (created in `create_Y.R`)

 ---------------------------------------------------------------------------------------------------------------------------------

### trade_net

#### ----- codes

* `preprocessing.R`: R codes to construct the node-level features X based on WDIs and the adjacency matrix Y based on trade values.

**NOTE**: .jl file is skipped here as it is almost the same as `rd_sp_citation.jl`

#### ----- data and results

* `trade_Y_n=142.csv`: the adjacency matrix Y (constructed based on trade data set http://www.cepii.fr/CEPII/en/bdd_modele/presentation.asp?id=37).

* `trade_X_n=142.csv`: the node-level features X (constructed based on WDI data set https://datacatalog.worldbank.org/dataset/world-development-indicators).

* `final_country_codes_n=142.csv`: the collection of n = 142 selected countries. (created in `preprocessing.R`)

* `Y_2017.csv`: aggregated total value of trade in cash for edge i->j in the year 2017. (created in `preprocessing.R`)

* `country_codes_V202001.csv`: country code for trade data set downloaded from http://www.cepii.fr/CEPII/en/bdd_modele/presentation.asp?id=37.

* `WDI.xlsx`: WDI data set downloaded from https://datacatalog.worldbank.org/dataset/world-development-indicators.
