# ETTAA - CCW

The main workhorse function is `fixITB`, which takes the following arguments with default arguments based on the ETTAA application.

* `data`: A `data.frame` object,
* `GP`: Scalar denoting the grace period for the trial emulation (default: 365.25),
* `truncation.time`: Scalar denoting the truncation/administrative censoring time (default: 7 years),
* `bump`: Scalar, if intervention/event time occurs at time zero, small amount to add for analysis (default: 0.5),
* `stablised.weights`: Logical, should stabilised weights be used (default: `TRUE`),
* `rmsts`: Numeric vector containing the RMSTs of interest (default: 1,3,5,7 years),
* `rmst.scale`: Scalar, for presentation purposes what should the RMST be scaled by? (Default: 1 i.e. ignored),
* `surgery.formula`: Character, right-hand-side of the formula for the surgery IPTW/IPCW models (default: Combination of the key confounders in ETTAA).
* `control.formula`: Character, right-hand-side of the formula for the control IPTW/IPCW models (default: `NULL`, which copies `surgery.formula`).
* `max.weight`: Scalar, ad-hoc maximal weight for analysis (default: `NULL` which is no maximum weight),
* `verbose`: Logical, should statements be printed showing progress of `fixITB`? (Default: `FALSE`),
* `outputs`: Character vector, must be one (or both) of `"composite", "CR"`, tells `fixITB` what ATE to compute, defaults to both.

A set of dummy data is supplied to demonstrate this in the file `demo.R`. Here, S3 methods for `fixITB` are presented (printing/plotting). 

The file `explanatoryPlot.R` creates Figure 1.
