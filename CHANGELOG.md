# Changelog

## [1.5.1] - 2025-02-28
### Changed
- There is no longer a check to see if output folder exists. The tool will run regardless of it existing or not.

## [1.5.0] - 2025-02-20
### Changed
- **QC/Trimming Removal:** Removed the integrated quality control and read trimming functionality from the pipeline. All QC operations are now centralized in the separate CGEqc tool.

### Fixed
- Minor bug fixes and improved error messages in downstream analysis steps.

## [1.4.0-alpha] - 2025-02-13
### ⚠️ Important
- **KMA Dependency:** This version requires KMA from the `cgelabs` branch (alpha version)
  - https://bitbucket.org/genomicepidemiology/kma/branch/cgelabs
  - Features like QC reporting and trimming will not work with the main KMA branch

### Added
- **Quality Control Features** (requires KMA cgelabs branch):
  - Automatic QC report generation for trimmed reads
  - Read quality and length distribution visualization
  - Three-tier quality grading system (GOOD/FAIR/POOR)
  - Configurable quality thresholds via central configuration file
  - Quality assessment metrics including sequencing depth, read quality, and GC content
- **Read Trimming** (requires KMA cgelabs branch):
  - Integrated KMA trim functionality with configurable parameters
  - Detailed trim statistics and parameter reporting

### Changed
- **Quality Assessment:**
  - Calibrated thresholds specifically for ONT bacterial sequencing
  - Expanded GC content range to 25-75% to reflect bacterial genome diversity
  - Enhanced depth of coverage calculation with clear methodology

### Fixed
- **User Experience:**
  - Improved feedback clarity in quality assessment reports
  - Added explanatory notes for technical metrics

## [1.3.0] – 2025-02-13
### Added
- **PDF Reporting:** Integrated PDF report generation using WeasyPrint and Jinja2 with custom HTML templates and CSS (assets and templates directories).
- **Styline:** Custom DTU styling for reports
- **AMR Visualization:** Implemented a drug class distribution plot for AMR genes using matplotlib.

### Changed
- **Pipeline Enhancements:** Updated `isolate_pipeline.py` to include robust error handling, and PDF report generation. 
- **Packaging:** Modified `setup.py` to include new dependencies (WeasyPrint, Jinja2, matplotlib) and package data (assets and templates). The Conda recipe in `make_conda.py` now reflects these new dependencies.
- **Version Bump:** Version updated from 1.2.0 to 1.3.0.

### Fixed
- **User Feedback:** Enhanced logging and progress messages to better inform users of pipeline steps and any issues encountered.

## [1.2.0] – Initial Release
- Basic bacterial isolate analysis pipeline featuring KMA alignment and text-based report generation.
- Core modules included for running KMA (`kma.py`) and parsing alignment results.
- Provided a simple command-line interface for running the analysis.
