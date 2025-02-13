"""
Configuration file for QC thresholds and parameters
"""

# Quality assessment thresholds
QC_THRESHOLDS = {
    "GOOD": {
        "min_quality": 15,
        "min_coverage": 50,
        "min_read_length": 3000,
    },
    "FAIR": {
        "min_quality": 12,
        "min_coverage": 20,
        "min_read_length": 2000,
    },
    "POOR": {
        "min_quality": 12,  # Less than this is poor
        "min_coverage": 20,  # Less than this is poor
        "min_read_length": 1000,  # Less than this is unusually short for ONT
    },
    "GC_CONTENT": {
        "min": 25, # bacterial genome range: 25 - 75%
        "max": 75
    }
}

# Default trim parameters
TRIM_DEFAULTS = {
    "min_length": 16,
    "max_length": 2147483647,
    "min_phred": 20,
    "min_internal_phred": 0,
    "min_average_quality": 10,  # -eq parameter
    "trim_5_prime": 0,
    "trim_3_prime": 0
}