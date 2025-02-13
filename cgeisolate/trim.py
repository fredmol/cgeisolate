import os
import logging
import subprocess
import sys

from cgeisolate.qc_report import create_qc_report
from cgeisolate.qc_config import TRIM_DEFAULTS

class TrimRunner():
    def __init__(self, input_file, output_dir, name):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.input_file = input_file
        self.output_dir = output_dir
        self.trimmed_name = f"{name}_trimmed"
        self.trimmed_output_path = os.path.join(self.output_dir, self.trimmed_name)
        self.parameters = TRIM_DEFAULTS

    def run(self):
        """Runs KMA trim on input file"""
        # Make sure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
        
        trim_cmd = (
            f"kma trim "
            f"-ml {self.parameters['min_length']} "
            f"-xl {self.parameters['max_length']} "
            f"-mp {self.parameters['min_phred']} "
            f"-mi {self.parameters['min_internal_phred']} "
            f"-eq {self.parameters['min_average_quality']} "
            f"-5p {self.parameters['trim_5_prime']} "
            f"-3p {self.parameters['trim_3_prime']} "
            f"-qc "
            f"-i {self.input_file} "
            f"-o {self.trimmed_output_path}"
        )
        
        self.logger.info(f"Running KMA trim with command: {trim_cmd}")
        
        ret = os.system(trim_cmd)
        if ret != 0:
            sys.exit(f"Error: KMA trim failed with return code {ret}")
            
        expected_output = f"{self.trimmed_output_path}.fq"
        if not os.path.exists(expected_output):
            sys.exit(f"Error: KMA trim failed to create output file {expected_output}")
            
        # Generate QC report if json exists
        json_output = f"{self.trimmed_output_path}.json"
        if os.path.exists(json_output):
            try:
                qc_report_path = create_qc_report(
                    json_output, 
                    self.output_dir, 
                    self.trimmed_name, 
                    self.parameters  # Pass the parameters here
                )
                self.logger.info(f"Generated QC report: {qc_report_path}")
            except Exception as e:
                self.logger.error(f"Failed to generate QC report: {str(e)}")
        
        return expected_output