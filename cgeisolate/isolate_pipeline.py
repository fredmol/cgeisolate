import os
import sys
import subprocess
import csv

from cgeisolate import kma

def isolate_pipeline(args):
    if args.folder is not None:
        if args.name is None:
            sys.exit('Please provide a name for the merged file')
        else:
            merge_fastq_files(args.folder)
            args.input = os.path.join(os.path.expanduser('~'), args.name)
            #args.output = os.path.join(os.path.expanduser('~'), args.name)

    if args.db_dir is None:
        if not os.path.exists('/opt/cge/cge_db'):
            sys.exit('Please install the cge_db. It should be located in /opt/cge/cge_db')
        else:
            args.db_dir = '/opt/cge/cge_db'

    os.system('mkdir -p ' + args.output)

    # Run KMA for bacteria alignment
    kma.KMARunner(args.input,
                  args.output + "/bacteria_alignment",
                  args.db_dir + '/bac_db/bac_db',
                  "-ID 75 -md 5 -ont -1t1 -mem_mode -t 8").run()

    highest_scoring_hit = get_highest_scoring_hit_template(args.output + "/bacteria_alignment.res")

    if 'Escherichia coli' in highest_scoring_hit:
        kma.KMARunner(args.input,
                      args.output + "/virulence",
                      args.db_dir + '/virulence_db/virulence_db',
                      "-ont -md 5").run()

    # Run KMA for AMR
    kma.KMARunner(args.input,
                  args.output + "/amr",
                  args.db_dir + '/resfinder_db/resfinder_db',
                  "-ont -md 5 -mem_mode -t 8").run()

    # Run KMA for plasmid
    kma.KMARunner(args.input,
                  args.output + "/plasmid",
                  args.db_dir + '/plasmid_db/plasmid_db',
                  "-ont -md 5 -mem_mode -t 8").run()

    # Run MLST
    cmd = 'kgt_mlst -i {} -o {} -db_dir {} -md 5'\
        .format(args.input, args.output + "/mlst", args.db_dir)
    os.system(cmd)

    report = create_report(args, highest_scoring_hit)
    with open(args.output + "/report.txt", 'w') as file:
        file.write(report)

    return 'isolate_pipeline'

def merge_fastq_files(source_directory, output_name):
    """
    Merge all fastq.gz files in the given directory and save the output with the specified name in the home directory.

    Args:
    source_directory (str): Path to the directory containing fastq.gz files.
    output_name (str): Name for the output file.
    """
    # Home directory path
    home_directory = os.path.expanduser('~')

    # Output file path with the specified name
    output_file = os.path.join(home_directory, f'{output_name}.fastq.gz')

    # Get a list of all fastq.gz files in the source directory
    fastq_files = [f for f in os.listdir(source_directory) if f.endswith('.fastq.gz')]

    # Open the output file in write mode
    with gzip.open(output_file, 'wb') as f_out:
        # Iterate over each file and append its content to the output file
        for file in fastq_files:
            file_path = os.path.join(source_directory, file)
            with gzip.open(file_path, 'rb') as f_in:
                # Copy the content of each file to the output file
                f_out.writelines(f_in)

    print(f"All files merged into {output_file}")


def get_highest_scoring_hit_template(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        highest_scoring_hit = None
        max_score = float('-inf')  # Initialize with the smallest possible number

        for row in reader:
            try:
                score = float(row['Score'])
                if score > max_score:
                    highest_scoring_hit = row
                    max_score = score
            except ValueError:
                # Handle the case where the score is not a number
                continue

    return highest_scoring_hit['#Template'] if highest_scoring_hit else None

def read_tab_separated_file(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        return list(reader)

def format_results_section(results, section_title):
    report = f"{section_title}:\n"
    report += "-" * 60 + "\n"
    for result in results:
        report += f"Template: {result['#Template']}\n"
        report += f"Identity: {result['Template_Identity'].strip()}, Coverage: {result['Template_Coverage'].strip()}, Depth: {result['Depth'].strip()}\n\n"
    return report

def create_report(args, highest_scoring_hit):
    # Load results
    gene_data = read_tab_separated_file(args.db_dir + '/phenotypes.txt')

    bacteria_results = read_tab_separated_file(args.output + "/bacteria_alignment.res")
    amr_results = read_tab_separated_file(args.output + "/amr.res")
    plasmid_results = read_tab_separated_file(args.output + "/plasmid.res")
    highest_scoring_hit = get_highest_scoring_hit_template(args.output + "/bacteria_alignment.res")

    phenotypes = set()
    for amr_result in amr_results:
        for gene in gene_data:
            if gene['Gene_accession no.'] == amr_result.get('#Template'):
                phenotypes.update(gene['Phenotype'].split(','))

    highest_scoring_hit_details = get_highest_scoring_hit_details(args.output + "/bacteria_alignment.res")

    report = "Pipeline Results Report\n"
    report += "=" * 60 + "\n"
    if highest_scoring_hit_details:
        report += f"Template: {highest_scoring_hit_details['#Template']}\n"
        report += f"Identity: {highest_scoring_hit_details['Template_Identity']}, "
        report += f"Coverage: {highest_scoring_hit_details['Template_Coverage']}, "
        report += f"Depth: {highest_scoring_hit_details['Depth']}\n\n"
    else:
        report += "No bacteria alignment hits found.\n\n"

    # AMR Results Section
    report += format_results_section(amr_results, "Antimicrobial Resistance (AMR) Findings")

    # Plasmid Results Section
    report += format_results_section(plasmid_results, "Plasmid Findings")

    # Expected Phenotypes Based on AMR Genes Section
    report += "Expected Phenotypes Based on AMR Genes:\n"
    report += "-" * 60 + "\n"
    if phenotypes:
        for phenotype in sorted(phenotypes):
            report += f"• {phenotype.strip()}\n"
    else:
        report += "No phenotypes expected based on AMR genes.\n"
    report += "\n"

    if 'Escherichia coli' in highest_scoring_hit:
        virulence_results = read_tab_separated_file(args.output + '/virulence.res')
        report += "Virulence Factors for Escherichia coli:\n"
        report += "-" * 60 + "\n"
        for result in virulence_results:
            report += f"Template: {result['#Template']}\n"
            report += f"Identity: {result['Template_Identity'].strip()}, Coverage: {result['Template_Coverage'].strip()}, Depth: {result['Depth'].strip()}\n\n"
    else:
        report += "No virulence factors analysis for Escherichia coli.\n"


    return report


def get_highest_scoring_hit_details(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        highest_scoring_hit = None
        max_score = float('-inf')

        for row in reader:
            try:
                score = float(row['Score'])
                if score > max_score:
                    highest_scoring_hit = row
                    max_score = score
            except ValueError:
                continue

    return highest_scoring_hit
