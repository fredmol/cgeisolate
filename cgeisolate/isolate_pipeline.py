import os
import sys
import subprocess
import csv

from cgeisolate import kma

def isolate_pipeline(args):
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

    report = create_report(args)
    print (report)

    return 'isolate_pipeline'

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

def create_report(args):
    # Load results
    bacterial_results = read_tab_separated_file(args.output + "/bacteria_alignment.res")
    amr_results = read_tab_separated_file(args.output + "/amr.res")
    plasmid_results = read_tab_separated_file(args.output + "/plasmid.res")
    highest_scoring_hit = get_highest_scoring_hit_template(args.output + "/bacteria_alignment.res")

    report = "Pipeline Results Report\n"
    report += "=" * 60 + "\n"

    # Top Scoring Bacterial Alignment Hit
    report += "Top Scoring Bacterial Alignment Hit:\n"
    report += "-" * 60 + "\n"
    report += f"Template: {highest_scoring_hit}\n\n"

    # AMR Results Section
    report += format_results_section(amr_results, "Antimicrobial Resistance (AMR) Findings")

    # Plasmid Results Section
    report += format_results_section(plasmid_results, "Plasmid Findings")

    # Expected Phenotypes Based on AMR Genes Section (assuming this data is within amr_results)
    report += "Expected Phenotypes Based on AMR Genes:\n"
    report += "-" * 60 + "\n"
    phenotypes = set()
    for amr_result in amr_results:
        phenotypes.update(amr_result['Phenotype'].split(','))  # Assuming 'Phenotype' is a field in amr_results
    for phenotype in sorted(phenotypes):
        report += f"• {phenotype.strip()}\n"
    report += "\n"

    return report

