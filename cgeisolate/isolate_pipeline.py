import os
import sys
import subprocess
import csv

from cgeisolate import kma


def isolate_pipeline(args):
    print("Starting the isolate pipeline...")

    if args.db_dir is None:
        if not os.path.exists('/var/lib/cge/database/cge_db/'):
            sys.exit('Please install the cge_db. It should be located in /var/lib/cge/database/cge_db')
        else:
            args.db_dir = '/var/lib/cge/database/cge_db/'
            print(f"Using CGE database directory: {args.db_dir}")

    if args.output is None:
        args.output = '/var/lib/cge/results/{}'.format(args.name)

    print(f"Creating output directory: {args.output}")
    os.system('mkdir -p ' + args.output)

    print(f"Running KMA for bacteria alignment on input: {args.input}")
    os.system('{}kma -t_db {} -i {} -o {} -ID 75 -md 5 -ont -1t1 -mem_mode -t 8'\
              .format(args.env_path_bin, args.db_dir + '/bac_db/bac_db', args.input, args.output + "/bacteria_alignment"))
    #kma.KMARunner(args.input,
    #              args.output + "/bacteria_alignment",
    #              args.db_dir + '/bac_db/bac_db',
    #              "-ID 75 -md 5 -ont -1t1 -mem_mode -t 8").run()

    highest_scoring_hit = get_highest_scoring_hit_template(args.output + "/bacteria_alignment.res")
    print(f"Highest scoring hit: {highest_scoring_hit}")

    if 'Escherichia coli' in highest_scoring_hit:
        print("Escherichia coli detected. Running KMA for virulence...")
        os.system('{}kma -t_db {} -i {} -o {} -ont -md 5'\
                  .format(args.env_path_bin, args.db_dir + '/virulence_db/virulence_db', args.input, args.output + "/virulence"))

    print("Running KMA for AMR...")
    os.system('{}kma -t_db {} -i {} -o {} -ont -md 5 -mem_mode -t 8'\
              .format(args.env_path_bin, args.db_dir + '/resfinder_db/resfinder_db', args.input, args.output + "/amr"))

    print("Running KMA for plasmid...")
    os.system('{}kma -t_db {} -i {} -o {} -ont -md 5 -mem_mode -t 8'\
                .format(args.env_path_bin, args.db_dir + '/plasmid_db/plasmid_db', args.input, args.output + "/plasmid"))

    print("Running MLST...")
    os.system('{}kgt_mlst -i {} -o {} -db_dir {} -md 5'\
                .format(args.env_path_bin, args.input, args.output + "/mlst", args.db_dir))

    print("Creating report...")
    report = create_report(args, highest_scoring_hit)
    with open(args.output + "/report.txt", 'w') as file:
        file.write(report)

    print("Pipeline completed successfully. Report generated and stored in: " + args.output + "/report.txt")
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
