import os
import sys
import subprocess
import csv
import shutil

from cgeisolate import kma
from cgeisolate import version

# New imports - Frederik Møller
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML, CSS
from pathlib import Path
from datetime import datetime
from io import BytesIO
import matplotlib.pyplot as plt
from io import BytesIO
import base64


def isolate_pipeline(args):
    print("Starting the isolate pipeline...")

    # Check if output folder already exists
    output_dir = '/var/lib/cge/results/{}'.format(args.name)
    if os.path.exists(output_dir):
        sys.exit(
            f"Error: Output directory '{output_dir}' already exists. Please choose a different name or delete the existing directory.")


    if args.db_dir is None:
        if not os.path.exists('/var/lib/cge/database/cge_db'):
            sys.exit('Please install the cge_db. It should be located in /var/lib/cge/database/cge_db')
        else:
            args.db_dir = '/var/lib/cge/database/cge_db'
            print(f"Using CGE database directory: {args.db_dir}")

    if args.output is None:
        args.output = '/var/lib/cge/results/{}'.format(args.name)

    print(f"Creating output directory: {args.output}")
    os.system('mkdir -p ' + args.output)

    # Copy the input FASTQ file to the output directory
    print(f"Copying input FASTQ file to the output directory: {args.input} -> {args.output}")
    shutil.copy(args.input, args.output)

    print(f"Running KMA for bacteria alignment on input: {args.input}")
    os.system('kma -t_db {} -i {} -o {} -ID 75 -md 5 -ont -1t1 -mem_mode -t 8 -ef'\
              .format(args.db_dir + '/bac_db/bac_db', args.input, args.output + "/bacteria_alignment"))


    # Make sure there is actually a hit, so job doesn't crash
    highest_scoring_hit = get_highest_scoring_hit_template(args.output + "/bacteria_alignment.res")

    if highest_scoring_hit is None:
        print("No highest scoring hit found in bacteria alignment. Skipping virulence genes.")
    else:
        print(f"Highest scoring hit: {highest_scoring_hit}")
        if 'Escherichia coli' in highest_scoring_hit:
            print("Escherichia coli detected. Running KMA for virulence...")
            os.system('kma -t_db {} -i {} -o {} -ont -md 5'\
                      .format(args.db_dir + '/virulence_db/virulence_db', args.input, args.output + "/virulence"))

    print("Running KMA for AMR...")
    os.system('kma -t_db {} -i {} -o {} -ont -md 5 -mem_mode -t 8'\
              .format(args.db_dir + '/resfinder_db/resfinder_db', args.input, args.output + "/amr"))

    print("Running KMA for plasmid...")
    os.system('kma -t_db {} -i {} -o {} -ont -md 5 -mem_mode -t 8'\
                .format(args.db_dir + '/plasmid_db/plasmid_db', args.input, args.output + "/plasmid"))

    print("Running MLST...")
    os.system('kgt_mlst -i {} -o {} -db_dir {} -md 5'\
                .format(args.input, args.output + "/mlst", args.db_dir))

    print("Creating report...")
    report = create_report(args, highest_scoring_hit)
    with open(args.output + "/report.txt", 'w') as file:
        file.write(report)
    
    
    ############################ PDF REPORT ###################################
    
    # Load all results
    gene_data = read_tab_separated_file(args.db_dir + '/phenotypes.txt')
    bacteria_results = read_tab_separated_file(args.output + "/bacteria_alignment.res")
    amr_results = read_tab_separated_file(args.output + "/amr.res")
    plasmid_results = read_tab_separated_file(args.output + "/plasmid.res")
    file_stats = get_file_stats(args.input, args.output)  
    
    # Get phenotypes
    phenotypes = set()
    for amr_result in amr_results:
        for gene in gene_data:
            if gene['Gene_accession no.'] == amr_result.get('#Template'):
                phenotypes.update(gene['Phenotype'].split(','))

    # Create text report
    report = create_report(args, highest_scoring_hit)
    with open(args.output + "/report.txt", 'w') as file:
        file.write(report)

    # Create PDF report
    print("Creating PDF report...")
    pdf_path = create_pdf_report(args, highest_scoring_hit, gene_data, bacteria_results, 
                               amr_results, plasmid_results, phenotypes, file_stats)
    if pdf_path:
        print(f"PDF report generated and stored in {pdf_path}")
        
    ############################ PDF REPORT DONE ##############################

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

    report = "Analysis report: {}\n".format(args.name)
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
    
    try:
        if 'Escherichia coli' in highest_scoring_hit:
            virulence_results = read_tab_separated_file(args.output + '/virulence.res')
            report += "Virulence Factors for Escherichia coli:\n"
            report += "-" * 60 + "\n"
            for result in virulence_results:
                report += f"Template: {result['#Template']}\n"
                report += f"Identity: {result['Template_Identity'].strip()}, Coverage: {result['Template_Coverage'].strip()}, Depth: {result['Depth'].strip()}\n\n"
        else:
            report += "No virulence factors analysis for Escherichia coli.\n"
    except:
        report += "No highest hit could be reported.\n"

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


############################## PDF REPORT #####################################


def create_pdf_report(args, highest_scoring_hit, gene_data, bacteria_results, amr_results, plasmid_results, phenotypes, file_stats):
    """Generate a PDF report for the isolate analysis"""
    # Setup paths
    package_dir = Path(__file__).parent
    assets_dir = package_dir / 'assets'
    logo_path = assets_dir / 'dtu_logo.png'
    
    # Default template data
    template_data = {
        'name': getattr(args, 'name', 'Unknown'),
        'date': datetime.now().strftime("%B %d, %Y"),
        'version': getattr(version, '__version__', 'Unknown'),
        'logo_data_url': '',
        'file_stats': file_stats if file_stats else {
            'file_name': 'Unknown',
            'file_size': 'N/A',
            'read_count': 'N/A',
            'kma_version': 'N/A'
        },
        'species_name': 'No species identified',
        'species_details': {
            'name': 'No species identified',
            'identity': 'N/A',
            'coverage': 'N/A',
            'depth': 'N/A'
        },
        'amr_count': 0,
        'phenotype_count': 0,
        'plasmid_count': 0,
        'mlst_result': 'N/A',
        'virulence_count': 'N/A',
        'virulence_label': 'Virulence Genes (E. coli only)',
        'amr_genes': [],
        'plasmids': [],
        'virulence_genes': [],
        'is_ecoli': False,
        'phenotypes': [],
        'grouped_phenotypes': [],
        'drug_class_plot': None
    }
    
    try:
        # Load and encode logo if available
        if logo_path.exists():
            with open(logo_path, 'rb') as f:
                logo_data = f.read()
                template_data['logo_data_url'] = f'data:image/png;base64,{base64.b64encode(logo_data).decode("utf-8")}'
        
        # Get and format species details
        hit_details = get_highest_scoring_hit_details(args.output + "/bacteria_alignment.res") if hasattr(args, 'output') else None
        if hit_details:
            try:
                species_details = {
                    'name': hit_details['#Template'],
                    'identity': f"{float(hit_details['Template_Identity']):.2f}",
                    'coverage': f"{float(hit_details['Template_Coverage']):.2f}",
                    'depth': f"{float(hit_details['Depth']):.2f}"
                }
                template_data['species_details'] = species_details
                template_data['species_name'] = ' '.join(species_details['name'].split()[1:3])
            except (KeyError, ValueError):
                pass
        
        # Determine if E. coli
        is_ecoli = 'Escherichia coli' in template_data['species_details']['name']
        template_data['is_ecoli'] = is_ecoli
        
        # Group phenotypes by drug class if data available
        if amr_results and gene_data:
            class_phenotypes = {}
            for result in amr_results:
                gene_id = result.get('#Template', '')
                for gene in gene_data:
                    if gene.get('Gene_accession no.') == gene_id:
                        drug_classes = [cls.strip() for cls in gene.get('Class', '').split(',') if cls.strip()]
                        phenotypes_list = [p.strip() for p in gene.get('Phenotype', '').split(',') if p.strip()]
                        
                        for drug_class in drug_classes:
                            if drug_class not in class_phenotypes:
                                class_phenotypes[drug_class] = set()
                            class_phenotypes[drug_class].update(phenotypes_list)
            
            template_data['grouped_phenotypes'] = [
                {'class': cls, 'phenotypes': sorted(list(phens))}
                for cls, phens in sorted(class_phenotypes.items())
            ]

        # Format AMR genes if available
        if amr_results:
            formatted_amr_genes = []
            for result in amr_results:
                try:
                    identity = float(result.get('Template_Identity', '0').strip())
                    coverage = float(result.get('Template_Coverage', '0').strip())
                    depth = float(result.get('Depth', '0').strip())
                    
                    if identity >= 100.0 and coverage >= 100.0:
                        match_quality = 'perfect-match'
                    elif identity >= 100.0 or coverage >= 100.0:
                        match_quality = 'good-match'
                    else:
                        match_quality = 'partial-match'
                    
                    formatted_amr_genes.append({
                        'template': result.get('#Template', 'Unknown'),
                        'identity': f"{identity:.2f}",
                        'coverage': f"{coverage:.2f}",
                        'depth': f"{depth:.2f}",
                        'match_quality': match_quality
                    })
                except (ValueError, KeyError):
                    continue
            
            template_data['amr_genes'] = formatted_amr_genes
            template_data['amr_count'] = len(formatted_amr_genes)
        
        # Format virulence genes if E. coli
        if is_ecoli:
            try:
                virulence_results = read_tab_separated_file(args.output + "/virulence.res")
                formatted_virulence_genes = []
                for result in virulence_results:
                    try:
                        identity = float(result.get('Template_Identity', '0').strip())
                        coverage = float(result.get('Template_Coverage', '0').strip())
                        depth = float(result.get('Depth', '0').strip())
                        
                        if identity >= 100.0 and coverage >= 100.0:
                            match_quality = 'perfect-match'
                        elif identity >= 100.0 or coverage >= 100.0:
                            match_quality = 'good-match'
                        else:
                            match_quality = 'partial-match'
                        
                        formatted_virulence_genes.append({
                            'template': result.get('#Template', 'Unknown'),
                            'identity': f"{identity:.2f}",
                            'coverage': f"{coverage:.2f}",
                            'depth': f"{depth:.2f}",
                            'match_quality': match_quality
                        })
                    except (ValueError, KeyError):
                        continue
                
                template_data['virulence_genes'] = formatted_virulence_genes
                template_data['virulence_count'] = str(len(formatted_virulence_genes))
                template_data['virulence_label'] = "Virulence Genes"
            except Exception:
                pass
        
        # Format plasmids if available
        if plasmid_results:
            formatted_plasmids = []
            for result in plasmid_results:
                try:
                    formatted_plasmids.append({
                        'template': result.get('#Template', 'Unknown'),
                        'identity': f"{float(result.get('Template_Identity', '0')):.2f}",
                        'coverage': f"{float(result.get('Template_Coverage', '0')):.2f}",
                        'depth': f"{float(result.get('Depth', '0')):.2f}"
                    })
                except (ValueError, KeyError):
                    continue
            
            template_data['plasmids'] = formatted_plasmids
            template_data['plasmid_count'] = len(formatted_plasmids)
        
        # Get additional information
        template_data['mlst_result'] = get_mlst_result(args.output) if hasattr(args, 'output') else 'N/A'
        template_data['phenotype_count'] = len(phenotypes) if phenotypes else 0
        template_data['phenotypes'] = sorted(phenotypes) if phenotypes else []
        
        # Generate drug class plot
        template_data['drug_class_plot'] = generate_drug_class_distribution(amr_results, gene_data)
        
        # Setup Jinja2 environment and render
        env = Environment(
            loader=FileSystemLoader(package_dir / 'templates'),
            autoescape=True
        )
        template = env.get_template('report.html')
        html_content = template.render(**template_data)
        
        # Create PDF
        pdf_path = Path(args.output) / f"{args.name}_report.pdf"
        HTML(string=html_content).write_pdf(
            pdf_path,
            stylesheets=[CSS(assets_dir / 'style.css')]
        )
        
        return pdf_path
    
    except Exception:
        # If anything fails during PDF generation, try one last time with minimal data
        try:
            env = Environment(
                loader=FileSystemLoader(package_dir / 'templates'),
                autoescape=True
            )
            template = env.get_template('report.html')
            html_content = template.render(**template_data)
            
            pdf_path = Path(args.output) / f"{args.name}_report.pdf"
            HTML(string=html_content).write_pdf(
                pdf_path,
                stylesheets=[CSS(assets_dir / 'style.css')]
            )
            
            return pdf_path
        except Exception:
            return None

def get_mlst_result(output_dir):
    """Get MLST result from the MLST output file"""
    try:
        with open(output_dir + "/mlst/mlst_results.tsv", 'r') as f:
            # Skip header line
            next(f)
            line = f.readline().strip()
            if line:
                st_number = line.split('\t')[1]
                return f"ST{st_number}"
            return "N/A"
    except:
        return "N/A"

def get_virulence_count(output_dir, is_ecoli):
    """Get virulence genes count if E. coli"""
    if not is_ecoli:
        return 0
    try:
        virulence_results = read_tab_separated_file(output_dir + "/virulence.res")
        return len(virulence_results)
    except:
        return 0
    
    

def get_file_stats(input_file, output_dir):
    """Get file stats using quick commands and mapstat info.
    
    If some values cannot be found (e.g., mapstat not found or missing entries),
    only those values are 'N/A'. The absence of one file does not revert all 
    values to defaults.
    """
    file_name = 'Unknown'
    file_size = 'N/A'
    fragment_count = None
    kma_version = None
    
    # Attempt to get file_name and file_size from input file
    if os.path.exists(input_file):
        file_name = os.path.basename(input_file)
        try:
            file_size_cmd = f"ls -lh {input_file} | cut -d' ' -f5"
            file_size = subprocess.check_output(file_size_cmd, shell=True).decode().strip()
            if not file_size:
                file_size = 'N/A'
        except:
            file_size = 'N/A'
    
    # Attempt to get fragment_count and kma_version from mapstat
    mapstat_file = os.path.join(output_dir, "bacteria_alignment.mapstat")
    if os.path.exists(mapstat_file):
        try:
            with open(mapstat_file, 'r') as f:
                for line in f:
                    if line.startswith("## fragmentCount"):
                        parts = line.strip().split('\t')
                        if len(parts) > 1:
                            try:
                                fragment_count = int(parts[1])
                            except ValueError:
                                fragment_count = None
                    elif line.startswith("## version"):
                        parts = line.strip().split('\t')
                        if len(parts) > 1:
                            kma_version = parts[1]
        except:
            # If parsing fails for some reason, just continue with what we have
            pass

    return {
        'file_name': file_name if file_name else 'Unknown',
        'file_size': file_size if file_size else 'N/A',
        'read_count': f"{fragment_count:,}" if fragment_count is not None else 'N/A',
        'kma_version': kma_version if kma_version else 'N/A'
    }

    
    
def generate_drug_class_distribution(amr_results, gene_data):
    """Generate stacked bar chart of AMR genes by drug class"""
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        from collections import defaultdict
        import numpy as np
        
        if not amr_results or not gene_data:
            # Create an empty plot with a message
            fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
            ax.text(0.5, 0.5, 'No AMR data available', 
                   ha='center', va='center', fontsize=12)
            ax.set_xticks([])
            ax.set_yticks([])
            
            buf = BytesIO()
            plt.savefig(buf, format='png', dpi=300, bbox_inches='tight',
                       pad_inches=0.1, facecolor='white', edgecolor='none')
            plt.close()
            buf.seek(0)
            return base64.b64encode(buf.read()).decode('utf-8')
        
        matplotlib.rcParams['figure.dpi'] = 300
        matplotlib.rcParams['figure.figsize'] = [8, 4]
        
        class_counts = defaultdict(lambda: {'perfect': 0, 'good': 0, 'partial': 0})
        
        for result in amr_results:
            try:
                gene_id = result.get('#Template', '')
                identity = float(result.get('Template_Identity', '0').strip())
                coverage = float(result.get('Template_Coverage', '0').strip())
                
                quality = 'partial'
                if identity >= 100.0 and coverage >= 100.0:
                    quality = 'perfect'
                elif identity >= 100.0 or coverage >= 100.0:
                    quality = 'good'
                
                for gene in gene_data:
                    if gene.get('Gene_accession no.') == gene_id:
                        classes = [cls.strip() for cls in gene.get('Class', '').split(',') if cls.strip()]
                        for drug_class in classes:
                            class_counts[drug_class][quality] += 1
            except (ValueError, AttributeError):
                continue
        
        # Only continue if we have data to plot
        if not class_counts:
            fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
            ax.text(0.5, 0.5, 'No valid AMR data to display', 
                   ha='center', va='center', fontsize=12)
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            classes = list(class_counts.keys())
            perfect_counts = [class_counts[c]['perfect'] for c in classes]
            good_counts = [class_counts[c]['good'] for c in classes]
            partial_counts = [class_counts[c]['partial'] for c in classes]
            
            fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
            colors = ['#6fbf50', '#bfddbe', '#9d9c9c']
            
            ax.bar(classes, perfect_counts, label='Perfect match', color=colors[0])
            ax.bar(classes, good_counts, bottom=perfect_counts, label='Good match', color=colors[1])
            ax.bar(classes, partial_counts, bottom=np.array(perfect_counts) + np.array(good_counts),
                    label='Partial match', color=colors[2])
            
            plt.xticks(rotation=45, ha='right', fontsize=9)
            plt.ylabel('Number of Genes', fontsize=10)
            ax.set_title('AMR Genes by Drug Class', fontsize=12, pad=20)
            ax.legend(fontsize=8, loc='upper right')
            ax.yaxis.grid(True, linestyle='--', alpha=0.3)
            
            plt.tight_layout()
        
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=300, bbox_inches='tight',
                   pad_inches=0.1, facecolor='white', edgecolor='none')
        plt.close()
        
        buf.seek(0)
        return base64.b64encode(buf.read()).decode('utf-8')
    except Exception:
        # Return a base64 encoded empty plot if anything fails
        fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
        ax.text(0.5, 0.5, 'Error generating AMR plot', 
               ha='center', va='center', fontsize=12)
        ax.set_xticks([])
        ax.set_yticks([])
        
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=300, bbox_inches='tight',
                   pad_inches=0.1, facecolor='white', edgecolor='none')
        plt.close()
        
        buf.seek(0)
        return base64.b64encode(buf.read()).decode('utf-8')
    
    
############################## PDF REPORT DONE ################################
