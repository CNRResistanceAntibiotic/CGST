#!/usr/bin/env python3
"""
Copyright 2020 Aurélien BIRER (abirer36@gmail.com)
https://github.com/Nilad/CGST.git

This script is the detection file part of CGST.

This file is part of CGST. CGST is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. CGST is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with CGST. If
not, see <http://www.gnu.org/licenses/>.
"""
from csv import DictReader, DictWriter, writer
import os
import shutil
from subprocess import Popen, PIPE, STDOUT
from collections import OrderedDict
from Bio import SeqIO

from cgst.log import section_header, explanation, log_process_with_output_file, log
from cgst.utils import read_output_mentalist


def mentalist_detection(r1, r2, database, work_dir_path, name, species_full, threads):
    """
    This function manage the MentaLiST detection
    :param r1: The R1 fastq file path
    :param r2: The R2 fastq file path
    :param database: The CGST database path
    :param work_dir_path: The working directory path
    :param name: The name of
    :param species_full: The scientific name of the strain
    :param threads: The number of threads to allocate
    """
    kmer_threshold = 5
    kmer_build = 31

    species_full = species_full.strip().rstrip()
    species = species_full.split(" ")[0][0].lower() + species_full.split(" ")[1].lower()

    exe = shutil.which("mentalist")
    exe_parse = shutil.which("parse_novel_alleles.py")
    exe_update = shutil.which("update_fasta_db.py")

    cgmlst_database_path = os.path.join(database, "cgMLST", f"{species}")
    fasta_db_path = db_path = cgmlst_dir_path = ""
    cg_db_dict = {}

    for file_1 in os.listdir(cgmlst_database_path):
        file_1_path = os.path.join(cgmlst_database_path, file_1)
        if os.path.isdir(file_1_path):
            cgmlst_dir_path = file_1_path
            for file_2 in os.listdir(cgmlst_dir_path):
                file_2_path = os.path.join(cgmlst_dir_path, file_2)
                if "_fasta" in file_2 and os.path.isdir(file_2_path):
                    fasta_db_path = file_2_path
                    db_name = file_2.split("_fasta")[0]+".db"
                    db_path = os.path.join(cgmlst_dir_path, db_name)
        cg_db_dict[file_1] = {'db_path': db_path, 'fasta_db_path': fasta_db_path}

    # make detection on each core-genome available for the species
    for name_db, value_hash in cg_db_dict.items():
        i = 1
        db_path = value_hash["db_path"]
        fasta_db_path = value_hash["fasta_db_path"]
        section_header('Launch MentaLiST:')
        explanation('Name Database  : {0}'.format(name_db))
        explanation('Database Path : {0}'.format(db_path))
        explanation('Fasta Path : {0}'.format(fasta_db_path))
        output_final = os.path.join(work_dir_path, f"{name}_{name_db}_output_final")

        # If the kmer index database do not exist -> create it
        if not os.path.exists(db_path):
            ###################################
            # Run MentaLiST Build DB
            section_header('Create MentaLiST Database')
            explanation('Before run a detection MentaLiST need to construct his own kmer-index database')
            # prepare
            cmd = f"{exe} build_db --db {db_path} -k {kmer_build} -d {fasta_db_path} --threads {threads}"
            log_message = f"Command used : \n {cmd}"
            # launch
            log_file_path = os.path.join(work_dir_path, "logBuildDB_{0}.txt".format(name_db))
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, executable='/bin/bash')
            log_process_with_output_file(process, log_message, log_file_path)

        ###################################
        # Run MentaLiST call
        section_header(f'Run MentaLiST Call : Round {i}')
        explanation('MentaLiST Detection ')
        # prepare
        output = os.path.join(work_dir_path, f"output_mentalist_{name_db}_{i}")
        cmd = f"{exe} call --db {db_path} --output_votes -o {output} -1 {r1} -2 {r2} --kt {kmer_threshold}"
        log_message = f"Command used : \n {cmd}\n"
        # launch
        process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, executable='/bin/bash')
        log_file_path = os.path.join(work_dir_path, f"logMentaLiST_{name_db}_{i}.txt")
        log_process_with_output_file(process, log_message, log_file_path)
        fasta_novel_st = os.path.join(work_dir_path, f"output_mentalist_{i}.novel.fa")

        while os.path.exists(fasta_novel_st) and os.stat(fasta_novel_st).st_size != 0:
            ###################################
            # check MentaLiST output for Novel and Multiple Votes
            section_header(f'Check MentaLiST Output : Round {i}')
            explanation('Check MentaLiST output for Novel Votes')
            check_mentalist_output(fasta_novel_st, output, fasta_db_path)

            ###################################
            # Run MentaLiST Parse Novel Fasta
            section_header(f'Run MentaLiST parse novel fasta : Round {i}')
            explanation('MentaLiST parse new novel variant ')
            # prepare
            fasta_novel_st = os.path.join(work_dir_path, f"output_mentalist_{name_db}_{i}.novel.fa")
            result_parse_path = os.path.join(work_dir_path, f"all_novel_alleles_{name_db}_{i}")
            cmd = f"{exe_parse} -f {fasta_novel_st} -o {result_parse_path}"
            log_message = f"Command used : \n {cmd}\n"
            # launch
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
            log_file_path = os.path.join(work_dir_path, f"logParseNovelAlleles_{name_db}_{i}.txt")
            log_process_with_output_file(process, log_message, log_file_path)

            ###################################
            # Run MentaLiST Update Fasta DB
            section_header(f'Run MentaLiST update DB fasta : Round {i}')
            explanation('MentaLiST update DB fasta with new novel variant ')
            # prepare
            result_parse_path = os.path.join(work_dir_path, f"all_novel_alleles_{name_db}_{i}")
            cmd = f"{exe_update} -db {fasta_db_path} -n {result_parse_path}.fa"
            log_message = f"Command used : \n {cmd}\n"
            # launch
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
            log_file_path = os.path.join(work_dir_path, f"logCreateNewSchemeWithNovel_{name_db}_{i}.txt")
            log_process_with_output_file(process, log_message, log_file_path)

            ###################################
            # Run MentaLiST Build DB
            section_header(f'Create MentaLiST Database : Round {i}')
            explanation('Before run a detection MentaLiST need to construct his own kmer-index database')
            # prepare
            cmd = f"{exe} build_db --db {db_path} -k {kmer_build} -d {fasta_db_path} --threads {threads}"
            log_message = f"Command used : \n {cmd}"
            # remove previous db
            os.remove(db_path)
            # launch
            log_file_path = os.path.join(work_dir_path, f"logBuildDB_{name_db}_{i}.txt")
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, executable='/bin/bash')
            log_process_with_output_file(process, log_message, log_file_path)

            ###################################
            # Update counter
            i += 1

            ###################################
            # Run MentaLiST call
            section_header(f'Run MentaLiST Call : Round {i}')
            explanation('MentaLiST Detection ')
            # prepare
            output = os.path.join(work_dir_path, f"output_mentalist_{name_db}_{i}")
            cmd = f"{exe} call --db {db_path} --output_votes -o {output} -1 {r1} -2 {r2} --kt {kmer_threshold}"
            log_message = f"Command used : \n {cmd}\n"
            # launch
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, executable='/bin/bash')
            log_file_path = os.path.join(work_dir_path, f"logMentaLiST_{name_db}_{i}.txt")
            log_process_with_output_file(process, log_message, log_file_path)
            fasta_novel_st = os.path.join(work_dir_path, f"output_mentalist_{name_db}_{i}.novel.fa")

            ###################################
            # check MentaLiST output for Novel and Multiple Votes
            section_header(f'Check MentaLiST Output : Round {i}')
            explanation('Check MentaLiST output for Novel Votes and edit false novel in previous output')
            check_mentalist_output(fasta_novel_st, output, fasta_db_path)

        ###################################
        # Process final result
        section_header('Process Final Output')

        # final files
        shutil.move(output, output_final)
        intermediate_file_dir = os.path.join(work_dir_path, "intermediate_files")

        if not os.path.exists(intermediate_file_dir):
            os.mkdir(intermediate_file_dir)

        for file in os.listdir(work_dir_path):
            file_path = os.path.join(work_dir_path, file)
            if not os.path.isdir(file_path):
                if not ("_output_final" in file or "statistics_" in file or "combination_result_"):
                    shutil.move(file_path, os.path.join(intermediate_file_dir, file))

        # add ST in modify output_mentalist
        st_dict = add_st_to_output_mentalist(output_final, name)

        # Statistic
        detection_result_dict = read_output_mentalist(output_final, name)

        count_none = count_loc = count_low_cov = count_perfect = count_multi = 0

        for locus, value_dict in detection_result_dict.items():
            if locus == "Sample" or locus == "ST" or locus == "clonal_complex" or locus in st_dict:
                continue
            elif value_dict[name] == "0?":
                count_none += 1
                count_loc += 1
            elif "-" in value_dict[name]:
                count_low_cov += 1
                count_loc += 1
            elif "+" in value_dict[name]:
                count_multi += 1
                count_loc += 1
            else:
                count_perfect += 1
                count_loc += 1

        explanation(f'Total Locus'
                    f' : {count_loc}')
        explanation(f'Perfect Locus'
                    f' : {count_perfect} -> {round((count_perfect / count_loc * 100), 2)}% of total locus')
        explanation(f'Multiple Locus'
                    f' : {count_multi} -> {round((count_multi / count_loc * 100), 2)}% of total locus')
        explanation(f'Low Coverage Locus'
                    f' : {count_low_cov} -> {round((count_low_cov / count_loc * 100), 2)}% of total locus')
        explanation(f'None Locus'
                    f' : {count_none} -> {round((count_none / count_loc * 100), 2)}% of total locus')
        stats_file = os.path.join(work_dir_path, "statistics_{0}.tsv".format(name_db))

        with open(stats_file, "w") as output_file:
            csv_writer = writer(output_file, delimiter="\t")
            csv_writer.writerow(["Total Locus", count_loc])
            csv_writer.writerow(["Resume", "Count", "Percentage on Total locus"])
            csv_writer.writerow(["Perfect Locus", count_perfect, (count_perfect / count_loc * 100)])
            csv_writer.writerow(["Multiple Locus", count_multi, (count_multi / count_loc * 100)])
            csv_writer.writerow(["Low Coverage Locus", count_low_cov, (count_low_cov / count_loc * 100)])
            csv_writer.writerow(["None Locus", count_none, (count_none / count_loc * 100)])

        # Load Known Combination
        known_comb_path = os.path.join(os.path.dirname(db_path), "combination_list.tsv")

        known_comb_dict = {}

        if os.path.exists(known_comb_path):
            with open(known_comb_path, "r") as file:
                reader = DictReader(file, delimiter="\t")
                for row in reader:
                    known_comb_dict[row["Name"]] = row["Combination"]
        else:
            explanation(f"combination file not found for {species_full} at {known_comb_path}")

        # Search For Known Combination
        final_file_path = os.path.join(work_dir_path, f"{name}_{name_db}_output_final")

        detection_result_dict = read_output_mentalist(final_file_path, name)

        comb_strain_list = []

        for locus, sample_dict in detection_result_dict.items():
            if locus == "Sample" or locus == "ST" or locus == "clonal_complex":
                continue
            else:
                comb_strain_list.append(f"{locus}:{sample_dict[name]}")
        if known_comb_dict:
            # search for each known combination
            combine_result_file = os.path.join(work_dir_path, f"combination_result_{name_db}.tsv")
            with open(combine_result_file, "w") as combine_file:
                csv_writer = writer(combine_file, delimiter="\t")
                csv_writer.writerow(["Combination Name", "Count Reference Locus", "Count Sample Locus", "Ratio",
                                     "Comment"])
                for name_comb, combs in known_comb_dict.items():
                    combination_known = combs.split(",")
                    result = list(set(comb_strain_list).intersection(combination_known))
                    ratio = round((len(result) / len(combination_known)) * 100, 2)
                    if ratio == 100:
                        comment = "Perfect"
                    elif ratio >= 98:
                        comment = "Very Close"
                    elif ratio >= 90:
                        comment = "Close"
                    elif ratio >= 80:
                        comment = "Like"
                    elif ratio >= 70:
                        comment = "Close Like"
                    else:
                        comment = "No relevant"
                    csv_writer.writerow([name_comb, len(combination_known), len(result), ratio, comment])
        section_header('Finish {0} Coregenome Analysis'.format(name_db))
    section_header('Finish ALL CoreGenome Analysis')


def ariba_detection(r1, r2, database, work_dir_path, name, species_full):
    """
    This function manage the ARIBA detection
    :param r1: The R1 fastq file path
    :param r2: The R2 fastq file path
    :param database: The CGST database path
    :param work_dir_path: The working directory path
    :param name: The name used
    :param species_full: The scientific name of the strain
    """
    exe = shutil.which("ariba")
    species_full = species_full.strip().rstrip()
    species = species_full.split(" ")[0][0].lower() + species_full.split(" ")[1].lower()
    db_path_species = os.path.join(database, "ariba", f"{species}")
    reference_dict = get_reference_ariba_mlst()

    for file in os.listdir(db_path_species):
        file_path = os.path.join(db_path_species, file)
        if os.path.isdir(file_path):
            db_path = os.path.join(file_path, "ref_db")
            ###################################
            # Run Ariba call
            section_header('Run Ariba Call')
            explanation('Ariba Detection ')
            # prepare
            output = os.path.join(work_dir_path, f"{name}_output_ariba_{reference_dict[file]}")
            cmd = f"{exe} run {db_path} {r1} {r2} {output}"
            log_message = f"Command used : \n {cmd}\n"
            # launch
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
            log_file_path = os.path.join(work_dir_path, f"logAriba_{reference_dict[file]}.txt")
            log_process_with_output_file(process, log_message, log_file_path)


def check_mentalist_output(fasta_novel_st, output, fasta_db_path):
    """
    This function manage the check of novels alleles detected by MentaLiST
    :param fasta_novel_st: The fasta file path with novel alleles
    :param output: The output file directory path
    :param fasta_db_path: The fasta database directory path
    """
    fasta_corrected = os.path.join(os.path.dirname(fasta_novel_st),
                                   f"{os.path.basename(fasta_novel_st)}_corrected")
    output_corrected = os.path.join(os.path.dirname(output),
                                    f"{os.path.basename(output)}_corrected")

    false_novel_allele_dict = {}
    with open(fasta_corrected, "w") as output_handle:
        for record in SeqIO.parse(fasta_novel_st, "fasta"):
            locus_list = record.id.split("_")
            # remove comment id to get the identifier file
            locus_list.pop()
            locus = "_".join(locus_list)
            for file in os.listdir(fasta_db_path):
                if locus in file:
                    locus_file_path = os.path.join(fasta_db_path, file)
                    pivot_remove = False
                    for record_db in SeqIO.parse(locus_file_path, "fasta"):
                        if str(record.seq) == str(record_db.seq):
                            log(f"False New alleles for locus {locus}."
                                f" Add for this locus this allele {str(record_db.id)}")
                            false_novel_allele_dict[locus] = str(record_db.id)
                            pivot_remove = True
                            break
                        else:
                            continue
                    # write or not duplicate novel
                    if pivot_remove:
                        continue
                    else:
                        SeqIO.write(record, output_handle, "fasta")
                else:
                    continue

    os.remove(fasta_novel_st)
    shutil.move(fasta_corrected, fasta_novel_st)
    output_corrected_dict = OrderedDict()
    with open(output, "r") as output_file:
        reader = DictReader(output_file, delimiter="\t")
        for row in reader:
            for locus, value in row.items():
                if "N" in value:
                    if locus in false_novel_allele_dict:
                        # change to the true allele
                        output_corrected_dict[locus] = false_novel_allele_dict[locus].split("_")[-1]
                    else:
                        output_corrected_dict[locus] = value
                        continue
                else:
                    output_corrected_dict[locus] = value
                    continue

    with open(output_corrected, "w") as output_corrected_file:
        csv_writer = writer(output_corrected_file, delimiter="\t")
        csv_writer.writerow(dict(output_corrected_dict))
        csv_writer.writerow(dict(output_corrected_dict).values())

    os.remove(output)
    shutil.move(output_corrected, output)


def get_reference_ariba_mlst():
    """
    This function manage the ARIBA MLST database
    :return: A dictionary of the ARIBA MLST database
    """
    reference_dict = {}
    path = os.path.abspath(os.path.dirname(__file__))
    reference_file = os.path.join(path, "reference_mlst.tsv")
    with open(reference_file, "r") as file:
        reader = DictReader(file, delimiter="\t")
        for row in reader:
            reference_dict[row["Schema"]] = row["Reference"]
    return reference_dict


def add_st_to_output_mentalist(file_path, name):
    """
    This function manage the addition of the ST value to MentaLiST output
    :param file_path: The output MentaLiST file path
    :param name: the name used
    :return A dict of ST values
    """
    file_tmp_path = file_path + "_tmp"
    with open(file_tmp_path, 'w', newline='')as file_tmp:
        with open(file_path, "r") as detection_file:
            st_dict = {}
            # add ST ariba
            for file in os.listdir(os.path.dirname(file_path)):
                if "output_ariba" in file:
                    schema = file.split(f"{name}_output_ariba_")[1]
                    st_value = ""
                    mlst_report_file = os.path.join(os.path.dirname(file_path), file, "mlst_report.tsv")
                    with open(mlst_report_file, "r") as mlst_file:
                        reader_mlst = DictReader(mlst_file, delimiter="\t")
                        for row_mlst in reader_mlst:
                            st_value = row_mlst["ST"]
                    st_dict[f"ST {schema}"] = st_value
            reader = DictReader(detection_file, delimiter="\t")
            row_headers = reader.fieldnames + [*st_dict]
            csv_writer = DictWriter(file_tmp, row_headers, delimiter="\t")
            csv_writer.writeheader()
            for row in reader:
                csv_writer.writerow({**row, **st_dict})
    os.remove(file_path)
    shutil.move(file_tmp_path, file_path)
    return st_dict


def main_detection(r1, r2, database, work_dir_path, name, species_full, force, only_mentalist, threads):
    """
    This function manage the detection functions
    :param r1: The R1 fastq file path
    :param r2: The R2 fastq file path
    :param database: The CGST database path
    :param work_dir_path: The working directory path
    :param name: The name used
    :param species_full: The scientific name of the strain
    :param force: Force the output
    :param only_mentalist: Boolean to run only MentaLiST
    :param threads: The number of threads to allocate
    """
    if not os.path.exists(work_dir_path):
        os.mkdir(work_dir_path)
    elif force and os.path.exists(work_dir_path):
        shutil.rmtree(work_dir_path)
        os.mkdir(work_dir_path)

    if not only_mentalist:
        ariba_detection(r1, r2, database, work_dir_path, name, species_full)
    mentalist_detection(r1, r2, database, work_dir_path, name, species_full, threads)
