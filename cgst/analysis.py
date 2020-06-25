#!/usr/bin/env python3
"""
Copyright 2020 Aurélien BIRER (abirer36@gmail.com)
https://github.com/Nilad/CGST.git

This script is the analysis file part of CGST.

This file is part of CGST. CGST is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. CGST is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with CGST. If
not, see <http://www.gnu.org/licenses/>.
"""
from csv import DictReader, writer
import multiprocessing
import os
import shutil
from subprocess import Popen, PIPE, STDOUT
from collections import Counter
from datetime import datetime

from itertools import combinations
from shutil import rmtree

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from cgst.log import section_header, log_process_with_output_file, explanation, tool_error_log
from cgst.utils import read_output_mentalist


def main_analysis(detection_dir, database, work_dir, species_full, force, threads):
    """
    This function that manage the analysis of CGST
    :param detection_dir: The directory MentaLiST output
    :param database: The CGST database path
    :param work_dir: The analysis working directory path
    :param species_full: The scientific name of the strain
    :param force: Force the output
    :param threads: The number of threads to allocate
    """
    level_gg_dict = {"Diff_5_alleles": {"lvl": 5, "list_sample": []},
                     "Diff_10_alleles": {"lvl": 10, "list_sample": []},
                     "Diff_25_alleles": {"lvl": 25, "list_sample": []},
                     "Diff_50_alleles": {"lvl": 50, "list_sample": []},
                     "Diff_100_alleles": {"lvl": 100, "list_sample": []},
                     "Diff_150_alleles": {"lvl": 150, "list_sample": []},
                     "Diff_200_alleles": {"lvl": 200, "list_sample": []},
                     "Diff_300_alleles": {"lvl": 300, "list_sample": []}}
    species_full = species_full.strip().rstrip()
    species = species_full.split(" ")[0][0].lower() + species_full.split(" ")[1].lower()

    ###################################

    if not os.path.exists(work_dir):
        os.mkdir(work_dir)

    ###################################
    # Load Detection result
    section_header(f'Load Detection Result {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    detection_result_dict = {}
    detection_result_dict_list = []
    for file in os.listdir(detection_dir):
        if "_output_final" in file:
            strain_name = file.split("_output_final")[0]
            file_path = os.path.join(detection_dir, file)
            detection_result_dict_tmp = read_output_mentalist(file_path, strain_name)
            detection_result_dict_list.append(detection_result_dict_tmp)

    if len(detection_result_dict_list) >= 5:
        pass
    else:
        tool_error_log(f"Need more strains that just {len(detection_result_dict_list)}")
        exit()

    for detect_dict in detection_result_dict_list:
        for key, value_dict in detect_dict.items():
            if key in detection_result_dict:
                value_detect_dict = detection_result_dict[key]
                for sample, value in value_dict.items():
                    value_detect_dict[sample] = value
            else:
                detection_result_dict[key] = value_dict
    ##########
    # Get Difference alleles
    section_header(f'Analysis {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    detection_result_only_diff_dict = {}
    distance_dict = {}
    sample_list = []
    for key, sample_dict in detection_result_dict.items():
        if key == "Sample" or "clonal_complex" in key or "ST" in key:
            if key == "Sample":
                sample_list = [*sample_dict]
            continue
        else:
            element_list = []
            pivot = False
            for sample_name in sample_list:
                if "-" in sample_dict[sample_name] or "+" in sample_dict[sample_name] \
                        or "0" == sample_dict[sample_name] or "?" in sample_dict[sample_name]:
                    pivot = True
                    break
                else:
                    element_list.append(sample_dict[sample_name])
            if not pivot:
                # test if the element are the same or not
                if len(set(element_list)) == 1:
                    continue
                else:
                    detection_result_only_diff_dict[key] = sample_dict
                    for sample_name_1 in sample_list:
                        for sample_name_2 in sample_list:

                            if sample_dict[sample_name_1] == sample_dict[sample_name_2]:
                                continue
                            else:
                                if sample_name_1 in distance_dict:
                                    value1_dict = distance_dict[sample_name_1]
                                    if sample_name_2 in value1_dict:
                                        value1_count = value1_dict[sample_name_2]
                                        value1_dict[sample_name_2] = value1_count + 1
                                    else:
                                        value1_dict[sample_name_2] = 1

                                else:
                                    distance_dict[sample_name_1] = {sample_name_2: 1}
    explanation(f"Number of different relevant locus : {len(detection_result_only_diff_dict)}")
    ###################################
    # Get Similarity Matrix
    section_header(f'Get Similarity Matrix {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    sample_set_list = list(combinations(sample_list, 2))
    similarity_dict = {}
    for locus, sample_dict in detection_result_only_diff_dict.items():
        for sample_set in sample_set_list:
            if sample_dict[sample_set[0]] == sample_dict[sample_set[1]]:
                continue
            else:
                if sample_set in similarity_dict:
                    value_similarity = similarity_dict[sample_set]
                    value_similarity += 1
                    similarity_dict[sample_set] = value_similarity
                else:
                    similarity_dict[sample_set] = 1

    ###################################

    combination_dir = os.path.join(work_dir, "combination")
    if not os.path.exists(combination_dir):
        os.mkdir(combination_dir)

    cluster_dir = os.path.join(work_dir, "cluster")
    if not os.path.exists(cluster_dir):
        os.mkdir(cluster_dir)

    phylotree_dir = os.path.join(work_dir, "phylotree")
    if not os.path.exists(phylotree_dir):
        os.mkdir(phylotree_dir)

    ###################################
    # Write Similarity Matrix
    section_header(f'Write Similarity Matrix {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    similarity_file = os.path.join(combination_dir, "similarity_matrix.tsv")
    with open(similarity_file, "w") as similarity:
        writer_csv = writer(similarity, delimiter='\t')
        writer_csv.writerow([""] + sample_list)
        for sample_1 in sample_list:
            val_list = []
            for sample_2 in sample_list:
                tuple_1 = (sample_1, sample_2)
                tuple_2 = (sample_2, sample_1)
                if tuple_1 in similarity_dict:
                    val_list.append(similarity_dict[tuple_1])
                    continue
                if tuple_2 in similarity_dict:
                    val_list.append(similarity_dict[tuple_2])
                    continue
                if sample_1 == sample_2:
                    val_list.append("0")
                    continue
                else:
                    val_list.append("0")
            writer_csv.writerow([sample_1] + val_list)
    ###################################
    # R Script
    section_header(f'Execute RScript {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    ex_r = shutil.which("Rscript")
    r_script_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "r_script.R")
    cmd = f"{ex_r} {r_script_file} --wd {combination_dir}"
    log_message = f"Command used : \n {cmd}\n"
    # launch
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    log_file_path = os.path.join(combination_dir, "logR.txt")
    log_process_with_output_file(process, log_message, log_file_path)
    ###########################
    # Create Group File
    section_header(f'Create Group {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    r_result_file = os.path.join(combination_dir, "groups.tsv")
    groups_dict = {}
    strains_groups_list = []
    with open(r_result_file, "r") as groups_file:
        reader = DictReader(groups_file, delimiter="\t")
        headers = reader.fieldnames
        for row in reader:
            strains_groups_list.append(row[""])
            for head in headers:
                if head:
                    if head in groups_dict:
                        val_list = groups_dict[head]
                        pivot = True
                        for i, group_d in enumerate(val_list):
                            if row[head] in group_d["group"]:
                                pivot = False
                                if row[""] not in group_d["strains"]:
                                    str_list = group_d["strains"]
                                    str_list.append(row[""])
                                    group_d["strains"] = str_list
                                    continue
                                else:
                                    continue
                        if pivot:
                            val_list.append({"group": row[head], "strains": [row[""]]})
                            continue
                        groups_dict[head] = val_list
                        continue
                    else:
                        groups_dict[head] = [{"group": row[head], "strains": [row[""]]}]
                        continue

    strains_groups_list = list(set(strains_groups_list))
    groups_dict["1"] = [{"group": "1", "strains": strains_groups_list}]

    ###################################
    # Add allele to group
    for key_class, list_value in groups_dict.items():
        for group in list_value:
            share_90_allele_list = []
            share_strict_allele_list = []
            for key, sample_dict in detection_result_dict.items():
                if key == "Sample" or "clonal_complex" in key or "ST" in key:
                    if "ST" in key:
                        st_tmp_list = []
                        for strain in group["strains"]:
                            st_tmp_list.append(sample_dict[strain])
                        group["name_group"] = f"{species}-{key}{':'.join(list(set(st_tmp_list)))}"
                        continue
                else:
                    allele_tmp_list = []
                    for strain in group["strains"]:
                        allele_tmp_list.append(sample_dict[strain])
                    #######
                    unique = list(Counter(allele_tmp_list).keys())  # equals to list(set(words))

                    value = list(Counter(allele_tmp_list).values())
                    if len(unique) == 1:
                        share_strict_allele_list.append({"locus": key, "allele": unique[0]})
                        share_90_allele_list.append({"locus": key, "allele": unique[0]})
                    else:
                        for i, val in enumerate(value):
                            if (val / len(group["strains"])) * 100 >= 90:
                                share_90_allele_list.append({"locus": key, "allele": unique[i]})
                                break
                            else:
                                continue
                        continue
            group["share_strict_allele"] = share_strict_allele_list
            group["share_90_allele"] = share_90_allele_list
            group["count_strict_allele"] = len(share_strict_allele_list)
            group["count_90_allele"] = len(share_90_allele_list)
    ###################################
    # Create level lists
    section_header(f'Create Lvl list {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    for name_lvl, lvl_dict in level_gg_dict.items():
        for sample_name, count_dict in distance_dict.items():
            group_lvl_list = []
            for sample_name_2 in sample_list:
                if sample_name == sample_name_2:
                    continue
                if sample_name_2 not in count_dict:
                    group_lvl_list.append(sample_name_2)
                    continue
                if count_dict[sample_name_2] <= lvl_dict["lvl"]:
                    group_lvl_list.append(sample_name_2)
                    continue
            if group_lvl_list:
                group_lvl_list.append(sample_name)
                group_lvl_list.sort()
                if "list_sample" in lvl_dict:
                    list_tmp = lvl_dict["list_sample"]
                    # check if the list is already present
                    if group_lvl_list in list_tmp:
                        continue
                    list_tmp.append(group_lvl_list)
                    lvl_dict["list_sample"] = list_tmp
                else:
                    lvl_dict["list_sample"] = [group_lvl_list]
    ###################################
    # Create Multi fasta
    section_header(f'Create Multiple Fasta Files {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

    cgmlst_database_path = os.path.join(database, "cgMLST", f"{species}")

    fasta_db_path = cgmlst_dir_path = ""

    for file_1 in os.listdir(cgmlst_database_path):
        if "cgmlst-org" == file_1 or "cnr" == file_1:
            cgmlst_dir_path = os.path.join(cgmlst_database_path, file_1)
            for file_2 in os.listdir(cgmlst_dir_path):
                file_2_path = os.path.join(cgmlst_dir_path, file_2)
                if ".db" not in file_2 and "_fasta" in file_2 and os.path.isdir(file_2_path):
                    fasta_db_path = file_2_path

    output_dir_msa = os.path.join(phylotree_dir, "msa")
    if not os.path.exists(output_dir_msa):
        os.mkdir(output_dir_msa)
    elif force and os.path.exists(output_dir_msa):
        rmtree(output_dir_msa)
        os.mkdir(output_dir_msa)

    pivot_first_loop = True
    sequence_dict = {}
    final_resume_aln_dict = {}
    count_locus_cg = 0
    for file in os.listdir(fasta_db_path):
        if ".fasta" in file:
            count_locus_cg += 1

    output_aln_file_list = []

    pool = multiprocessing.Pool(processes=int(threads))
    list_jobs = []

    for locus_name, sample_dict in detection_result_only_diff_dict.items():
        fasta_file = ""
        for file in os.listdir(fasta_db_path):
            if file == f"{locus_name}.fasta":
                fasta_file = os.path.join(fasta_db_path, f"{locus_name}.fasta")
                break
        output_fasta_file = os.path.join(output_dir_msa, os.path.basename(fasta_file))
        with open(output_fasta_file, "w") as output_fasta:
            record_dict = SeqIO.index(fasta_file, "fasta")
            for sample_name, number_allele in sample_dict.items():
                seq = record_dict[f"{locus_name}_{number_allele}"]
                seq.id = sample_name
                SeqIO.write(seq, output_fasta, "fasta")
            record_dict.close()

        ###################################
        # MAFFT - MSA
        output_aln_file = os.path.join(output_dir_msa, f"{os.path.basename(fasta_file).split('.')[0]}.aln")
        output_aln_file_list.append(output_aln_file)

        list_jobs.append([output_dir_msa, output_fasta_file, locus_name, output_aln_file])

    pool.starmap(mafft, list_jobs)

    stop = 0
    for output_aln_file in output_aln_file_list:
        ###################################
        # EXPLOIT OUTPUT MAFFT
        with open(output_aln_file, "r") as handle:
            record_aln_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

        for id_seq in sorted(record_aln_dict, key=lambda id_s: len(record_aln_dict[id_s].seq), reverse=True):
            if pivot_first_loop:
                start = 1
                stop = len(record_aln_dict[id_seq].seq)
            else:
                start = stop + 1
                stop = start + len(record_aln_dict[id_seq].seq)
            final_resume_aln_dict[record_aln_dict[id_seq].description.split(" ")[-1]] = {"start": start, "stop": stop, "length": stop - start}
            break
        for sample in sample_list:
            if sample in sequence_dict:
                sequence = sequence_dict[sample]
                if sample not in record_aln_dict:
                    sequence = sequence + "-" * final_resume_aln_dict[record_aln_dict[sample].description.split(" ")[-1]]["length"]

                else:
                    sequence = sequence + record_aln_dict[sample].seq
                sequence_dict[sample] = sequence
                continue
            else:
                if sample not in record_aln_dict:
                    sequence = "-" + "-" * final_resume_aln_dict[record_aln_dict[sample].description.split(" ")[-1]]["length"]
                else:
                    sequence = record_aln_dict[sample].seq
                sequence_dict[sample] = sequence
                continue
        pivot_first_loop = False
    ###################################
    # Write Groups
    all_groups_alleles_file = os.path.join(combination_dir, "groups_alleles.tsv")
    with open(all_groups_alleles_file, "w") as all_groups_alleles:
        writer_group = writer(all_groups_alleles, delimiter='\t')
        writer_group.writerow(
            ["Class", "Group", "Strains", "Count Strains", "Name Group", "Alleles strict", "Count Alleles strict",
             "Strict Coverage CG", "Alleles 90%", "Count Alleles 90%", "90% Coverage CG"])
        for key_class, list_value in groups_dict.items():
            for ele_d in list_value:
                strict_string_list = []
                for ele_strict in ele_d["share_strict_allele"]:
                    strict_string_list.append(f"{ele_strict['locus']}:{ele_strict['allele']}")
                ninety_string_list = []
                for ele_ninety in ele_d["share_90_allele"]:
                    ninety_string_list.append(f"{ele_ninety['locus']}:{ele_ninety['allele']}")
                writer_group.writerow(
                    [key_class, ele_d["group"], ele_d["strains"], len(ele_d["strains"]), ele_d["name_group"],
                     ",".join(strict_string_list), ele_d["count_strict_allele"],
                     (ele_d["count_strict_allele"] / count_locus_cg) * 100, ",".join(ninety_string_list),
                     ele_d["count_90_allele"], (ele_d["count_90_allele"] / count_locus_cg) * 100])

    ###################################
    # Create MSA
    section_header(f'Create Alignment file {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    output_aln_final_file = os.path.join(phylotree_dir, "core-genome.aln")
    with open(output_aln_final_file, "w") as output_aln_final_handle:
        for sample in sample_list:
            record = SeqRecord(sequence_dict[sample], id=sample)
            SeqIO.write(record, output_aln_final_handle, "fasta")

    output_resume_final_file = os.path.join(phylotree_dir, "resume_core-genome.tsv")
    with open(output_resume_final_file, "w") as output_resume_final_handle:
        writer_resume = writer(output_resume_final_handle, delimiter="\t")
        writer_resume.writerow(["Gene", "Start", "Stop", "Length"])
        for gene, value_dict in final_resume_aln_dict.items():
            writer_resume.writerow([gene, value_dict["start"], value_dict["stop"], value_dict["length"]])
    ###################################
    # Gubbins
    section_header(f'Delete Recombination with Gubbins {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

    ex_gubbins = shutil.which("run_gubbins")

    gubbins_work_dir = os.path.join(phylotree_dir, "gubbins")

    if not os.path.exists(gubbins_work_dir):
        os.mkdir(gubbins_work_dir)
    elif force and os.path.exists(gubbins_work_dir):
        rmtree(gubbins_work_dir)
        os.mkdir(gubbins_work_dir)
    cmd = f"{ex_gubbins} -p {gubbins_work_dir}/gubbins --threads {threads} {output_aln_final_file}"
    log_message = f"Command used : \n {cmd}\n"
    # launch
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    log_file_path = os.path.join(gubbins_work_dir, "logGubbins.txt")
    log_process_with_output_file(process, log_message, log_file_path)
    ###################################
    # RAXML-ng
    section_header(f'Phylotree with RaXML-ng {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    ex_raxml_ng = shutil.which("raxml-ng")
    raxml_work_dir = os.path.join(phylotree_dir, "raxml-ng")
    raxml_prefix = os.path.join(raxml_work_dir, "raxml-ng")

    if not os.path.exists(raxml_work_dir):
        os.mkdir(raxml_work_dir)
    elif force and os.path.exists(raxml_work_dir):
        rmtree(raxml_work_dir)
        os.mkdir(raxml_work_dir)

    gubbins_snp_phylip = os.path.join(gubbins_work_dir, "gubbins.filtered_polymorphic_sites.phylip")
    cmd = f"{ex_raxml_ng} --all --msa {gubbins_snp_phylip} --prefix {raxml_prefix} --model GTR+FO+IO" \
          f" --bs-trees autoMRE --threads 4"
    log_message = f"Command used : \n {cmd}\n"
    # launch
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    log_file_path = os.path.join(raxml_work_dir, "logRaXML-ng.txt")
    log_process_with_output_file(process, log_message, log_file_path)

    ###################################
    # Write Analysis
    section_header(f'Write Final Analysis {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    all_report_file = os.path.join(work_dir, "all_report.tsv")
    with open(all_report_file, "w") as all_report:
        writer_report_all = writer(all_report, delimiter='\t')
        for key, sample_dict in detection_result_dict.items():
            if key == "Sample":
                writer_report_all.writerow(["Sample"] + [*sample_dict])
                pass
            else:
                write_list = [key]
                for sample_name in sample_list:
                    write_list.append(sample_dict[sample_name])
                writer_report_all.writerow(write_list)
    lvl_report_file = os.path.join(cluster_dir, "lvl_report.tsv")
    with open(lvl_report_file, "w") as lvl_report:
        writer_report = writer(lvl_report, delimiter='\t')
        writer_report.writerow(["Name level", "Groups"])
        for name_lvl, lvl_dict in level_gg_dict.items():
            writer_report.writerow([name_lvl, ";".join(str(v) for v in lvl_dict["list_sample"])])


def mafft(output_dir_msa, output_fasta_file, locus_name, output_aln_file):
    """
    This function call MAFFT tool
    :param output_dir_msa: The MSA directory path
    :param output_fasta_file: The fasta file path
    :param locus_name: The locus name
    :param output_aln_file: The MSA output directory path
    """
    ex_mafft = shutil.which("mafft")
    cmd = f"{ex_mafft} {output_fasta_file} > {output_aln_file}"
    # log_message = f"Command used : \n {cmd}\n"
    # launch
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    # log_file_path = os.path.join(output_dir_msa, f"logMAFFT_{locus_name}.txt")
    # log_process_with_output_file(process, log_message, log_file_path)
    while True:
        if process.poll() is not None:
            break

