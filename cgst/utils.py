#!/usr/bin/env python3
"""
Copyright 2020 Aurélien BIRER (abirer36@gmail.com)
https://github.com/Nilad/CGST.git

This script is the utils file part of CGST.

This file is part of CGST. CGST is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. CGST is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with CGST. If
not, see <http://www.gnu.org/licenses/>.
"""
import csv
import shutil
from subprocess import Popen, PIPE, STDOUT

import pandas

from cgst.log import section_header, explanation, log_process, tool_output_log


def mentalist_path_and_version(mentalist_path):
    """
    This function manage the MentaLiST executable path
    :param mentalist_path: The path of the executable MentaLiST
    :return: the MentaLiST path, the version of MentaLiST, status value
    """
    found_mentalist_path = shutil.which(mentalist_path)
    if found_mentalist_path is None:
        return mentalist_path, '', 'not found'
    command = [found_mentalist_path]
    process = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT, executable='/bin/bash')
    out, _ = process.communicate()
    out = out.decode().lower()
    if 'mentalist' not in out or '[-v]' not in out:
        return found_mentalist_path, '-', 'bad'
    return found_mentalist_path, mentalist_version(found_mentalist_path), 'good'


def mentalist_version(tool_path):
    """
    This function return the MentaLiST version
    :param tool_path: The MentaLiST executable path
    :return: the version of MentaLiST
    """
    command = "{0} -v".format(tool_path)
    process = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT, executable='/bin/bash')
    out, _ = process.communicate()
    out = out.decode().lower().strip()
    if len(out) == 0 or '.' not in out:
        return '-'
    if out.startswith('m'):
        return out[9:]
    else:
        return out


def read_output_mentalist(file_path, strain_name):
    """
    This function read and return a dictionary of the final output of MentaLiST
    :param file_path: The final output of MentaLiST
    :param strain_name: strain name
    :return: a dictionary of the MentaLiST output
    """
    detection_result_dict = {}
    with open(file_path, "r") as detection_file:
        reader = csv.DictReader(detection_file, delimiter="\t")
        for row in reader:
            for key, value in row.items():
                if key in detection_result_dict:
                    dict_value = detection_result_dict[key]
                    dict_value[strain_name] = value
                else:
                    detection_result_dict[key] = {strain_name: value}
    return detection_result_dict


def available_species_mentalist(species):
    """
    This function display the species available by MentaLiST
    :param species: Name of species to get the id
    :return: The list of species id
    """
    exe = shutil.which("mentalist")
    ###################################
    # Run MentaLiST Available Species
    section_header('MentaLiST Available cgMLST Genome')
    explanation('Display the cgMLST by MentaLiST function that search in the cgmlst.org database')
    # prepare
    cmd = f"{exe} list_cgmlst"
    # launch
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, executable='/bin/bash')
    id_list = log_process(process, species)

    ###################################
    section_header('CNR "Résistance aux antibiotiques" Available cgMLST Genome')
    explanation('Display the cgMLST that search in the core_genomes github repository')

    url = 'https://raw.githubusercontent.com/CNRResistanceAntibiotic/core_genomes/master/reference.csv'
    df = pandas.read_csv(url, error_bad_lines=False)

    for row in df["Species"]:
        tool_output_log(row)

    return id_list
