#!/usr/bin/env python3
"""
Copyright 2020 Aurélien BIRER (abirer36@gmail.com)
https://github.com/Nilad/CGST.git

This script is the build database step of CGST.

This file is part of CGST. CGST is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. CGST is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with CGST. If
not, see <http://www.gnu.org/licenses/>.
"""
import shutil
import os
from subprocess import Popen, PIPE, STDOUT

import pandas


from cgst.log import log, section_header, explanation, log_process
from cgst.utils import available_species_mentalist


def build_database_mentalist(cgst_database_path, species_full, threads):
    """
    This function build MentaLiST cgMLST database
    :param cgst_database_path: The CGST database path
    :param species_full: The scientific name of the strain
    :param threads: The number of threads to allocate
    """
    exe = shutil.which("mentalist")
    species_full = species_full.strip().rstrip()
    species = species_full.split(" ")[0][0].lower() + species_full.split(" ")[1].lower()

    if not os.path.exists(cgst_database_path):
        os.mkdir(cgst_database_path)
    cgmlst_database_path = os.path.join(cgst_database_path, "cgMLST")
    if not os.path.exists(cgmlst_database_path):
        os.mkdir(cgmlst_database_path)

    cgmlst_species_database_path = os.path.join(cgmlst_database_path, species)
    if not os.path.exists(cgmlst_species_database_path):
        os.mkdir(cgmlst_species_database_path)
    else:
        log(f"cgMLST for {species_full} already present")

    #####
    # Check if cgmlst in CNR repo
    url = 'https://raw.githubusercontent.com/CNRResistanceAntibiotic/core_genomes/master/reference.csv'
    df = pandas.read_csv(url, error_bad_lines=False)

    url_folder = ""
    for i, row in enumerate(df["Species"]):
        if row == species_full:
            log("Species {0} found in CNR GitHub".format(species_full))
            url_folder = "https://github.com/CNRResistanceAntibiotic/core_genomes/trunk/{0}/{1}"\
                .format(df["Name"][i], df["Sub-folder"][i])
        else:
            continue

    if url_folder:
        cmd = f"svn export {url_folder} {cgmlst_species_database_path} --force"
        print(cmd)
        # launch
        process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, executable='/bin/bash')
        while True:
            if process.stdout is not None:
                output = process.stdout.readline().decode("utf-8").rstrip()
            if process.stderr is not None:
                error = process.stderr.readline().decode("utf-8").rstrip()
            if output == '' and process.poll() is not None:
                break
    else:
        id_species_list = available_species_mentalist(species_full)
        if not id_species_list:
            log()
            log(f"Any ID found for the species : {species_full}")
            log()
            exit()
        else:
            for id_species in id_species_list:
                log(f"ID ({id_species.lower()}) found for {species_full} ")
        ###################################
        # Run MentaLiST Build Species Database
        section_header(f'MentaLiST Build database cgMLST for species : {species_full}')
        explanation('Download and install cgMLST for a species using MentaLiST function')
        # prepare

        fasta_database_name = f"{species}_cgmlst_fasta"
        fasta_database_path = os.path.join(cgmlst_species_database_path, fasta_database_name)

        database_name = f"{species}_cgmlst.db"
        database_path = os.path.join(cgmlst_species_database_path, database_name)

        cmd = f"{exe} download_cgmlst -k 31 -o {fasta_database_path} -s {id_species_list[0].lower()} --db {database_path}" \
              f" --threads {threads}"
        # launch
        process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, executable='/bin/bash')
        log_process(process, "")


def build_database_ariba(cgst_database_path, species_full):
    """
    This function build ARIBA MLST database
    :param cgst_database_path: The CGST database path
    :param species_full: The scientific name of the strain
    """
    exe = shutil.which("ariba")
    species_full = species_full.strip().rstrip()
    species = species_full.split(" ")[0][0].lower() + species_full.split(" ")[1].lower()
    id_species_list = available_species_ariba(species_full.lower())
    if not id_species_list:
        log()
        log(f"Any ID found for the species : {species_full}")
        log()
        exit()
    else:
        for id_species in id_species_list:
            log(f"ID ({id_species}) found for {species_full} ")
    ###################################
    # Run Ariba Build Species Database
    section_header(f'Ariba Build database for species : {species_full}')
    explanation('Download and install MLST for a species using Ariba function')
    # prepare

    if not os.path.exists(cgst_database_path):
        os.mkdir(cgst_database_path)
    mlst_database_path = os.path.join(cgst_database_path, "ariba")
    if not os.path.exists(mlst_database_path):
        os.mkdir(mlst_database_path)
    mlst_species_database_path = os.path.join(mlst_database_path, species)
    if not os.path.exists(mlst_species_database_path):
        os.mkdir(mlst_species_database_path)

    for id_species in id_species_list:
        id_species_name = id_species.lower().replace(" ", "-").replace("#", "")
        out_pubmlst_get_path = os.path.join(mlst_species_database_path, f"mlst_{id_species_name}")

        if os.path.exists(out_pubmlst_get_path):
            shutil.rmtree(out_pubmlst_get_path)

        cmd = f"{exe} pubmlstget \"{id_species}\" {out_pubmlst_get_path}"
        # launch
        process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, executable='/bin/bash')
        log_process(process, "")


def available_species_ariba(species):
    """
    This function display ARIBA schema MLST available
    :param species: The scientific name of the strain
    :return: the id of schema
    """
    exe = shutil.which("ariba")
    ###################################
    # Run Ariba Available Species
    section_header('Ariba Available MLST Schemas')
    explanation('Display the MLST by Ariba function that search in the pubmlst database')
    # prepare
    cmd = f"{exe} pubmlstspecies"
    # launch
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    id_list = log_process(process, species)
    return id_list


def build_database(cgst_database_path, species_full, threads):
    """
    This function manage the building database functions
    :param cgst_database_path: The CGST database path
    :param species_full: The scientific name of the strain
    :param threads: The number of threads to allocate
    """
    build_database_mentalist(cgst_database_path, species_full, threads)
    build_database_ariba(cgst_database_path, species_full)
