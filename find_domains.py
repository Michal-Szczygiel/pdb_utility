#!/usr/bin/python

import re
import argparse
import os
import glob
from collections import Counter
from bioinf_tools import *


# Funkcja pomocnicza - formatująca
def format_output(ligand : Ligand, residue : Residue, chain_name : str, pfam_domain : Pfam_domain) -> str:
    output = f'{ligand}'

    if residue != None:
        output += f' {residue}'

    if chain_name != None:
        output += ' {: <8}'.format(chain_name)

    if pfam_domain != None:
        output += ' Pfam domain -> {: <10} {: <23} ({: <5} - {: >5})'.format(pfam_domain.pfam_accession, pfam_domain.pfam_id, pfam_domain.start, pfam_domain.end)

    return output


# -------------------------------------------------------------------------------------------------------------


CONFIG_FILE = None                                          # Ścieżka do pliku konfiguracyjnego

PDB_DIRECTORY_PATH = ''                                     # Ścieżka do katalogu z plikami PDB
XML_DIRECTORY_PATH = ''                                     # Ścieżka do katalogu z plikami XML

LOG_FILE = None                                             # Plik wynikowy - generalne podsumowanie
GENERAL_STATISTICS_FILE = None                              # Plik wynikowy - statystyka wszystkich znalezionych domen
MATCHED_STATISTICS_FILE = None                              # Plik wynikowy - statystyka dopasowanych domen
WHERE_LIGAND_FILE = None                                    # Plik wynikowy - lokalizacja ligandów w domenach

LIGANDS = None                                              # Lista domen uwzględnianych w przeszukiwaniu

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--config', type=str)
    arg_parser.add_argument('--batch_s', type=int)
    arg_parser.add_argument('--verbose', action='store_true')
    args = arg_parser.parse_args()

    CONFIG_FILE = args.config

    try:
        # ----------------------------------------------------------------------------------------------------------------------
        # Parsowanie pliku konfiguracyjnego. Wczytywanie ścieżek niezbędnych do dalszej pracy skryptu oraz nazw ligandów.
        # ----------------------------------------------------------------------------------------------------------------------

        config_file_content = open(CONFIG_FILE, 'r').read()

        pdb_dir_parsed = re.findall(r'pdb_directory:.*\n', config_file_content)
        if len(pdb_dir_parsed) == 1:
            PDB_DIRECTORY_PATH = pdb_dir_parsed[0].split(':')[1].strip()

        xml_dir_parsed = re.findall(r'xml_directory:.*\n', config_file_content)
        if len(xml_dir_parsed) == 1:
            XML_DIRECTORY_PATH = xml_dir_parsed[0].split(':')[1].strip()

        ligands_parsed = re.findall(r'ligands: *\{[a-zA-Z0-9 \n]*\}', config_file_content)
        if len(ligands_parsed) == 1:
            LIGANDS = ligands_parsed[0].split(':')[1].replace('{', '').replace('}', '').split()

        log_file_path = re.findall(r'log_file:.*\n', config_file_content)
        if len(log_file_path) == 1:
            LOG_FILE = log_file_path[0].split(':')[1].strip()
            
        general_statistics_file = re.findall(r'general_statistics_file:.*\n', config_file_content)
        if len(log_file_path) == 1:
            GENERAL_STATISTICS_FILE = general_statistics_file[0].split(':')[1].strip()

        matched_statistics_file = re.findall(r'matched_statistics_file:.*\n', config_file_content)
        if len(log_file_path) == 1:
            MATCHED_STATISTICS_FILE = matched_statistics_file[0].split(':')[1].strip()

        where_ligand_file = re.findall(r'where_ligand_file:.*\n', config_file_content)
        if len(log_file_path) == 1:
            WHERE_LIGAND_FILE = where_ligand_file[0].split(':')[1].strip()

    except:
        print(f'Config file: {CONFIG_FILE} not found.')

    if os.path.isdir(PDB_DIRECTORY_PATH) and os.path.isdir(XML_DIRECTORY_PATH) and LIGANDS:
        # Listing wszystkich wczytanych danych wczystancyh z pliku konfiguracyjnego.

        print(f'PDB_DIRECTORY_PATH: {PDB_DIRECTORY_PATH}')
        print(f'XML_DIRECTORY_PATH: {XML_DIRECTORY_PATH}')
        print(f'LOG_FILE: {LOG_FILE}')
        print(f'GENERAL_STATISTICS_FILE: {GENERAL_STATISTICS_FILE}')
        print(f'MATCHED_STATISTICS_FILE: {MATCHED_STATISTICS_FILE}')
        print(f'WHERE_LIGAND_FILE: {WHERE_LIGAND_FILE}')
        print(f'LIGANDS: {LIGANDS}')
        
        PDB_FILES_LIST =  glob.glob(f'{PDB_DIRECTORY_PATH}/*.pdb')

        BATCH_SIZE = args.batch_s
        if BATCH_SIZE == None or BATCH_SIZE > len(PDB_FILES_LIST):
            BATCH_SIZE = len(PDB_FILES_LIST)

        if LOG_FILE != None:
            result_file = open(LOG_FILE, 'w')

        skipped_pdb_files = []

        all_domains_statistics : list[Pfam_domain] = []
        matched_domains_statistics = []
        matched_domains_statistics_detailed = []
        where_ligand : list[tuple[str, str, int, str, int, str]] = []             #pfam_id, ligand_name, ligand_seq_number, res_name, res_seq_number, chain_name

        # ----------------------------------------------------------------------------------------------------------------------
        # Główna pętla skryptu analizującego. Poniższe procedury zostają wykonane dla każdej struktury ze zbioru.
        # ----------------------------------------------------------------------------------------------------------------------
        for (counter, pdb_file_path) in enumerate(PDB_FILES_LIST):
            if counter > BATCH_SIZE - 1:
                break

            # Wywołania funkcji z modułu bioinf_tools
            molecule = parse_PDB_file(pdb_file_path, LIGANDS)
            molecule.collect_data_from_Pfam(XML_DIRECTORY_PATH)

            # Sprawdzenie poprawności wczytanych danych, kontynuacja tylko w przypadku braku błędów.
            if molecule.data_correctnes_flag == True:
                all_domains = molecule.PFAM_DOMAINS
                ligands_environment = molecule.get_ligands_environment()
                ligands_binding_domains = molecule.get_ligands_binding_domains()
                ligands_binding_domains_detailed = molecule.get_ligands_binding_domains_detailed()
                
                # Wykonanie tylko jeżeli istnieje potrzeba zapisu danych o składzie aminokwasowym miejsc wiązania wskazanych ligandów w kontekście domen.
                if WHERE_LIGAND_FILE != None:
                    for (ligand, residue, chain_name, pfam_domain) in ligands_environment:
                        if pfam_domain != None:
                            where_ligand.append( (
                                pfam_domain.pfam_id,
                                ligand.residue_name,
                                ligand.residue_sequence_number,
                                residue.residue_name,
                                residue.residue_sequence_number,
                                chain_name
                            ) )
                
                all_domains_statistics += [(domain.pfam_accession, domain.pfam_id) for domain in all_domains]
                matched_domains_statistics += [(domain.pfam_accession, domain.pfam_id, ligand.residue_name) for (ligand, domain) in ligands_binding_domains]
                matched_domains_statistics_detailed += [(domain.pfam_accession, domain.pfam_id, ligand.residue_name) for (ligand, domain) in ligands_binding_domains_detailed]

                # Wykonanie jeżeli istnieje potrzeba zapisu szczegółowych danych o przetwarzanych strukturach.
                if LOG_FILE != None:
                    result_file.write(f'PDB_ID: {molecule.structure_PDB_ID}\n\nDOMAIN ARCHITECTURE:\n')

                    if all_domains:
                        for domain in all_domains:
                            result_file.write(f'{domain}\n')
                    else:
                        result_file.write(f'< no domains found >\n')

                    if args.verbose:
                        result_file.write('\nLIGANDS ENVIRONMENT:\n')

                        if ligands_environment:
                            for (ligand, residue, chain_name, pfam_domain) in ligands_environment:
                                result_file.write(f'{format_output(ligand, residue, chain_name, pfam_domain)}\n')
                        else:
                            result_file.write(f'< no ligands found >\n')
                    else:
                        result_file.write('\nMATCHED DOMAINS:\n')

                        if ligands_binding_domains:
                            for (ligand, domain) in ligands_binding_domains:
                                result_file.write(f'{domain}    {ligand}\n')
                        else:
                            result_file.write(f'< no domains found >\n')

                    result_file.write('\n------------------------------------------------------------------------------------------------------------------------------------------------------\n\n')

                print(Colors.GREEN + 'Done: {: <32} ({})'.format(pdb_file_path, counter + 1) + Colors.END)

            else:
                skipped_pdb_files.append(pdb_file_path)
                print(Colors.RED + 'Skipped: {: <32} ({})'.format(pdb_file_path, counter) + Colors.END)

        if skipped_pdb_files and LOG_FILE != None:
            result_file.write('SKIPPED PDB FILES:\n')

            for skipped in skipped_pdb_files:
                result_file.write(f'{skipped}\n')

        if LOG_FILE != None:
            result_file.close()

        print(Colors.YELLOW + f'Processed {BATCH_SIZE} structures ({len(skipped_pdb_files)} skipped)' + Colors.END)

        # -------------------------------------------------------------------------------------------------------------

        # Wykonanie jeżeli istnieje potrzeba zapisu generalnej statystyki dotyczącej wszystkich znalezionych domen białkowych. 
        if GENERAL_STATISTICS_FILE != None:
            general_statistics = open(GENERAL_STATISTICS_FILE, 'w')
            for ((pfam_accession, pfam_id), occurrences) in Counter(all_domains_statistics).most_common():
                general_statistics.write('Pfam accession: {: <15} Pfam ID: {: <23} Occurrences: {: <10}\n'.format(pfam_accession, pfam_id, occurrences))

            general_statistics.close()

        # Wykonanie jeżeli istnieje potrzeba zapisu statystyki dotyczącej wszystkich domen w obrębie których znaleziono wskazane ligandy. 
        if MATCHED_STATISTICS_FILE != None:
            alfabetic = lambda record: ' '.join([record[0][1], record[0][2]])

            matched_statistics = open(MATCHED_STATISTICS_FILE, 'w')
            hit_number = [number for ((_, _, _), number) in sorted(Counter(matched_domains_statistics_detailed).items(), key=alfabetic)]
            for (((pfam_accession, pfam_id, ligand_name), occurrences), number) in zip(sorted(Counter(matched_domains_statistics).items(), key=alfabetic), hit_number):
                matched_statistics.write('Pfam accession: {: <15} Pfam ID: {: <23} Ligand name: {: <10} Occurrences: {: <10} Molecule/domain: {: <10}\n'
                    .format(pfam_accession, pfam_id, ligand_name, occurrences, round(number / occurrences, 2)))

            matched_statistics.close()
        
        # Wykonanie tylko jeżeli istnieje potrzeba zapisu danych o składzie aminokwasowym miejsc wiązania wskazanych ligandów w kontekście domen.
        if WHERE_LIGAND_FILE != None:
            def first(record):
                return record[0]

                                        #pfam_id, ligand, residues
            residue_statistics : list[tuple[str, str, list[str]]] = []
            residue_statistics_quick_compare : list[tuple[str, str]] = []
            where_ligand.sort(key=first)
                #pfam_id, ligand_name, ligand_seq_number, res_name, res_seq_number, chain_name
            for (pfam_id, ligand_name, _, res_name, _, _) in where_ligand:
                try:
                    index = residue_statistics_quick_compare.index( (pfam_id, ligand_name) )
                    if res_name not in residue_statistics[index][2]:
                        residue_statistics[index][2].append(res_name)
                except:
                    residue_statistics_quick_compare.append( (pfam_id, ligand_name) )
                    residue_statistics.append( (pfam_id, ligand_name, [res_name]) )

            where_ligand_statistics = open(WHERE_LIGAND_FILE, 'w')

            for (pfam_id, ligand_name, res_list) in residue_statistics:
                where_ligand_statistics.write('{: <23} {: <4}    [{}]\n'.format(pfam_id, ligand_name, ' '.join(sorted(res_list))))

            where_ligand_statistics.close()

