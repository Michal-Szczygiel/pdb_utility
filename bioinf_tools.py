#!/usr/bin/python

try:
    import platform
    import xml.dom.minidom as mnd
    import urllib3
except:
    print('Necessary libraries not found :/')

if platform.system() == 'Linux':
    class Colors:
        RED = '\33[31m'
        GREEN = '\33[32m'
        BLUE = '\33[34m'
        YELLOW = '\33[33m'
        END = '\033[0m'

else:
    class Colors:
        RED = ''
        GREEN = ''
        BLUE = ''
        YELLOW = ''
        END = ''


# ------------------------------------------------------------------------------------------------------


class Unp_record:
    """Klasa służąca do przechowywania odniesień do bazy danych UniProt znalezionych w plikach PDB."""

    def __init__(self, chain_name, unp_accession, unp_id):
        try:
            self.chain_name : str = str(chain_name).strip()
            self.unp_accession : str = str(unp_accession).strip()
            self.unp_id : str = str(unp_id).strip()
        except:
            print(Colors.RED + 'UNP_record: ValueError!' + Colors.END)
            raise ValueError

    def __repr__(self):
        return f'UNP Record -> chain: {self.chain_name}, UNP accession: {self.unp_accession}, ID code: {self.unp_id}'


# ------------------------------------------------------------------------------------------------------


class Pfam_domain:
    """Klasa służąca do przechowywania informacji o domenach białkowych (Pfam)."""

    def __init__(self, chain_name, pfam_accession, pfam_id, start, end):
        try:
            self.chain_name : str = str(chain_name).strip()
            self.pfam_accession : str = str(pfam_accession).strip()
            self.pfam_id : str = str(pfam_id).strip()
            self.start : int = int(start)
            self.end : int = int(end)
        except:
            print(Colors.RED + 'Pfam_domain: ValueError!' + Colors.END)
            raise ValueError

    def __repr__(self):
        return 'Pfam_domain -> {: <10} {: <23} {: <8} ({: <5} - {: >5})'.format(self.pfam_accession, self.pfam_id, self.chain_name, self.start, self.end)


# ------------------------------------------------------------------------------------------------------


class Atom:
    """Klasa służąca do przechowywania informacji o pojedynczych atomach pobranych ze struktury białka (nazwa, symbol, położenie w R^3 [A])."""

    def __init__(self, atom_name, element_symbol, x_coordinate, y_coordinate, z_coordinate):
        try:
            self.atom_name : str = str(atom_name).strip()
            self.element_symbol : str = str(element_symbol).strip()
            self.x_coordinate : float = float(x_coordinate)
            self.y_coordinate : float = float(y_coordinate)
            self.z_coordinate : float = float(z_coordinate)
        except:
            print(Colors.RED + 'Atom: ValueError!' + Colors.END)
            raise ValueError


# ------------------------------------------------------------------------------------------------------


class Termination_symbol:
    """Znacznik końca łańcucha aminokwasowego."""

    def __init__(self, residue_name, residue_sequence_number):
        try:
            self.residue_name : str = str(residue_name).strip()
            self.residue_sequence_number : int = int(residue_sequence_number)
        except:
            print(Colors.RED + 'Termination_symbol: ValueError!' + Colors.END)
            raise ValueError


# ------------------------------------------------------------------------------------------------------


class Residue:
    """Klasa zawierająca grupę wczytanych atomów (klasa: Atom), stanowiąca reprezentację reszty aminokwasowej."""

    def __init__(self, residue_name, residue_sequence_number):
        try:
            self.residue_name : str = str(residue_name).strip()
            self.residue_sequence_number : int = int(residue_sequence_number)
            self.ATOMS : list[Atom] = []
        except:
            print(Colors.RED + 'Residue: ValueError!' + Colors.END)
            raise ValueError

    def __repr__(self):
        return 'Residue -> {} {: <4}'.format(self.residue_name, self.residue_sequence_number)

    def push_atom(self, next_atom : Atom):
        self.ATOMS.append(next_atom)


# ------------------------------------------------------------------------------------------------------


class Ligand:
    """Klasa zawierająca grupę wczytanych atomów (klasa: Atom), stanowiąca reprezentację cząsteczki liganda."""

    def __init__(self, residue_name, residue_sequence_number):
        try:
            self.residue_name : str = str(residue_name).strip()
            self.residue_sequence_number : int = int(residue_sequence_number)
            self.ATOMS : list[Atom] = []
        except:
            print(Colors.RED + 'Ligand: ValueError!' + Colors.END)
            raise ValueError

    def __repr__(self):
        return 'Ligand -> {} {: <8}'.format(self.residue_name, self.residue_sequence_number)

    def push_atom(self, next_atom : Atom):
        self.ATOMS.append(next_atom)


# ------------------------------------------------------------------------------------------------------


class Chain:
    """Klasa zawierająca grupę wczytanych reszt aminokwasowych (klasa: Residue), stanowiąca reprezentację łańcucha białkowego."""

    def __init__(self, chain_name):
        try:
            self.chain_name : str = str(chain_name).strip()
            self.RESIDUES : list[Residue] = []
        except:
            print(Colors.RED + 'Chain: ValueError!' + Colors.END)
            raise ValueError

    def __repr__(self):
        chain_arrow = f'{Colors.GREEN}CHAIN: {self.chain_name} --> {Colors.END}'
        residues_names = ''

        for i in range(len(self.RESIDUES)):
            residues_names += self.RESIDUES[i].residue_name
            if ((i + 1) % 10 == 0):
                residues_names += '\n'
                residues_names += ''.ljust(len(chain_arrow) - 9)
            else:
                residues_names += ' '

        return f'{chain_arrow}{residues_names}'

    def push_residue(self, next_residue : Residue):
        self.RESIDUES.append(next_residue)


# ------------------------------------------------------------------------------------------------------


class Structure:
    """Klasa zawierająca grupę wczytanych łańcuchów aminokwasowych (klasa: Chain), stanowiąca reprezentację całego białka."""

    def __init__(self, structure_PDB_ID):
        try:
            self.structure_PDB_ID : str = str(structure_PDB_ID).strip()
            self.data_correctnes_flag : bool = True
            self.CHAINS : list[Chain] = []
            self.LIGANDS : list[Ligand] = []   
            self.UNP_RECORDS : list[Unp_record] = []
            self.PFAM_DOMAINS : list[Pfam_domain] = []
        except:
            print(Colors.RED + 'Structure: ValueError!' + Colors.END)
            raise ValueError

    def collect_data_from_Pfam(self, backup_folder):
        """Funkcja pobierająca informacje o architekturach domenowych z bazy danych Pfam (pod warunkiem, że danych tych w postaci plików *.xml nie ma w lokalnym repozytorium)."""

        try:
            if self.data_correctnes_flag == True:
                http = urllib3.PoolManager()
                for unp_record in self.UNP_RECORDS:
                    try:
                        xml_file = open(f'{backup_folder}/{self.structure_PDB_ID}_{unp_record.chain_name}_{unp_record.unp_accession}.xml', 'r')
                        decoded_content = xml_file.read()
                        xml_file.close()
                    except:
                        resource = http.request('GET', f'http://pfam.xfam.org/protein?output=xml&acc={unp_record.unp_accession}')
                        decoded_content = resource.data.decode('utf-8')
                        xml_file = open(f'{backup_folder}/{self.structure_PDB_ID}_{unp_record.chain_name}_{unp_record.unp_accession}.xml', 'w')
                        xml_file.write(decoded_content)
                        xml_file.close()

                    parser = mnd.parseString(decoded_content)
                    matches = parser.getElementsByTagName('match')
                    for match in matches:
                        for child in match.childNodes:
                            if child.nodeType == child.ELEMENT_NODE and child.tagName == 'location':
                                pfam_record = Pfam_domain(
                                    unp_record.chain_name,
                                    match.getAttribute('accession'), 
                                    match.getAttribute('id'),
                                    child.getAttribute('start'),
                                    child.getAttribute('end'))
                                self.PFAM_DOMAINS.append(pfam_record)

        except:
            self.data_correctnes_flag = False
            print(Colors.RED + f'Could not get domain information for structure: {self.structure_PDB_ID}' + Colors.END)
    
    def get_ligands_environment(self) -> list[tuple[Ligand, Residue, str, Pfam_domain]]:
        """Funkcja zwracająca listę znalezionych miejsc kontaktu ligand - białko."""

        distance = lambda atom_a, atom_b: (
            (atom_a.x_coordinate - atom_b.x_coordinate) ** 2 
            + (atom_a.y_coordinate - atom_b.y_coordinate) ** 2 
            + (atom_a.z_coordinate - atom_b.z_coordinate) ** 2
        ) ** 0.5

        comparison_tuple = lambda ligand, residue, chain: (
            ligand.residue_name,
            ligand.residue_sequence_number,
            residue.residue_name,
            residue.residue_sequence_number,
            chain.chain_name
        )

        fits = lambda domain, residue: domain.start <= residue.residue_sequence_number and domain.end >= residue.residue_sequence_number

        ligands_environment : list[tuple[Ligand, Residue, str, Pfam_domain]] = []
        ligands_environment_quick_compare = []

        for ligand in self.LIGANDS:
            for chain in self.CHAINS:
                for residue in chain.RESIDUES:
                    for atom_a in [atom for atom in residue.ATOMS if atom.atom_name[0] != 'H']:
                        for atom_b in [atom for atom in ligand.ATOMS if atom.atom_name[0] != 'H']:
                            if distance(atom_a, atom_b) < 5.0:
                                comp_tuple = comparison_tuple(ligand, residue, chain)
                                if comp_tuple not in ligands_environment_quick_compare:
                                    ligands_environment_quick_compare.append(comp_tuple)

                                    found_domain = False
                                    for pfam_domain in self.PFAM_DOMAINS:
                                        if pfam_domain.chain_name == chain.chain_name and fits(pfam_domain, residue):
                                            ligands_environment.append( (ligand, residue, chain.chain_name, pfam_domain) )
                                            found_domain = True

                                    if found_domain == False:
                                        ligands_environment.append( (ligand, residue, chain.chain_name, None) )

        for ligand in self.LIGANDS:
            if (ligand.residue_name, ligand.residue_sequence_number) not in [(m_ligand.residue_name, m_ligand.residue_sequence_number) for (m_ligand, _, _, _) in ligands_environment]:
                ligands_environment.append( (ligand, None, None, None) )

        return ligands_environment

    def get_ligands_binding_domains(self) -> list[tuple[Ligand, Pfam_domain]]:
        """Funkcja zwracająca domeny białkowe wiążące wskazane ligandy (nie uwzględniono wielokrotnych miejsc wiązania dla tego samego liganda w obrębie tej samej domeny)."""

        ligands_environment = self.get_ligands_environment()
        matched_domains : list[tuple[Ligand, Pfam_domain]] = []
        matched_domains_quick_compare = []

        for (ligand, _, _, pfam_domain) in ligands_environment:
            if pfam_domain != None:
                domain = (
                    ligand.residue_name,
                    pfam_domain.chain_name, 
                    pfam_domain.pfam_accession, 
                    pfam_domain.pfam_id, 
                    pfam_domain.start, 
                    pfam_domain.end
                )
                if domain not in matched_domains_quick_compare:
                    matched_domains_quick_compare.append(domain)
                    matched_domains.append( (ligand, pfam_domain) )

        return matched_domains

    def get_ligands_binding_domains_detailed(self) -> list[tuple[Ligand, Pfam_domain]]:
        """Funkcja zwracająca domeny białkowe wiążące wskazane ligandy (uwzględniono wielokrotne miejsca wiązania dla tego samego liganda w obrębie tej samej domeny)."""

        ligands_environment = self.get_ligands_environment()
        matched_domains : list[tuple[Ligand, Pfam_domain]] = []
        matched_domains_quick_compare = []

        for (ligand, _, _, pfam_domain) in ligands_environment:
            if pfam_domain != None:
                domain = (
                    ligand.residue_name,
                    ligand.residue_sequence_number,
                    pfam_domain.chain_name, 
                    pfam_domain.pfam_accession, 
                    pfam_domain.pfam_id, 
                    pfam_domain.start, 
                    pfam_domain.end
                )
                if domain not in matched_domains_quick_compare:
                    matched_domains_quick_compare.append(domain)
                    matched_domains.append( (ligand, pfam_domain) )

        return matched_domains


# ------------------------------------------------------------------------------------------------------   


def parse_PDB_file(pdb_file_path, ligands : list[str] = []) -> Structure:
    """Funkcja parsująca pliki PDB. Zwraca wczytaną cząsteczkę w postaci obiektu klasy Structure."""

    try:
        structure = Structure(pdb_file_path[-8 : -4].upper())
        pdb_file = open(pdb_file_path, 'r')
        line_counter = -1

        try:
            for DATA in pdb_file:
                line_counter += 1
                RECORD_TYPE = DATA[0 : 6]
                LIGAND_TYPE = DATA[17 : 20]
                DATABASE = DATA[26 : 32]

                if RECORD_TYPE == 'HETATM' and LIGAND_TYPE.strip() in ligands:
                    ATOM_NAME = DATA[12 : 16]
                    ELEM_SYMBOL = DATA[76 : 78]
                    X_COORD = DATA[30 : 38]
                    Y_COORD = DATA[38 : 46]
                    Z_COORD = DATA[46 : 54]

                    RES_NAME = DATA[17 : 20]
                    RES_SEQ_NUM = DATA[22 : 26]
                    
                    next_atom = Atom(ATOM_NAME, ELEM_SYMBOL, X_COORD, Y_COORD, Z_COORD)

                    if structure.LIGANDS:
                        if structure.LIGANDS[-1].residue_sequence_number == int(RES_SEQ_NUM):
                            structure.LIGANDS[-1].push_atom(next_atom)
                        else:
                            next_ligand = Ligand(RES_NAME, RES_SEQ_NUM)
                            structure.LIGANDS.append(next_ligand)
                            structure.LIGANDS[-1].push_atom(next_atom)
                    else:
                        next_ligand = Ligand(RES_NAME, RES_SEQ_NUM)
                        structure.LIGANDS.append(next_ligand)
                        structure.LIGANDS[-1].push_atom(next_atom)

                elif RECORD_TYPE == 'ATOM  ' or RECORD_TYPE == 'HETATM':
                    ATOM_NAME = DATA[12 : 16]
                    ELEM_SYMBOL = DATA[76 : 78]
                    X_COORD = DATA[30 : 38]
                    Y_COORD = DATA[38 : 46]
                    Z_COORD = DATA[46 : 54]

                    CHAIN_IDENT = DATA[21 : 22]
                    RES_NAME = DATA[17 : 20]
                    RES_SEQ_NUM = DATA[22 : 26]

                    next_atom = Atom(ATOM_NAME, ELEM_SYMBOL, X_COORD, Y_COORD, Z_COORD)
                    if structure.CHAINS:
                        if structure.CHAINS[-1].chain_name == str(CHAIN_IDENT).strip():
                            if structure.CHAINS[-1].RESIDUES[-1].residue_sequence_number == int(RES_SEQ_NUM):
                                structure.CHAINS[-1].RESIDUES[-1].push_atom(next_atom)
                            else:
                                next_residue = Residue(RES_NAME, RES_SEQ_NUM)
                                structure.CHAINS[-1].push_residue(next_residue)
                                structure.CHAINS[-1].RESIDUES[-1].push_atom(next_atom)
                        else:
                            next_chain = Chain(CHAIN_IDENT)
                            next_residue = Residue(RES_NAME, RES_SEQ_NUM)
                            structure.CHAINS.append(next_chain)
                            structure.CHAINS[-1].push_residue(next_residue)
                            structure.CHAINS[-1].RESIDUES[-1].push_atom(next_atom)
                    else:
                        next_chain = Chain(CHAIN_IDENT)
                        next_residue = Residue(RES_NAME, RES_SEQ_NUM)
                        structure.CHAINS.append(next_chain)
                        structure.CHAINS[-1].push_residue(next_residue)
                        structure.CHAINS[-1].RESIDUES[-1].push_atom(next_atom)

                elif RECORD_TYPE == 'TER   ':
                    RES_NAME = DATA[17 : 20]
                    RES_SEQ_NUM = DATA[22 : 26]

                    structure.CHAINS[-1].push_residue(Termination_symbol(RES_NAME, RES_SEQ_NUM))

                elif RECORD_TYPE == 'DBREF ' and DATABASE == 'UNP   ':
                    CHAIN_IDENT = DATA[12 : 13]
                    DB_ACCESSION = DATA[33 : 41]
                    ID_CODE = DATA[42 : 54]
                    structure.UNP_RECORDS.append(Unp_record(CHAIN_IDENT, DB_ACCESSION, ID_CODE))

            pdb_file.close()

            try:
                chains_to_remove = []
                for index in range(len(structure.CHAINS)):
                    if type(structure.CHAINS[index].RESIDUES[-1]) != Termination_symbol:
                        chains_to_remove.append(index)
                    else:
                        structure.CHAINS[index].RESIDUES.pop(-1)

                chains_to_remove.reverse()
                for index in chains_to_remove:
                    structure.CHAINS.pop(index)

            except:
                structure.data_correctnes_flag = False
                print(Colors.RED + f'The file: \'{pdb_file_path}\' contains an error in the marking of amino acid chains :O' + Colors.END)

        except:
            structure.data_correctnes_flag = False
            print(Colors.RED + f'The file \'{pdb_file_path}\' has an error in line: {line_counter} :O' + Colors.END)

    except FileNotFoundError:
        structure.data_correctnes_flag = False
        print(Colors.RED + f'File: \'{pdb_file_path}\' not found :/' + Colors.END)

    return structure

