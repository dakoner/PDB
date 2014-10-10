debug = 0


residue_skip_list = ["HOH", "WAT"]
nucleic = [ 'A', 'G', 'T', 'C', 'U']

nucleic_three_to_one = {
    'ADE': 'A',
    'GUA': 'G',
    'THY': 'T',
    'CYT': 'C',
    'URA': 'U',
    }

nucleic_one_to_three = {
    'A': 'ADE',
    'G': 'GUA',
    'T': 'THY',
    'C': 'CYT',
    'U': 'URA'
    }

ligand_resnames = [
    'HOH', 'WAT'
    ]

protein_one_to_three = {
    'A': 'ALA',
    'C': 'CYS',
    'D': 'ASP',
    'E': 'GLU',
    'F': 'PHE',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'K': 'LYS',
    'L': 'LEU',
    'M': 'MET',
    'N': 'ASN', 
    'P': 'PRO',
    'Q': 'GLN',
    'R': 'ARG',
    'S': 'SER',
    'T': 'THR',
    'V': 'VAL',
    'W': 'TRP',
    'Y': 'TYR',
    'X': 'UNK'
    }

# This table is taken from the RAF release notes, and includes the
# undocumented mapping "UNK" -> "X"
protein_three_to_one = {
    'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
    'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
    'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
    'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
    '2AS':'D', '3AH':'H', '5HP':'E', 'ACL':'R', 'AIB':'A',
    'ALM':'A', 'ALO':'T', 'ALY':'K', 'ARM':'R', 'ASA':'D',
    'ASB':'D', 'ASK':'D', 'ASL':'D', 'ASQ':'D', 'AYA':'A',
    'BCS':'C', 'BHD':'D', 'BMT':'T', 'BNN':'A', 'BUC':'C',
    'BUG':'L', 'C5C':'C', 'C6C':'C', 'CCS':'C', 'CEA':'C',
    'CHG':'A', 'CLE':'L', 'CME':'C', 'CSD':'A', 'CSO':'C',
    'CSP':'C', 'CSS':'C', 'CSW':'C', 'CXM':'M', 'CY1':'C',
    'CY3':'C', 'CYG':'C', 'CYM':'C', 'CYQ':'C', 'DAH':'F',
    'DAL':'A', 'DAR':'R', 'DAS':'D', 'DCY':'C', 'DGL':'E',
    'DGN':'Q', 'DHA':'A', 'DHI':'H', 'DIL':'I', 'DIV':'V',
    'DLE':'L', 'DLY':'K', 'DNP':'A', 'DPN':'F', 'DPR':'P',
    'DSN':'S', 'DSP':'D', 'DTH':'T', 'DTR':'W', 'DTY':'Y',
    'DVA':'V', 'EFC':'C', 'FLA':'A', 'FME':'M', 'GGL':'E',
    'GLZ':'G', 'GMA':'E', 'GSC':'G', 'HAC':'A', 'HAR':'R',
    'HIC':'H', 'HIP':'H', 'HMR':'R', 'HPQ':'F', 'HTR':'W',
    'HYP':'P', 'IIL':'I', 'IYR':'Y', 'KCX':'K', 'LLP':'K',
    'LLY':'K', 'LTR':'W', 'LYM':'K', 'LYZ':'K', 'MAA':'A',
    'MEN':'N', 'MHS':'H', 'MIS':'S', 'MLE':'L', 'MPQ':'G',
    'MSA':'G', 'MSE':'M', 'MVA':'V', 'NEM':'H', 'NEP':'H',
    'NLE':'L', 'NLN':'L', 'NLP':'L', 'NMC':'G', 'OAS':'S',
    'OCS':'C', 'OMT':'M', 'PAQ':'Y', 'PCA':'E', 'PEC':'C',
    'PHI':'F', 'PHL':'F', 'PR3':'C', 'PRR':'A', 'PTR':'Y',
    'SAC':'S', 'SAR':'G', 'SCH':'C', 'SCS':'C', 'SCY':'C',
    'SEL':'S', 'SEP':'S', 'SET':'S', 'SHC':'C', 'SHR':'K',
    'SOC':'C', 'STY':'Y', 'SVA':'S', 'TIH':'A', 'TPL':'W',
    'TPO':'T', 'TPQ':'A', 'TRG':'K', 'TRO':'W', 'TYB':'Y',
    'TYQ':'Y', 'TYS':'Y', 'TYY':'Y', 'AGM':'R', 'GL3':'G',
    'SMC':'C', 'ASX':'B', 'CGU':'E', 'CSX':'C', 'GLX':'Z',
    'UNK':'X',

## listed as automatically or manually curated on ASTRAL site
    "143": "C",
    "2MR": "R",
    "AAR": "R",
##     "ABA": "C", ## ambiguity
##     "ABA": "N", ## ambiguity
    "ACA": "A",
    "ACY": "G",
    "AEI": "T",
    "ARO": "R",
    "ASI": "N",
    "BFD": "D",
    "BTA": "E",
    "BTC": "G",
    "CAF": "C",
    "CAS": "C",
    "CAY": "C",
    "CCY": "S",
    "CLB": "S",
    "CLD": "S",
    "CMT": "C",
    "CRG": "S",
##     "CRO": "G",## weird cyclization used for different residues
##     "CRO": "S",## weird cyclization used for different residues
    "CRQ": "Q",
    "CSB": "C",
    "CSE": "C",
    "CSH": "S",
    "CSR": "C",
    "CSY": "S",
    "CSZ": "C",
    "CZZ": "C",
    "DOH": "D",
    "EHP": "F",
    "ESD": "S",
    "FGL": "C",
    "FTR": "W",
    "GLH": "E",
    "H2P": "H",
    "HPH": "F",
##     "IAS": "D",## ambiguity betwen 1at6 and other PDB files
##     "IAS": "N",## ambiguity betwen 1at6 and other PDB files
    "IIC": "S",
    "INI": "K",
    "LYN": "K",
    "M3L": "K",
    "MGN": "Q",
    "MHO": "M",
    "MLY": "K",
    "MLZ": "K",
    "MME": "M",
    "MSO": "M",
    "MTY": "F",
    "NIY": "Y",
    "NPH": "C",
    "OCY": "C",
    "PBI": "S",
    "PIA": "A",
    "PRS": "P",
    "PVL": "S",
    "PYX": "C",
    "S1H": "S",
    "SBD": "S",
    "SBL": "S",
    "SEB": "S",
    "SEG": "S",
    "SME": "M",
    "SNC": "C",
    "SNN": "D",
    "SUI": "D",
    "TRF": "W",
    "TRN": "W",
    "TYI": "Y",
    "TYN": "Y",
    "YOF": "Y",
}

def three_to_one(res):
    if protein_three_to_one.has_key(res):
        return protein_three_to_one[res]
    elif res in ligand_resnames:
        return "?"
    elif nucleic_three_to_one.has_key(res):
        return nucleic_three_to_one[res]
    elif res in nucleic:
        return res
    else:
##        raise RuntimeError, "PDB file has unknown residue '%s'" % res
        return 'X'
