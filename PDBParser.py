import PDB, Common, string, sys, os

def addResidueToResidueList(curr_residue_list, curr_res_num, curr_res_name, curr_atom_list):
    if Common.debug:
        print "Adding atoms to residue '%s' %s: %s" % ( \
            curr_res_name,
            curr_res_num,
            map(lambda x: x.getName(), curr_atom_list))
    r = PDB.Residue(number = curr_res_num,
                name = curr_res_name,
                atoms = curr_atom_list,
                chain = None)

    ## store the chain to the residue accumulation list
    curr_residue_list.append(r)
##     if Common.debug: print "Current residue list is", map(lambda x: (x.getName(), x.getNumber()), curr_residue_list)


    ## store the back-pointer for each atom to its residue
    for atom in curr_atom_list:
        atom.setResidue(r)


def addChainToChainList(chains, curr_chain, curr_residue_list):
    chains[curr_chain] = PDB.Chain(curr_chain, curr_residue_list)

    ## save backpointers from each residue to its chain
    for residue in curr_residue_list:
        residue.setChain(chains[curr_chain])



class PDBParser:
    def __init__(self, handle):
        self._handle = handle

    def parse(self):
        extra_records = []
        ## data structures that are saved after parsing
        models = {}
        seqres_chains = {}

        ## running data structures used to accumulate per-chain, per-residue, and per-atom data
        curr_seqres_chain = None

        chains = {}
        curr_model_number = None
        curr_res_num = None
        curr_res_name = None
        curr_atom_list = []
        curr_residue_list = []
        curr_chain = None
        chain_ter_seen = {}

        while 1:
            line = self._handle.readline()
            if Common.debug: print "PDBParser line: %s" % line,

            ## end of file
            if line == '':
                if Common.debug:
                    print "PDBParser end of file"
                    print "curr_atom_list is", curr_atom_list
                    print "curr_residue_list is", curr_residue_list
                    print "curr_chain is", curr_chain
                if curr_atom_list != [] and curr_res_name not in Common.residue_skip_list:
                    addResidueToResidueList(curr_residue_list, curr_res_num, curr_res_name, curr_atom_list)
                if curr_residue_list != [] and curr_res_name not in Common.residue_skip_list:
                    addChainToChainList(chains, curr_chain, curr_residue_list)
                ## in a single-model file without MODEL and ENDMDL records, we must set this:
                if curr_model_number == None:
                    if Common.debug: print "PDBParser end of file, autosetting model 1"
                    models[1] = PDB.Model(chains)
                break

            ## if we see a new model record, record its number
            elif line[:5] == 'MODEL':
                ## the official format:
##                 curr_model_number = int(line[10:14])
                ## fix/hack due to TAB stored in pdbstyle files:
                curr_model_number = int(string.split(line)[1])
                if Common.debug: print "PDBParse: model entry", curr_model_number
                continue

            ## at the end of a new model, clear away state data.
            ## The TER record before ENDMDL would have dealt with adding the model to the chain

            elif line[:6] == 'ENDMDL':
                if Common.debug: print "PDBParse end of model entry", curr_model_number
                models[curr_model_number] = PDB.Model(chains)
                ## reset the chain_ter_seen for the new model
                chains = {}
                curr_res_num = None
                curr_res_name = None
                curr_atom_list = []
                curr_residue_list = []
                curr_chain = None
                chain_ter_seen = {}
                continue


            elif line[:6] == 'SEQRES':
                chain_id = line[11]
                chain_data = string.split(line[19:70])

                if Common.debug: print "PDBParser SEQRES:", chain_id, chain_data

                if Common.debug: print "SEQRES", chain_id, chain_data
                if curr_seqres_chain == None or curr_seqres_chain != chain_id:
                    if Common.debug: print "Setup new chain", chain_id
                    curr_seqres_chain = chain_id
                    seqres_chains[curr_seqres_chain] = ""
                if Common.debug: print "Store SEQRES chain '%s' (%s)" % (chain_id, chain_data)
                for res in chain_data:
                    seqres_chains[curr_seqres_chain] += Common.three_to_one(res)

            elif line[:4] == 'ATOM' or line[:6] == 'HETATM':
                atom_num = int(line[6:11])
                atom_name = string.strip(line[12:16])
                alternate = line[16]
                res_name = string.strip(line[17:20])
                chain_id = line[21]
                res_num = string.strip(line[22:27])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                try:
                    occupancy = float(line[55:60])
                except:
                    occupancy = 0.0
                try:
                    bfactor = float(line[60:66])
                except:
                    bfactor = 0.0

                try:
                    element = line[77]
                except IndexError:
                    element = ' '

##                 if Common.debug: print "PDBParser ATOM:", atom_num, atom_name, res_name, chain_id, res_num, x, y, z, occupancy, bfactor, element
                ## if we've seen a TER record for this chain already, ignore any further entries with that
                ## chain_id, because these are just associated ions which will not contribute to the sequence
                ## this is broken if people re-use chain IDs, but that would be... stupid!
                if chain_ter_seen.has_key(chain_id):
                    continue

                ## new chain; store previously accumulated residues
                if curr_chain != chain_id and curr_res_name not in Common.residue_skip_list:
                    if Common.debug: print "New chain", chain_id, "curr chain", curr_chain
                    if curr_chain != None:
                        ## build a new chain object
                        ## store all the accumulated residues into this chain
                        if curr_atom_list != None and curr_res_num != None and curr_res_name != None:
                            addResidueToResidueList(curr_residue_list, curr_res_num, curr_res_name, curr_atom_list)
                            if Common.debug: print "Adding residues to chain %s: %s" % (curr_chain, map(lambda x: (x.getName(), x.getNumber()), curr_residue_list))

                            addChainToChainList(chains, curr_chain, curr_residue_list)
                    ## store the new current chain
                    curr_chain = chain_id
                    ## clear the accumulated residue
                    curr_residue_list = []
                    ## clear the current residue number
                    curr_res_num = None
                    curr_res_name = None

                ## new residue; store the previously accumulated atoms
                if curr_res_num != res_num:

                    ## finished parsing a residue, so build a new
                    ## residue object and store all the accumulated
                    ## atoms into this residue

                    if curr_res_num != None and curr_res_name not in Common.residue_skip_list:
                        addResidueToResidueList(curr_residue_list, curr_res_num, curr_res_name, curr_atom_list)

                    if Common.debug: print "New residue", res_num

                    ## rememebr the new current residue
                    curr_res_num = res_num
                    curr_res_name = res_name

                    ## clear the accumulated atom list
                    curr_atom_list = []

                ## build a new atom object
                if Common.debug: print "New atom: '%s', '%s', '%s', '%s', '%s'" % (atom_num, atom_name, alternate, res_name, res_num)
                isHetatm = 0
                if line[:6] == 'HETATM':
                    isHetatm = 1
                atom = PDB.Atom(serial = atom_num,
                            name = atom_name,
                            alternate = alternate,
                            residue = None,
                            coords = (x,y,z),
                            occupancy = occupancy,
                            bfactor = bfactor,
                            element = element,
                            isHetatm = isHetatm)
                curr_atom_list.append(atom)

            ## when we see a TER record, save the chain it came from
            ## from because some PDB files put ions 'associated' with the chain
            ## after the TER record:
            ## ATOM   3660  OXT GLN B 502     129.488  87.534  67.598  1.00168.80           O
            ## TER    3661      GLN B 502
            ## HETATM 3662  S   SO4 B   1      99.307  73.882  58.307  1.00 63.92           S

            elif line[:3] == 'TER':
                if Common.debug: print "PDBParser: TER on chain", chain_id
                try:
                    chain_id = line[21]
                    chain_ter_seen[chain_id] = 1
                except:
                    chain_ter_seen[curr_chain] = 1
                addResidueToResidueList(curr_residue_list, curr_res_num, curr_res_name, curr_atom_list)
                addChainToChainList(chains, curr_chain, curr_residue_list)
                curr_res_num = None
                curr_res_name = None
                curr_atom_list = []
                curr_residue_list = []

            elif line[:4] == 'END ' or line[:6] == 'CONECT':
                pass
            else:
                extra_records.append(line)

        ## Build the sequences after everything else is done
        ## (I did this in case the SEQRES came after ATOM records, even though that's illegal)
        for chain_id in chains.keys():
            ## store the per-chain seqres data in the chain
            if seqres_chains.has_key(chain_id):
                if Common.debug: print "Storing SEQRES chain '%s' (%s)" % (chain_id,  seqres_chains[chain_id])
                chains[chain_id].setSeqres(seqres_chains[chain_id])

        return PDB.PDB(models, extra_records)

def main():
##     Common.debug = 1
    import getopt
    pdb = None
    chain = None
    seqres = 0
    outputPrefix = "PDB"
    output_fasta = 0
    output_pdb = 0

    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:c:o:sfd",["pdb=", "chain=","outputPrefix=","seqres","fasta","pdb"])
    except getopt.GetoptError, what:
	raise RuntimeError, "usage"
    for o, a in opts:
        if o in ('-p', 'pdb'):
            pdbFile = a
        if o in ('-c', 'chain'):
            chain = a
        if o in ('-o', 'outputPrefix'):
            outputPrefix = a
        if o in ('-s', 'seqres'):
            seqres = 1
        if o in ('-f', 'fasta'):
            output_fasta = 1
        if o in ('-d', 'pdb'):
            output_pdb = 1

    if not pdbFile or not chain:
        raise RuntimeError, "missing pdbfile or chain"

    f = open(pdbFile)
    p = PDBParser(f).parse()
    ms = p.getModels().keys()
    ms.sort()
    m = p.getModel(ms[0])

    if output_pdb:
        pdb_output = open("%s.pdb" % outputPrefix, "w")
        if chain == "-":
            p.formatPDB(pdb_output, None, not seqres)
        else:
            p.formatPDB(pdb_output, [chain], not seqres)

    if output_fasta:
        fasta_output = open("%s.fa" % outputPrefix, "w")
        if chain == '-':
            fasta_output.write(">%s all chains\n" % (os.path.basename(pdbFile)))
            cs = m.getChains().keys()
            for chain in cs:
                c = m.getChain(chain)
                if seqres:
                    seq = c.getSeqres()
                else:
                    seq = c.getSequence()
                fasta_output.write("%s/" % seq)
            fasta_output.write("*\n")
        else:
            fasta_output.write(">%s chain %s\n" % (os.path.basename(pdbFile), chain))
            c = m.getChain(chain)
            if seqres:
                fasta_output.write("%s\n" % c.getSeqres())
            else:
                fasta_output.write("%s\n" % c.getSequence())

def test():
    p = PDBParser(open("/lab/db/pdb/hash/oc/pdb1oca.ent")).parse()
    p.formatPDB(sys.stdout, None, 1)
if __name__ == '__main__':
    main()
