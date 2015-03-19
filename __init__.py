## a mini PDB parser which provides all the functionality we need
import string, sys, os
import Common


## generic PDBError exception, never used directly
class PDBError(Exception):
    pass

class ResidueNumberError(PDBError):
    def __init__(self, number):
        self.number = number
    def __str__(self):
        return "Non-existent residue number %d\n" % number

class ResidueIntegrityError(PDBError):
    def __init__(self, number):
        self.number = number
    def __str__(self):
        return "Integrity constraint failed: more than one residue numbered %d in chain\n" % number

    
class Atom:
    def __init__(self, serial=1, name=None, alternate=" ", residue=None, coords=(0.,0.,0.), occupancy=1.0, bfactor=1.0, element=' ', isHetatm=0):
        self._serial = serial
        self._name = name
        self._alternate = alternate
        self._residue = residue
        self._coords = coords
        self._occupancy = occupancy
        self._bfactor = bfactor
        self._element = element
        self._isHetatm = isHetatm

    def getSerial(self):
        return self._serial

    def getName(self):
        return self._name

    ## note: changing the residue should change the residue's atom ownership
    def setResidue(self, residue):
        self._residue = residue

    def getResidue(self):
        return self._residue

    def setCoords(self, coords):
        self._coords = coords

    def getCoords(self):
        return self._coords
    
    def setOccupancy(self, occupancy):
        self._occupancy = occupancy

    def getOccupancy(self):
        return self._occupancy

    def setBfactor(self, bfactor):
        self._bfactor = bfactor

    def getBfactor(self):
        return self._bfactor

    def getRecord(self):
        record_type = 'ATOM'
        if self._isHetatm:
            record_type = 'HETATM'
        residue_name = self._residue.getName()
        residue_number = self._residue.getNumber()
        if len(self._name) == 4 or (self._name[0] in string.digits and self._name[1] == 'H'):
            atom_name = "%-4s" % self._name
        else:
            atom_name = " %-3s" % self._name
        if residue_number[-1] in string.ascii_letters:
            resnum = residue_number[:-1]
            insertion_code = residue_number[-1]
        else:
            resnum = residue_number
            insertion_code = " "
        if len(residue_name) == 4:
            code = "%-6s%5d %4s%c%4s%c%4s%c   %8.3f%8.3f%8.3f %5.2f%6.2f" % (record_type, self._serial, atom_name, insertion_code, self._residue.getName(),  self._residue.getChain().getName(), resnum, insertion_code, self._coords[0], self._coords[1], self._coords[2], self._occupancy, self._bfactor)
        else:
            code = "%-6s%5d %4s%c%3s %c%4s%c   %8.3f%8.3f%8.3f %5.2f%6.2f" % (record_type, self._serial, atom_name, self._alternate, self._residue.getName(),  self._residue.getChain().getName(), resnum, insertion_code, self._coords[0], self._coords[1], self._coords[2], self._occupancy, self._bfactor)
            
        return code
             
    def __str__(self):
        return  "Atom %d, %s @ (%5.2f %5.2f %5.2f)" % (self._serial, self._name, self._coords[0], self._coords[1], self._coords[2])

        
class Residue:
    def __init__(self, number=1, name=None, chain=None, atoms=None):
        self._number = number
        self._name = name
        self._chain = chain
        self._atoms = atoms

    def getNumber(self):
        return self._number

    def getName(self):
        return self._name

    def getAtoms(self):
        return self._atoms

    def getAtom(self, name):
        for atom in self._atoms:
            if atom.getName() == name:
                return atom
        return None

    def setChain(self, chain):
        self._chain = chain

    def getChain(self):
        return self._chain

    def __str__(self):
        return "Residue '%s', '%s'" % (self._number, self._name)

class Chain:
    def __init__(self, name, residues=None, seqres=None, model=None):
        self._name = name
        self._residues = residues
        self._sequence = None
        self._seqres = None
        self._model = None

        self.buildSequence()
        
    def buildSequence(self):
        ## build the chain's sequence from the residue list
        sequence = [" "] * len(self._residues)
        for i in range(len(self._residues)):
            residue = self._residues[i]
            try:
                res = Common.three_to_one(residue.getName())
            except RuntimeError:
                res = "?"
            sequence[i] = res
        self._sequence = string.join(sequence, '')

    def setModel(self, model):
        self._model = model

    def getModel(self):
        return _model
    
    def getName(self):
        return self._name

    def getResidues(self):
        return self._residues
    
        
    
    def getResidue(self, number):
        x = filter(lambda x: x.getNumber() == number, self._residues)
        if len(x) == 0:
            msg = "Residue numbered %d not in chain" % number
            raise ResidueNumberError, number
        elif len(x) > 1:
            raise ResidueIntegrityConstraint, number
        return x[0]
        
    def getSequence(self):
        return self._sequence

    def setSeqres(self, seqres):
        self._seqres = seqres

    def getSeqres(self):
        return self._seqres

    def __str__(self):
        return "Chain %s\nSequence: %s\nSEQRES Sequence: %s" % (self.getName(), self.getSequence(), self.getSeqres())

class Model:
    def __init__(self, chains=None):
        self._chains = chains

    def getChain(self, chain):
        return self._chains[chain]
    
    def getChains(self):
        return self._chains
    
class PDB:
    def __init__(self, models, extra_records):
        self._models = models
        self._extra_records = extra_records

    def getModels(self):
        return self._models

    def getModel(self, model_number):
        return self._models[model_number]

    def __str__(self):
        s = ""
        for key in self._models.keys():
            model = self._models[key]
            chains = model.getChains()
            chainkeys = chains.keys()
            chainkeys.sort()
            for chainkey in chainkeys:
                chain = chains[chainkey]
                s += chain.__str__() + "\n"
                for residue in chain.getResidues():
                    s += residue.__str__() + "\n"
                    for atom in residue.getAtoms():
                        s += atom.__str__() + "\n"
        return s


    ## destructively transform the second pdb onto the first using the rotation matrix from the structure match
    def transformPDB(self,matrix):
        ## for every model,
        for key in self._models.keys():
            model = self._models[key]
            chains = model.getChains()
            chainkeys = chains.keys()
            chainkeys.sort()
            ## for every chain, 
            for chainkey in chainkeys:
                chain = chains[chainkey]
                residues = chain.getResidues()
                ## for every residue, 
                for res in residues:
                    atoms = res.getAtoms()
                    ## for every atom,
                    for atom in atoms:
                        c = atom.getCoords()
                        ## transform the coordinates 
                        newc = Numeric.innerproduct(Numeric.array((c[0], c[1], c[2], 1.)), matrix)
                        atom.setCoords((newc[0], newc[1], newc[2]))

    ## write a PDB, including SEQRES records, to a file
    def formatPDB(self, handle, chainkeys = None, use_atom_seqres=1):
        for record in self._extra_records:
            handle.write(record)
        modelkeys = self._models.keys()
        ## print SEQRES only once, for model #1
        model = self._models[modelkeys[0]]
        chains = model.getChains()
        if chainkeys == None:
            chainkeys = chains.keys()
            chainkeys.sort()
            
        for chainkey in chainkeys:
            chain = chains[chainkey]
            if use_atom_seqres:
                residues = chain.getResidues()
                for i in range(0, len(residues), 13):
                    handle.write("SEQRES%4d %c%5d " % ((i/13)+1, chain.getName(), len(residues)))
                    for j in range(i, i+13):
                        if j == len(residues):
                            break
                        handle.write(" %3s" % residues[j].getName())
                    handle.write("\n")
            else:
                seq = chain.getSeqres()
                for i in range(0, len(seq), 13):
                    handle.write("SEQRES%4d %c%5d " % ((i/13)+1, chain.getName(), len(seq)))
                    for j in range(i, i+13):
                        if j == len(seq):
                            break
                        handle.write(" %3s" % protein_one_to_three[seq[j]])
                    handle.write("\n")
 
        for key in modelkeys:
            model = self._models[key]
            if len(modelkeys) > 1:
                handle.write("MODEL %d\n" % key)
            for chainkey in chainkeys:
                chain = chains[chainkey]
                residues = chain.getResidues()
                for res in residues:
                    atoms = res.getAtoms()
                    for atom in atoms:
                        handle.write(atom.getRecord())
                        handle.write("\n")
                handle.write("TER\n")
            if len(modelkeys) > 1:
                handle.write("ENDMDL\n")
        handle.write("END\n")
        handle.close()
