import os, re, urllib, string, sys
import Common

rcsb_entry = "http://www.rcsb.org/pdb/cgi/export.cgi/%s.pdb?format=PDB&pdbId=%s&compression=None"
pdb_re = re.compile("(?P<prefix>pdb)*(?P<code>[0-9][0-9a-z][0-9a-z][0-9a-z])(?P<suffix>\.ent|\.pdb)*")
pdbstyle_re = re.compile("d(?P<code>[0-9][0-9a-z][0-9a-z][0-9a-z][0-9_a-z][0-9_])(?P<suffix>\.ent)*")
pdbstyle_entry = "http://astral.stanford.edu/pdbstyle.cgi?id=%s"

class PDBCache:

    def __init__(self, pdbdir, entry, mode="rw"):
        self.filename = None
        self.read = 0
        self.write = 0
        if mode.find("r") != -1:
            self.read = 1
        if mode.find("w") != -1:
            self.write = 1
            if not os.path.exists(pdbdir):
                os.makedirs(pdbdir)

        if self.write and not self.read:
            self.filename = self.getPDB(pdbdir, entry)
        elif self.read:
            ## easy out
            filename = "%s/%s" % (pdbdir, entry)
            if os.path.exists(filename):
                if Common.debug:
                    print "Existing cached file", filename
                self.filename = filename
            x = pdb_re.match(entry)
            if x:
                if Common.debug:
                    print "Appears to be a PDB entry"
                ## "canonicalize" entries, check for cached file
                if not x.groupdict()["prefix"] or not x.groupdict()["suffix"]:
                    if Common.debug:
                        print "PDB code missing prefix 'pdb' or suffix '.ent' or '.pdb':", x.groupdict()["code"]
                    entry = "pdb%s.ent" % x.groupdict()["code"]

                ## try a one-level filesystem
                filename = "%s/%s" % (pdbdir, entry)
                if Common.debug:
                    print "Looking up cached PDB", filename
                if os.path.exists(filename):
                    if Common.debug:
                        print "Existing cached PDB", filename
                    self.filename = filename
                    return

                ## try a two-level filesystem
                filename = "%s/%s/%s" % (pdbdir, x.groupdict()["code"][1:3], entry)
                if Common.debug:
                    print "Looking up cached PDB", filename
                if os.path.exists(filename):
                    if Common.debug:
                        print "Existing cached two-level PDB", filename
                    self.filename = filename
                    return

                ## full filename, retrieve from RCSB
                if self.write:
                    self.filename = self.getPDB(pdbdir, entry)

            x = pdbstyle_re.match(entry)
            if x:
                if Common.debug:
                    print "Appears to be a pdbstyle code"
                ## "canonicalize" entries
                if not x.groupdict()["suffix"]:
                    if Common.debug:
                        print "pdbstyle entry missing suffix", x.groupdict()["code"]
                    entry = "d%s.ent" % x.groupdict()["code"]

                ## try a one-level filesystem
                filename = "%s/%s" % (pdbdir, entry)
                if os.path.exists(filename):
                    if Common.debug:
                        print "Existing cached pdbstyle", filename
                    self.filename = filename
                    return

                ## try a two-level filesystem
                filename = "%s/%s/%s" % (pdbdir, x.groupdict()["code"][1:3], entry)
                if Common.debug:
                    print "Looking up cached PDB", filename
                if os.path.exists(filename):
                    if Common.debug:
                        print "Existing cached two-level PDB", filename
                    self.filename = filename
                    return

                ## full filename, retrieve from astral.stanford.edu
                if self.write:
                    self.filename = self.getPDB(pdbdir, entry)

        if not self.filename:
            if self.write:
                raise RuntimeError, "Could not locate entry %s" % entry
            elif self.read:
                raise RuntimeError, "Entry %s does not exist in cache" % entry

    def getPDB(self, pdbdir, entry):
        if Common.debug:
            print "Retrieving from the network", entry
        x = pdb_re.match(entry)
        if x:
            if Common.debug:
                print "Entry is a PDB", x.groupdict()["code"]
            upper = string.upper(x.groupdict()["code"])
            filename = "%s/pdb%s.ent" % (pdbdir, x.groupdict()["code"])
            f = open(filename, "w")
            data = urllib.urlopen(rcsb_entry % (upper, upper)).read()
            f.write(data)
            return filename

        x = pdbstyle_re.match(entry)
        if x:
            if Common.debug:
                print "Entry is a pdbstyle", x.groupdict()["code"]
            filename = "%s/d%s.ent" % (pdbdir, x.groupdict()["code"])
            request = "d%s" % x.groupdict()["code"]
            f = open(filename, "w")
            data = urllib.urlopen(pdbstyle_entry % request).read()
            f.write(data)
            return filename

        raise RuntimeError, "Unrecognized entry format", entry



def test():
    print PDBCache("/lab/db/pdb/hash", "pdb1oct.ent", "r").filename
    print PDBCache("/lab/db/pdb/hash", "1oct", "r").filename
    print PDBCache("/lab/db/pdb/hash", "pdb1oct", "r").filename
    print PDBCache("/scratch/dek/pdbstyle-1.61", "d1octc1", "rw").filename

    ## possible in cache, look it up
    print PDBCache(".", "pdb1oct.ent", "rw").filename
    os.remove("pdb1oct.ent")

    print PDBCache(".", "pdb1oct.ent", "rw").filename

    ## PDB missing pdb prefix
    print PDBCache(".", "1oct.ent", "rw").filename

    ## PDB missing pdb suffix
    print PDBCache(".", "pdb1oct", "rw").filename

    ## PDB just code
    print PDBCache(".", "1oct", "rw").filename

    ## possible in cache, look it up
    print PDBCache(".", "d1octc1.ent", "rw").filename
    os.remove("d1octc1.ent")

    print PDBCache(".", "d1octc1.ent", "rw").filename

    ## PDB missing pdb prefix
    print PDBCache(".", "d1octc1", "rw").filename

if __name__ == '__main__':
##    test()
    PDBCache(sys.argv[1], sys.argv[2], sys.argv[3])
