import sys
import PDB
f = open(sys.argv[1])
l = f.readlines()
for line in l:
    if line.startswith("HETATM"):
        tag, atom_num, atom_name, res_name, chain_id, res_num, x, y, z, charge, radius = line.split()
        print "%-6s%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f %5.2f%6.2f" % (
                tag, int(atom_num), atom_name, ' ', res_name, ' ', int(res_num), ' ', float(x), float(y), float(z), float(charge), float(radius))
    else:
        print line.strip()
