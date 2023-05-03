import argparse
import math
import sys

def calculate_pdb_chain_mean_minimum_distances(pdb_file_path=None):
    if pdb_file_path:
        with open(pdb_file_path, 'r') as f:
            lines = f.readlines()
    else:
        lines = sys.stdin.readlines()
    coord_dict = {}
    for line in lines:
# go only through lines where residues are
        if line.startswith('ATOM'):
    # this is the position of the chain name in the line
            chain = line[21]
            residue = line[17:20].strip()
            x_coord = float(line[30:38].strip())
            y_coord = float(line[38:46].strip())
            z_coord = float(line[46:54].strip())
            if chain not in coord_dict:
    # makes a dictionary with chain names as keys
                coord_dict[chain] = []
            if not coord_dict[chain] or residue != coord_dict[chain][-1][0]:
    # adds coordinates in a tuple as values of the dictionary for each chain
                coord_dict[chain].append([residue, [(x_coord, y_coord, z_coord)]])
            else:
                coord_dict[chain][-1][1].append((x_coord, y_coord, z_coord))

    for chain_name, res_coord in coord_dict.items():
        min_distances = []
        for i in range(len(res_coord)):
            for j in range(i+1, len(res_coord)):
                res_1 = res_coord[i][1]
                res_2 = res_coord[j][1]
                min_distance = float("inf")
                for coord1 in res_1:
                    for coord2 in res_2:
                    # assign coordinates 
                        x1, y1, z1 = coord1
                        x2, y2, z2 = coord2
                    # use the formula for the euclidian distances
                        distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
                    # search and set the minimal distance
                        if distance < min_distance:
                            min_distance = distance
                min_distances.append(min_distance)
        mean_minimal_residue = round(sum(min_distances)/len(min_distances), 4)
        print(f"{chain_name}:{mean_minimal_residue}")
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_file', nargs='?', help='PDB file path')
    args = parser.parse_args()

    calculate_pdb_chain_mean_minimum_distances(args.pdb_file)
