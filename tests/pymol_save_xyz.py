from pymol import cmd

def save_xyz(filename, selection='all'):
    with open(filename, 'w') as f:
        f.write('[\n')
        atoms = cmd.get_model(selection)
        for atom in atoms.atom:
            #f.write(f"{atom.name:>4} {atom.coord[0]:>8.3f} {atom.coord[1]:>8.3f} {atom.coord[2]:>8.3f}\n")
            f.write(f"[{atom.coord[0]:>8.3f}, {atom.coord[1]:>8.3f}, {atom.coord[2]:>8.3f}],\n")
        f.write(']')

save_xyz('C:\\Users\\Hector\\OneDrive\\Non Sys\\Projects\\Python\\Projects\\rigiddock\\tests\\txt_output\\b_xyz.txt', 'all')
