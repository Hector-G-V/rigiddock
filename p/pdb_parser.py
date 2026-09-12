# Parses PDB files using the mmcif package.
# For help, see https://github.com/rcsb/py-mmcif_demo/blob/master/parseSimple.ipynb
"""
In general, every function must always return predictably-structured objects. 
This is creating inter-function dependency.
"""

import os
import constants as c

import urllib.request  #for file download
from mmcif.io.IoAdapterCore import IoAdapterCore


def parseURL():
    url = c.URL 
    
    # ... put here the code to download a cif file from the PDB,
    # then parse the file. Code for this is in the GitHub for 
    # the package mmcif, under demo.


def parseLocal(file_path):
    
    # filename: the file from which to extract data
    # Returns the coordinates of atoms in the protein, extracted from the PDB file
    # Work on building the comments later


    # Retrieve the atom_site data from the pdb file

    io = IoAdapterCore()
    list_data_container = io.readFile(file_path)
    # read file, generate list of data containers
    
    data_container = list_data_container[0]
    # select the 1st data container
    # I'm not sure there are ever more than one container

    atom_site = data_container.getObj('atom_site')
    # get all data in the entity_poly_seq category


    # ----------------------------------------------------------------------
    # Get column indices
    # ----------------------------------------------------------------------

    label_comp_id_idx = atom_site.getAttributeIndex("label_comp_id")
    label_alt_id_idx = atom_site.getAttributeIndex("label_alt_id")
    occupancy_idx = atom_site.getAttributeIndex("occupancy")

    label_asym_id_idx = atom_site.getAttributeIndex("label_asym_id")
    label_seq_id_idx = atom_site.getAttributeIndex("label_seq_id")
    label_atom_id_idx = atom_site.getAttributeIndex("label_atom_id")

    # Optional but recommended for unusual residue numbering
    ins_code_idx = atom_site.getAttributeIndex("pdbx_PDB_ins_code")

    # ----------------------------------------------------------------------
    # Pass 1: remove waters
    # ----------------------------------------------------------------------

    rows = [
        row for row in atom_site.data
        if row[label_comp_id_idx] != "HOH"
    ] # Filter out rows where _atom_site.label_comp_id is 'HOH'

    # ----------------------------------------------------------------------
    # Pass 2: Determine the highest-occupancy alternate conformer
    # ----------------------------------------------------------------------

    """
    Algorithm for Pass 2 & 3:
        For atoms without an alternate location 
        (label_alt_id is . or ? or empty), keep them unchanged.

        For atoms with alternate locations (A, B, etc.), 
        compare the rows representing the same atom.
        
        Keep the row with the highest occupancy.

        If occupancies are tied, keep the first one 
        encountered (Python's > comparison naturally does this).

    The algo creates gaps in _atom_site.id. This change is valid in mmcif format. 
    Recall that _atom_site.id is just a unique identifier, not a sorted sequence.
    """
    best_alt = {}

    for row in rows:

        alt_id = row[label_alt_id_idx]

        # No alternate conformation
        if alt_id in (".", "?", "", None):
            continue

        key = (
            row[label_asym_id_idx],
            row[label_seq_id_idx],
            row[ins_code_idx],
            row[label_comp_id_idx],
            row[label_atom_id_idx],
        )

        occupancy = float(row[occupancy_idx])

        # Keep the first row encountered if occupancies are equal
        if key not in best_alt or occupancy > best_alt[key][0]:
            best_alt[key] = (occupancy, row)

    # ----------------------------------------------------------------------
    # Pass 3: Filter the rows while preserving their original order
    # ----------------------------------------------------------------------

    filtered_rows = []

    for row in rows:

        alt_id = row[label_alt_id_idx]

        # No alternate conformation
        if alt_id in (".", "?", "", None):
            filtered_rows.append(row)
            continue

        key = (
            row[label_asym_id_idx],
            row[label_seq_id_idx],
            row[ins_code_idx],
            row[label_comp_id_idx],
            row[label_atom_id_idx],
        )

        if row is best_alt[key][1]:
            filtered_rows.append(row)

    # Replace the atom_site table
    atom_site.setRowList(filtered_rows)


    # ----------------------------------------------------------------------
    # Isolate the Cartesian coordinates
    # ----------------------------------------------------------------------

    x = atom_site.getAttributeValueList('Cartn_x')
    y = atom_site.getAttributeValueList('Cartn_y')
    z = atom_site.getAttributeValueList('Cartn_z')
    # get list of values by attr Cartn_x, Cartn_y, Cartn_z

    x = list(map(float, x))
    y = list(map(float, y))
    z = list(map(float, z))
    # convert from string to float

    data = []
    for i in range(len(x)):
        data.append([x[i], y[i], z[i]])
    # Organize the data into points (x,y,z) for output


    #return data
    return data