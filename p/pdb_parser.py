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


    # Remove the waters

    label_comp_id_index = atom_site.getAttributeIndex('label_comp_id')
    # Get index of _atom_site.label_comp_id column. 3-letter AA code, HOH, etc
    
    filtered_rows = [row for row in atom_site.data if row[label_comp_id_index] != 'HOH']
    # Filter out rows where _atom_site.label_comp_id is 'HOH'
    
    atom_site.setRowList(filtered_rows)
    # Overwrite the atom_site data with the filtered rows
    

    # Isolate the Cartesian coordinates

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


    return data
    # Change to this very soon in the future. For now, keep the same
