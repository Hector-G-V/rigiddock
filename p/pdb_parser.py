# Parses PDB files using the mmcif package.
# For help, see https://github.com/rcsb/py-mmcif_demo/blob/master/parseSimple.ipynb
"""
In general, every function must always return predictably-structured objects. 
This is creating inter-function dependency.
"""

import os
import string
import constants as c
import urllib.request  #for file download

import numpy as np

from mmcif.io.IoAdapterCore import IoAdapterCore
from mmcif.api.PdbxContainers import DataContainer
from mmcif.api.DataCategory import DataCategory


def parseURL():
    url = c.URL 
    
    # ... put here the code to download a cif file from the PDB,
    # then parse the file. Code for this is in the GitHub for 
    # the package mmcif, under demo.


def preprocessing(file_path):

    # Execute the preprocessing options given by the user. At the moment:
    #   Remove water
    #   Remove alternate locations
    
    # file_path: the file from which to extract data
    # Returns the data container
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

    return atom_site


def atom_site_xyz(atom_site):

    # Isolate the Cartesian coordinates from the data container
    # filename: the file from which to extract data
    # Returns the coordinates of atoms in the protein, extracted from the PDB file
    # Work on building the comments later

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


def parseLocal(file_path):
    
    # Performs the preprocessing and xyz coordinate retrieval from the pdb file
    # filename: the file from which to extract data
    # Returns the coordinates of atoms in the protein, extracted from the PDB file
    # Work on building the comments later

    # Remove waters and alternate sites (for now)
    preprocess = preprocessing(file_path)

    # Extract xyz coordinates
    coordinates = atom_site_xyz(preprocess)

    return coordinates


def write_docked_file(
    filepath_A,
    filepath_B,
    transformed_coordinates_B,
    output_file,
):
    """
    Create a new mmCIF containing proteins A and B.

    Each input may contain multiple chains.

    Protein A:
        - Coordinates unchanged.
        - All chains retained.
        - Assigned unique output chain IDs.
        - Assigned unique entity IDs.

    Protein B:
        - Coordinates replaced by transformed_coordinates_B.
        - All chains retained.
        - Assigned unique output chain IDs that do not overlap A's.
        - Assigned unique entity IDs.

    Atom IDs are renumbered sequentially.

    The output uses a common _atom_site schema so that 
    A and B do not need to have identical mmCIF schemas.

    Parameters
    ----------
    filepath_A : str
        Filepath to protein A.

    filepath_B : str
        Filepath to protein B.

    transformed_coordinates_B : numpy.ndarray
        Array of shape (N, 3) containing the transformed
        coordinates for all retained atoms in B, in the same
        order as the rows returned by preprocessing().

    output_file : str
        Path of the output mmCIF file.
    """

    # ------------------------------------------------------------------
    # 1. Preprocess the input structures.
    # ------------------------------------------------------------------

    # Preprocessing: remove HOH, alternate sites, etc.
    # Returns preprocessed mmCIF container for proteins
    atom_site_A = preprocessing(filepath_A)
    atom_site_B = preprocessing(filepath_B)

    # ------------------------------------------------------------------
    # 2. Validate transformed coordinates
    # ------------------------------------------------------------------

    transformed_coordinates_B = np.asarray(
        transformed_coordinates_B
    )

    if transformed_coordinates_B.shape != (
        atom_site_B.getRowCount(),
        3,
    ):
        raise ValueError(
            "transformed_coordinates_B must have shape "
            f"({atom_site_B.getRowCount()}, 3)"
        )

    if not np.isfinite(transformed_coordinates_B).all():
        raise ValueError(
            "transformed_coordinates_B contains non-finite coordinates."
        )

    # ------------------------------------------------------------------
    # 3. Validate required input columns
    # ------------------------------------------------------------------

    required_attributes = [
        "type_symbol",
        "label_atom_id",
        "label_comp_id",
        "label_asym_id",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
    ]

    for name, atom_site in [
        ("A", atom_site_A),
        ("B", atom_site_B),
    ]:
        missing = [
            attr for attr in required_attributes
            if not atom_site.hasAttribute(attr)
        ]

        if missing:
            raise ValueError(
                f"Input {name} is missing required _atom_site "
                f"attributes: {missing}"
            )

    # ------------------------------------------------------------------
    # 4. Collect original chain IDs in their order of appearance
    # ------------------------------------------------------------------

    def get_chain_ids(atom_site):
        chain_idx = atom_site.getAttributeIndex("label_asym_id")

        # dict preserves insertion order in modern Python.
        return list(dict.fromkeys(
            row[chain_idx] for row in atom_site.data
        ))

    chains_A = get_chain_ids(atom_site_A)
    chains_B = get_chain_ids(atom_site_B)

    if not chains_A or not chains_B:
        raise ValueError(
            "Both docking partners must contain at least one atom "
            "and one chain."
        )

    # ------------------------------------------------------------------
    # 5. Assign unique output chain IDs
    #
    # Use A-Z, a-z, 0-9, then extend to multi-character IDs if needed.
    # Chain IDs are identifiers; their names do not affect coordinates.
    # ------------------------------------------------------------------

    def generate_chain_ids(n):
        alphabet = string.ascii_uppercase + string.ascii_lowercase + string.digits
        base = len(alphabet)

        ids = []

        for i in range(n):
            if i < base:
                ids.append(alphabet[i])
            else:
                # Extend to multi-character IDs.
                value = i
                chars = []

                while value >= base:
                    chars.append(alphabet[value % base])
                    value = value // base - 1

                chars.append(alphabet[value])
                ids.append("".join(reversed(chars)))

        return ids

    all_chains = chains_A + chains_B
    output_chain_ids = generate_chain_ids(len(all_chains))

    chain_map_A = dict(zip(chains_A, output_chain_ids[:len(chains_A)]))
    chain_map_B = dict(zip(chains_B, output_chain_ids[len(chains_A):]))

    # Assign one entity ID per output chain.
    # This deliberately keeps the MVP mapping simple.
    entity_map_A = {
        chain: str(i + 1)
        for i, chain in enumerate(chains_A)
    }

    entity_map_B = {
        chain: str(len(chains_A) + i + 1)
        for i, chain in enumerate(chains_B)
    }

    # ------------------------------------------------------------------
    # 6. Source IDs for provenance
    # ------------------------------------------------------------------

    source_A = os.path.basename(filepath_A)
    source_B = os.path.basename(filepath_B)

    # ------------------------------------------------------------------
    # 7. Create output DataContainer and title
    # ------------------------------------------------------------------

    output_container = DataContainer("rdock_docked")

    struct = DataCategory(
        "struct",
        ["entry_id", "title"],
    )

    struct.append([
        "rdock_docked",
        (
            "Protein-protein rigid-body docking model generated "
            f"by RDock from {source_A} and {source_B}"
        ),
    ])

    output_container.append(struct)

    # ------------------------------------------------------------------
    # 8. Create entity category
    # ------------------------------------------------------------------

    entity = DataCategory(
        "entity",
        ["id", "type", "pdbx_description"],
    )

    for chain in chains_A:
        entity.append([
            entity_map_A[chain],
            "polymer",
            f"Protein A, source chain {chain}",
        ])

    for chain in chains_B:
        entity.append([
            entity_map_B[chain],
            "polymer",
            f"Protein B, source chain {chain}",
        ])

    output_container.append(entity)

    # ------------------------------------------------------------------
    # 9. Create struct_asym category
    # ------------------------------------------------------------------

    struct_asym = DataCategory(
        "struct_asym",
        ["id", "entity_id"],
    )

    for chain in chains_A:
        struct_asym.append([
            chain_map_A[chain],
            entity_map_A[chain],
        ])

    for chain in chains_B:
        struct_asym.append([
            chain_map_B[chain],
            entity_map_B[chain],
        ])

    output_container.append(struct_asym)

    # ------------------------------------------------------------------
    # 10. Define a common _atom_site schema
    # ------------------------------------------------------------------

    atom_site_attributes = [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "label_alt_id",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
    ]

    output_atom_site = DataCategory(
        "atom_site",
        atom_site_attributes,
    )

    # ------------------------------------------------------------------
    # Helper to retrieve optional attributes
    #
    # Some mmCIF files may not contain every optional attribute.
    # Return "." when an attribute is absent.
    # ------------------------------------------------------------------

    def get_value(atom_site, row, attribute):
        if atom_site.hasAttribute(attribute):
            return row[atom_site.getAttributeIndex(attribute)]

        return "."

    # ------------------------------------------------------------------
    # 11. Add Protein A atoms
    # ------------------------------------------------------------------

    output_atom_id = 1

    for row in atom_site_A.data:

        original_chain = get_value(
            atom_site_A, row, "label_asym_id"
        )

        output_row = [
            get_value(atom_site_A, row, "group_PDB"),
            str(output_atom_id),
            get_value(atom_site_A, row, "type_symbol"),
            get_value(atom_site_A, row, "label_atom_id"),
            get_value(atom_site_A, row, "label_comp_id"),
            chain_map_A[original_chain],
            entity_map_A[original_chain],
            get_value(atom_site_A, row, "label_seq_id"),
            get_value(atom_site_A, row, "pdbx_PDB_ins_code"),
            get_value(atom_site_A, row, "label_alt_id"),
            get_value(atom_site_A, row, "Cartn_x"),
            get_value(atom_site_A, row, "Cartn_y"),
            get_value(atom_site_A, row, "Cartn_z"),
            get_value(atom_site_A, row, "occupancy"),
            get_value(atom_site_A, row, "B_iso_or_equiv"),
        ]

        output_atom_site.append(output_row)
        output_atom_id += 1

    # ------------------------------------------------------------------
    # 12. Add Protein B atoms with transformed coordinates
    # ------------------------------------------------------------------

    for row, (x, y, z) in zip(
        atom_site_B.data,
        transformed_coordinates_B,
    ):

        original_chain = get_value(
            atom_site_B, row, "label_asym_id"
        )

        output_row = [
            get_value(atom_site_B, row, "group_PDB"),
            str(output_atom_id),
            get_value(atom_site_B, row, "type_symbol"),
            get_value(atom_site_B, row, "label_atom_id"),
            get_value(atom_site_B, row, "label_comp_id"),
            chain_map_B[original_chain],
            entity_map_B[original_chain],
            get_value(atom_site_B, row, "label_seq_id"),
            get_value(atom_site_B, row, "pdbx_PDB_ins_code"),
            get_value(atom_site_B, row, "label_alt_id"),
            f"{x:.3f}",
            f"{y:.3f}",
            f"{z:.3f}",
            get_value(atom_site_B, row, "occupancy"),
            get_value(atom_site_B, row, "B_iso_or_equiv"),
        ]

        output_atom_site.append(output_row)
        output_atom_id += 1

    # ------------------------------------------------------------------
    # 13. Add atom_site and write the output file
    # ------------------------------------------------------------------

    output_container.append(output_atom_site)

    io = IoAdapterCore()

    io.writeFile(
        output_file,
        [output_container],
    )
