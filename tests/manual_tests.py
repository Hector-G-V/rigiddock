# There are no functions here that can be called. Instead, these are scripts 
# that can be copy-pasted to run.py to test different stuff.


###
# Tests
###

# Manual test - old code:
#tests.utils.plot_output(a_c, b_c_pdb, e_angles, v_shift) # Works, but uses old code

# Manual test:
"""vs_a, a_bar = p.p_bar.out(a_c, c.RHO, constants)
vs_a, a_bar = p.p_bar.add_thickness(a_bar, vs_a, c.RHO)
vs_a, a_bar = p.p_bar.add_thickness(a_bar, vs_a, c.RHO) # Testing disc stage
vs_b, b_bar = p.p_bar.out(b_c_pdb, c.DEL_B, constants)
b_r = rotation.utils.e_rotate(e_angles, b_c_pdb) + v_shift
vs_br, br_bar = p.p_bar.out(b_r, c.DEL_B, constants)
tests.utils.plot_protein_voxels([a_bar, b_bar, br_bar],[vs_a, vs_b, vs_br])
tests.utils.plot_protein_voxels_s([a_bar, b_bar, br_bar],[vs_a, vs_b, vs_br])
# """

# Plot tests
"""tests.utils.plot_protein_voxels_l([a, a-vs_a]) # Now relegated to testing code; p is no longer forward-facing
tests.utils.plot_protein_voxels([a_bar, b_bar],[vs_a,vs_b])
tests.utils.plot_protein_voxels_s([a_bar, b_bar],[vs_a,vs_b])
tests.utils.plot_protein_voxels_s([a_bar, b_bar, b_bar],[vs_a, vs_b, vs_a - v_k])
#   Note: vs_b excluded because b_grid is already at Octant 1"""

# Scan stage tests
"""arr = tests.utils.to_python('scan_out.txt') # Retrieve scan stage data from txt
tests.utils.plot_1D([i[0] for i in arr]) # Plot scan stage data 
# """

# Compare PyMol output to Python manual test. So far, the same: Write .cif to .py
"""arr = tests.utils.to_python('b_xyz.txt')
tests.utils.to_txt(arr)
# """


# ##
# Plot real proteins!
# ##

"""# Shrink by plotting only a_c voxels that are within radius c.R

# Compute the threshold squared distance
threshold = (1.8 + 2.2) ** 2 # Approx R + added surface layer thickness

# Compute squared norms
a_c_squared = np.sum(a_c**2, axis=1)[:, np.newaxis]  # (n, 1)
b_c_pdb_squared = np.sum(b_c_pdb**2, axis=1)  # (m,)

# Compute pairwise squared distances
distances_squared = a_c_squared + b_c_pdb_squared - 2 * np.dot(a_c, b_c_pdb.T)

# Ensure non-negative distances (due to precision issues)
distances_squared = np.maximum(distances_squared, 0)

# Determine which a_c rows satisfy the condition with any b_c_pdb point
satisfies_condition = np.any(distances_squared <= threshold, axis=1)

# Filter the rows of a_c that meet the condition
new_a_c = a_c[satisfies_condition]

print()
print('len new_a_c: ', len(new_a_c),'\n')
print('smallest distance: ', np.sqrt(np.min(distances_squared)),'\n')
print(np.sqrt(distances_squared[satisfies_condition]))

# Now plot
constants = utils.set_constants(2)
vs_a, a_bar = p.p_bar.out(new_a_c, c.RHO, constants)
vs_a, a_bar = p.p_bar.add_thickness(a_bar, vs_a, c.RHO) # Add thickness to 'a'
vs_a, a_bar = p.p_bar.add_thickness(a_bar, vs_a, c.RHO) # Add thickness to 'a'
vs_b, b_bar = p.p_bar.out(b_c_pdb, c.DEL_B, constants)
tests.utils.plot_protein_voxels_s([a_bar, b_bar],[vs_a,vs_b])
# """
