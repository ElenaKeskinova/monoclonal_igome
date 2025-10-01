from pymol import cmd
import os

file = "D:/Elena/ban/monoclonal_igome/dock_motifier_mimotopes/allpep_models/allpeps_21c.txt"
infolder = "D:/Elena/ban/monoclonal_igome/dock_motifier_mimotopes/allpep_models/21c"


with open(file, "r") as f:
    lines = f.readlines()

for i, seq in enumerate(lines, start=1):
    seq = seq.strip()
    

    # Build peptide object
    obj_name = seq
    cmd.fab(seq, obj_name)

    # Save as PDB
    cmd.save(f"{infolder}/{obj_name}.pdb", obj_name)
    cmd.delete("all")