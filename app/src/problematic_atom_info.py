from rdkit.Chem.rdMolTransforms import GetBondLength

def problematic_atom_info(rdkit_at,
                          at_ba,
                          rdkit_mol):
    # check hydrogens for clash

    message = f"The bonding environment {', '.join(at_ba.split('/')[1])} of this atom is non-standard for proteins. "


    for bond in rdkit_at.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        a1_symbol = a1.GetSymbol()
        a2_symbol = a2.GetSymbol()
        if "H" in [a1_symbol, a2_symbol]:
            bond_length = GetBondLength(rdkit_mol.GetConformer(), a1.GetIdx(), a2.GetIdx())
            if bond_length < 0.9:
                if a1_symbol == "H":
                    clashed_hydrogen = a1
                else:
                    clashed_hydrogen = a2
                message += f"Hydrogen {clashed_hydrogen.GetPDBResidueInfo().GetResidueName()} {clashed_hydrogen.GetPDBResidueInfo().GetResidueNumber()}{clashed_hydrogen.GetPDBResidueInfo().GetName().rstrip()} is probably too close to this atom ({round(bond_length, 2)}Å). The problem occurs during protonation using the external tool pdb2pqr."
                break
    else:
        message += "This atom is probably wrongly predicted by AlphaFold2."

    print(message)
    return message

