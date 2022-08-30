import requests

def download_pdb(code: str) -> bool:
    """print Step 1/x: """
    URL = f"https://alphafold.ebi.ac.uk/files/AF-{code}-F1-model_v3.pdb"
    response = requests.get(URL)
    open(f"AF-{code}-F1-model_v3.pdb", "wb").write(response.content)
    return True

download_pdb("Q5VSL9")
