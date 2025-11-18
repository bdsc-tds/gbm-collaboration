"""Utilities for interacting with Ensembl API."""

import requests
from typing import Optional, Dict


def get_canonical_transcript(gene_id: str) -> str:
    """Get the canonical transcript for a gene ID.
    
    Parameters
    ----------
    gene_id : str
        Ensembl gene ID
        
    Returns
    -------
    str
        Canonical transcript ID or error message
    """
    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?content-type=application/json"
    response = requests.get(url)
    
    if response.status_code == 200:
        gene_data = response.json()
        return gene_data.get("canonical_transcript", "No canonical transcript found")
    else:
        return f"Error: {response.status_code}"


def get_gene_coordinates(gene_name: str, species: str = "homo_sapiens") -> Optional[Dict[str, any]]:
    """Get gene coordinates from Ensembl.
    
    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., "BRAF", "EGFR")
    species : str
        Species name, default "homo_sapiens"
        
    Returns
    -------
    dict or None
        Dictionary with gene_name, ensembl_id, chromosome, start, end, strand
        Returns None if gene not found
    """
    try:
        # Handle gene names with special characters - try the first part if it contains "/"
        query_gene = gene_name.split("/")[0].strip() if "/" in gene_name else gene_name
        
        url = (
            f"https://rest.ensembl.org/xrefs/symbol/{species}/{query_gene}"
            "?content-type=application/json"
        )
        headers = {"Content-Type": "application/json"}
        response = requests.get(url, headers=headers)
        if response.status_code != 200 or not response.json():
            raise ValueError(f"Gene {gene_name} not found in Ensembl")
        
        gene_id = response.json()[0]["id"]
        gene_url = (
            f"https://rest.ensembl.org/lookup/id/{gene_id}"
            "?content-type=application/json"
        )
        gene_data = requests.get(gene_url, headers=headers).json()
        return {
            "gene_name": gene_name,
            "ensembl_id": gene_id,
            "chromosome": "chr" + gene_data["seq_region_name"],
            "start": gene_data["start"],
            "end": gene_data["end"],
            "strand": "forward" if gene_data["strand"] == 1 else "reverse",
        }
    except Exception as e:
        print(f"Warning: Could not get coordinates for {gene_name}: {e}")
        return None

