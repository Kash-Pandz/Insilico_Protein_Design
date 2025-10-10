import torch
import numpy as np
import pandas as pd
from transformers import EsmForMaskedLM, AutoTokenizer
from Bio import SeqIO

def compute_pll(seq: str, model, tokenizer, device) -> float:
    """Compute pseudo-log-likelihood for a protein sequence."""
    seq = seq.upper()
    valid_aas = set(list("ACDEFGHIKLMNPQRSTVWY"))
    seq = "".join([aa if aa in valid_aas else "X" for aa in seq])
    
    pll = 0.0

    for i, aa in enumerate(seq):
        masked_seq = seq[:i] + tokenizer.mask_token + seq[i + 1:]
        inputs = tokenizer(masked_seq, return_tensors="pt").to(device)

        with torch.no_grad():
            outputs = model(**inputs)
        logits = outputs.logits[0]

        mask_idx = (inputs.input_ids[0] == tokenizer.mask_token_id).nonzero(as_tuple=True)[0].item()
        aa_id = tokenizer.convert_tokens_to_ids(aa)
        log_prob = torch.log_softmax(logits[mask_idx], dim=-1)[aa_id]
        pll += float(log_prob)
    
    return pll

def rank_seqs(fasta_path, model_name="facebook/esm2_t33_650M_UR50D"):
    """Rank protein sequences by average PLL using ESM2."""
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmForMaskedLM.from_pretrained(model_name)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device).eval()

    records = list(SeqIO.parse(fasta_path, "fasta"))
    results = []

    for record in records:
        seq = str(record.seq)
        pll = compute_pll(seq, model, tokenizer, device)
        pll_avg = pll / len(seq)
        results.append({
            "id": record.id,
            "sequence": seq,
            "pll": pll,
            "pll_avg": pll_avg
        })

    df = pd.DataFrame(results)
    df["rank"] = df["pll_avg"].rank(ascending=False, method="dense").astype(int)
    df = df.sort_values("pll_avg", ascending=False).reset_index(drop=True)

    return df
