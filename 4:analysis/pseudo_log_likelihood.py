import torch
import numpy as np
from transformers import EsmForMaskedLM, AutoTokenizer
from Bio import SeqIO


def compute_pll(seq: str, model, tokenizer, device) -> float:
    """Compute the pseudo-log-likelihood scores for each generated sequence"""
    pll = 0.0
    seq = seq.upper()

    for i, aa in enumerate(seq):
        masked_seq = seq[:i] + tokenizer.mask_token + seq[i + 1:]
        inputs = tokenizer(masked_seq, return_tensors="pt").to(device)

        with torch.no_grad():
            outputs = model(**inputs)
        logits = outputs.logits[0]

        mask_idx = (inputs.input_ids[0] == tokenizer.mask_token_id).nonzero(as_tuple=True)[0].item()
        aa_id = tokenizer.convert_tokens_to_ids(aa)
        log_prob = torch.log_softmax(logits[mask_idx], dim=-1)[aa_id]
        pll += log_prob.item()
    
    return pll


def rank_seqs(fasta_path, model_name="facebook/esm2_t33_650M_UR50D"):
    """Rank the generated sequences based on pseudo-log-likelihood scores."""
    
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

    ranked = sorted(results, key=lambda x: x["pll_avg"], reverse=True)
    return ranked
