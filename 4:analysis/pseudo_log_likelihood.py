import torch
import pandas as pd
from transformers import EsmForMaskedLM, AutoTokenizer
from Bio import SeqIO

def compute_pll(seq, model, tokenizer, device, batch_size=64):
    seq = seq.upper()
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    seq = "".join([aa if aa in valid_aas else "X" for aa in seq])
    masked_seqs = [seq[:i] + tokenizer.mask_token + seq[i + 1:] for i in range(len(seq))]
    pll = 0.0

    for i in range(0, len(masked_seqs), batch_size):
        batch = masked_seqs[i:i + batch_size]
        inputs = tokenizer(batch, return_tensors="pt", padding=True, truncation=True).to(device)
        with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16):
            outputs = model(**inputs)
        logits = outputs.logits
        for j in range(len(batch)):
            mask_idx = (inputs.input_ids[j] == tokenizer.mask_token_id).nonzero(as_tuple=True)[0].item()
            aa = seq[i + j]
            aa_id = tokenizer.convert_tokens_to_ids(aa)
            log_prob = torch.log_softmax(logits[j, mask_idx].float(), dim=-1)[aa_id]
            pll += float(log_prob)
    return pll

def rank_seqs(fasta_path, model_name="facebook/esm2_t33_650M_UR50D", batch_size=64):
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmForMaskedLM.from_pretrained(model_name, torch_dtype=torch.bfloat16, device_map="auto")
    try:
        model = torch.compile(model)
    except Exception:
        pass
    device = torch.device("cuda")
    model.eval()
    torch.backends.cudnn.benchmark = True

    records = list(SeqIO.parse(fasta_path, "fasta"))
    results = []
    for record in records:
        seq = str(record.seq)
        pll = compute_pll_batched_h100(seq, model, tokenizer, device, batch_size)
        pll_avg = pll / len(seq)
        results.append({"id": record.id, "sequence": seq, "pll": pll, "pll_avg": pll_avg})

    df = pd.DataFrame(results)
    df["rank"] = df["pll_avg"].rank(ascending=False, method="dense").astype(int)
    df = df.sort_values("pll_avg", ascending=False).reset_index(drop=True)
    return df

if __name__ == "__main__":
    fasta_file = "generated_sequences.fasta"
    df = rank_seqs(fasta_file, batch_size=64)
    df.to_csv("pll_ranked_sequences.csv", index=False)

 
