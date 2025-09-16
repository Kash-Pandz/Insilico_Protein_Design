import os

def get_res(res_input):
    """ 
    Returns a space-separated string of residues from a string or txt file.

    """
    if not res_input:
        return ""

    if os.path.isfile(res_input):
        with open(res_input, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        res_str = ' '.join(lines)
    else:
        res_str = res_input

    # Create residue ranges
    expanded = []
    for part in res_str.split():
        if '-' in part:
            chain = part[0]
            start, end = map(int, part[1:].split('-'))
            expanded.extend([f"{chain}{i}" for i in range(start, end + 1)])
        else:
            expanded.append(part)
          
    return ' ' .join(expanded)
   

   
