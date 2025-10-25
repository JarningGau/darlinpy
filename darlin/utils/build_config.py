import json

def build_carlin_config(sequence, output_path="carlin_config.json"):
    """
    Build CARLIN config JSON where:
      - open/close penalties are generated dynamically from sequence lengths
      - match_score.* = len(sequence) - 1  (per the updated MATLAB rule)
    """
    def arange_inclusive(start, stop, step):
        vals, x = [], start
        if step > 0:
            while x <= stop + 1e-9:
                vals.append(round(x, 3))
                x += step
        else:
            while x >= stop - 1e-9:
                vals.append(round(x, 3))
                x += step
        return vals

    # lengths
    seg_len = len(sequence["segments"][0])
    if any(len(s) != seg_len for s in sequence["segments"]):
        raise ValueError("All segments must have the same length.")
    pam_len = len(sequence["pam"])
    prefix_len = len(sequence["prefix"])
    postfix_len = len(sequence["postfix"])

    # constants derived from MATLAB
    cutsite_len = 7
    consite_len = seg_len - cutsite_len
    if consite_len <= 0:
        raise ValueError("segment length must be greater than cutsite length (7).")

    # penalties
    open_penalty = {
        "init": 10,
        "cutsites": [8, 7.5, 7.0, 6.5] + [6] * (cutsite_len - 4),
        "consites": arange_inclusive(10, 8.5, -0.125),
        "pam": arange_inclusive(8.6, 8.6 + 0.2 * (pam_len - 1), 0.2),
        "prefix": [10] * prefix_len,
        "postfix": arange_inclusive(8.6, 8.6 + 0.2 * (postfix_len - 1), 0.2),
    }

    close_penalty = {
        "cutsites": [8, 7.5, 7.0, 6.5] + [6.5] * (cutsite_len - 4),
        "consites": arange_inclusive(10, 8.5, -0.125),
        "pam": arange_inclusive(8.6, 8.6 + 0.2 * (pam_len - 1), 0.2),
        "prefix": [10] * prefix_len,
        "postfix": arange_inclusive(8.6, 8.6 + 0.2 * (postfix_len - 1), 0.2),
    }

    # match scores per new MATLAB rule: len(seq) - 1
    match_score = {
        "Primer5": len(sequence["Primer5"]) - 1,
        "Primer3": len(sequence["Primer3"]) - 1,
        "SecondarySequence": len(sequence["SecondarySequence"]) - 1,
    }

    config = {
        "sequence": sequence,
        "match_score": match_score,
        "open_penalty": open_penalty,
        "close_penalty": close_penalty,
    }

    with open(output_path, "w") as f:
        json.dump(config, f, indent=4)

    return output_path

# Build a file using the user's earlier sequence as an example
# sequence = {
#     "segments": [
#         "GACTGCACGACAGTCGACGA",
#         "GACACGACTCGCGCATACGA",
#         "GACTACAGTCGCTACGACGA",
#         "GCGAGCGCTATGAGCGACTA",
#         "GATACGATACGCGCACGCTA",
#         "GAGAGCGCGCTCGTCGACTA",
#         "GCGACTGTACGCACACGCGA",
#         "GATAGTATGCGTACACGCGA",
#         "GAGTCGAGACGCTGACGATA",
#         "GATACGTAGCACGCAGACGA",
#     ],
#     "pam": "TGGAGTC",
#     "prefix": "CGCCG",
#     "postfix": "TGGGAGCT",
#     "Primer5": "GAGCTGTACAAGTAAGCGGC",
#     "Primer3": "CGACTGTGCCTTCTAGTTGC",
#     "SecondarySequence": "AGAATTCTAACTAGAGCTCGCTGATCAGCCT",
# }

# out_path = build_carlin_config(sequence, output_path="/mnt/data/carlin_config_v2.json")