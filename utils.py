import os
import pandas as pd
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import shutil
from Bio.PDB import PDBParser, PPBuilder
from Bio import PDB


# get key from original pdbs
def get_key(name):
    key = name.split('.')[0]
    return key

# get key from AFig results
def get_key_AF(name: str) -> str:
    return re.sub(r'_[^_]+_af2pred\.pdb$', '', name)

# get TMscore
def get_tmscore(pdb1, pdb2):
    """运行 USalign 并返回第二个 TM-score（标准值）"""
    result = subprocess.run(
        ["USalign", pdb1, pdb2, "-mm", "1", "-ter", "0"],
        capture_output=True, text=True
    )
    # 正则匹配所有 TM-score 行
    scores = re.findall(r"TM-score=\s*([\d.]+)", result.stdout)
    if scores:
        # 返回第二个（normalized by Structure_2）
        return float(scores[-1])
    else:
        return None

