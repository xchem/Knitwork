import mrich
from mrich import print

import numpy as np
from pathlib import Path
from itertools import product
from scipy.spatial.distance import cdist

from rdkit.Chem import Mol, rdShapeHelpers
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory


def pair_overlap(molA: Mol, molB: Mol):
    overlapA = 1 - rdShapeHelpers.ShapeProtrudeDist(molA, molB, allowReordering=False)
    overlapB = 1 - rdShapeHelpers.ShapeProtrudeDist(molB, molA, allowReordering=False)
    return max([overlapA, overlapB])


def pair_min_distance(molA: Mol, molB: Mol):
    confA = molA.GetConformer()
    confB = molB.GetConformer()
    coords1 = np.array(
        [list(confA.GetAtomPosition(i)) for i in range(molA.GetNumAtoms())]
    )
    coords2 = np.array(
        [list(confB.GetAtomPosition(i)) for i in range(molB.GetNumAtoms())]
    )
    dists = cdist(coords1, coords2)
    # print(coords1)
    # print(coords2)
    # print(dists)
    # raise NotImplementedError
    return np.min(dists)


def load_sig_factory(
    fdef_file: str | Path,
    max_point_count: int,
    bins: list[list[int]],
):

    # Get FDef file path

    fdef_path = Path(fdef_file)

    if not fdef_path.exists():
        fdef_path = Path(__file__).parent / fdef_file

    assert fdef_path.exists(), f"Can't find fdef_file: {fdef_file}"

    # feature factory

    feature_factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)

    # sig_factory

    sig_factory = SigFactory(feature_factory, maxPointCount=max_point_count)
    sig_factory.SetBins(bins)
    sig_factory.Init()
    sig_factory.GetSigSize()

    return sig_factory


def calc_pharm_fp(mol, sig_factory, as_str=True):
    """
    Calculate pharmacophore fingerprint using RDKit

    :param mol:
    :param sigFactory:
    :param as_str:
    :return:
    """
    fp = Generate.Gen2DFingerprint(mol, sig_factory)
    if as_str:
        fp = list(fp)
        return ";".join(map(str, fp))
    else:
        return list(fp)
