from rdkit.Chem import Mol, rdShapeHelpers
from itertools import product
from scipy.spatial.distance import cdist
import numpy as np
import mrich
from mrich import print


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

def load_sigFactory(fdef_file=config.FINGERPRINT_FDEF, max_point_count=config.FINGERPRINT_MAXPOINTCOUNT,
                    bins=config.FINGERPRINT_BINS):
    """
    Load signature factory for pharmacophore fp calculation

    :param fdef_file:
    :param max_point_count:
    :param bins:
    :return:
    """
    featFactory = ChemicalFeatures.BuildFeatureFactory(fdef_file)
    sigFactory = SigFactory(featFactory, maxPointCount=max_point_count)
    sigFactory.SetBins(bins)
    sigFactory.Init()
    sigFactory.GetSigSize()
    return sigFactory

def calc_pharm_fp(mol, sigFactory, asStr=True):
    """
    Calculate pharmacophore fingerprint using RDKit

    :param mol:
    :param sigFactory:
    :param asStr:
    :return:
    """
    fp = Generate.Gen2DFingerprint(mol, sigFactory)
    if asStr:
        fp = list(fp)
        return ";".join(map(str, fp))
    else:
        return list(fp)
