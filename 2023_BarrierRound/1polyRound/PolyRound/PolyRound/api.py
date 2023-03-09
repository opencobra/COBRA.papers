# ©2020-​2021 ETH Zurich, Axel Theorell

import os
import sys

path = os.path.dirname(__file__)
sys.path.insert(0, path)
from PolyRound.mutable_classes.polytope import Polytope
from PolyRound.static_classes.constraint_removal_reduction import PolytopeReducer
from typing import Dict
from PolyRound.static_classes.lp_utils import ChebyshevFinder
from PolyRound.static_classes.rounding.maximum_volume_ellipsoid import (
    MaximumVolumeEllipsoidFinder,
)
from PolyRound.default_settings import (
    default_hp_flags,
    default_0_width,
    default_numerics_threshold,
    default_accepted_tol_violation,
)
import numpy as np
from cobra.core.model import Model
from PolyRound.static_classes.parse_sbml_stoichiometry import StoichiometryParser
from PolyRound.static_classes.hdf5_csv_io import HDF5, CSV

from PolyRound.settings import PolyRoundSettings


class PolyRoundApi:
    @staticmethod
    def simplify_polytope(
        polytope: Polytope, settings: PolyRoundSettings = PolyRoundSettings(),
    ) -> Polytope:
        """
        Remove redundant constraints and refunction inequality constraints to equality constraints in case of dimension
        width less than thresh
        @param polytope:
        @param settings:
        @return:
        """
        polytope = polytope.copy()
        polytope.normalize()
        removed, refunctioned = 1, 1
        while (removed != 0 or refunctioned != 0) and polytope.A.size > 0:
            polytope, removed, refunctioned = PolytopeReducer.constraint_removal(
                polytope, settings,
            )
        if polytope.A.shape[0] == 0:
            raise ValueError(
                "All inequality constraints are redundant, implying that the polytope is a single point."
            )
        return polytope

    @staticmethod
    def transform_polytope(
        polytope: Polytope, settings: PolyRoundSettings = PolyRoundSettings(),
    ) -> Polytope:
        """
        Express polytope in a (shifted) orthogonal basis in the null space of the equality constraints to remove all
        equality constraints
        @param polytope:
        @param settings:
        @return:
        """
        if polytope.inequality_only:
            raise ValueError(
                "Polytope already transformed (only contains inequality constraints)"
            )
        polytope = polytope.copy()
        x, dist = ChebyshevFinder.chebyshev_center(polytope, settings)
        if polytope.border_distance(x) <= 0:
            raise ValueError("Chebyshev center outside polytope before transforming")
        if settings.verbose:
            print("chebyshev distance is : " + str(dist))
            pre_b_dist = polytope.border_distance(x)
            print("border distance pre-transformation is: " + str(pre_b_dist))
        # put x at zero!
        polytope.apply_shift(x)
        if settings.verbose:
            x_0 = np.zeros(x.shape)
            b_dist_at_zero = polytope.border_distance(x_0)
            print("border distance zero-transformation is: " + str(b_dist_at_zero))
        stoichiometry = polytope.S.values
        transformation = PolytopeReducer.null_space(
            stoichiometry, eps=default_numerics_threshold
        )
        polytope.apply_transformation(transformation)
        if settings.verbose:
            u = np.zeros((transformation.shape[1], 1))
            norm_check = np.linalg.norm(np.matmul(stoichiometry, transformation))
            print("norm of the null space is: " + str(norm_check))
            b_dist = polytope.border_distance(u)
            print("border distance after transformation is: " + str(b_dist))
            # test if we can reproduce the original x
            trans_x = polytope.back_transform(u)
            x_rec_diff = np.max(trans_x - np.squeeze(x))
            print("the deviation of the back transform is: " + str(x_rec_diff))
        return polytope

    @staticmethod
    def round_polytope(
        polytope: Polytope, settings: PolyRoundSettings = PolyRoundSettings(),
    ) -> Polytope:
        """
        Round polytope using the maximum volume ellipsoid approach
        @param polytope:
        @param settings:
        @return:
        """
        # check if there are Nans
        bool = False
        bool += np.isinf(polytope.A.values).any()
        bool += np.isinf(polytope.b.values).any()
        if bool:
            raise ValueError("Polytope assigned for rounding contains inf")

        # create a blank polytope so that we can make isolated checks on the rounding transform
        blank_polytope = Polytope(polytope.A, polytope.b)
        MaximumVolumeEllipsoidFinder.iterative_solve(blank_polytope, settings)
        # MaximumVolumeEllipsoidFinder.iterative_solve(
        #     o_polytope, backend, hp_flags=hp_flags, verbose=verbose, sgp=sgp
        # )
        # check if the transformation is full dimensional
        _, s, _ = np.linalg.svd(blank_polytope.transformation)
        if not np.min(s) > default_0_width / default_accepted_tol_violation:
            raise ValueError("Rounding transformation not full dimensional")
        # check if 0 is a solution
        if not blank_polytope.b.min() > 0:
            raise ValueError("Zero point not inside rounded polytope")
        polytope.apply_shift(blank_polytope.shift.values)
        polytope.apply_transformation(blank_polytope.transformation)

        # assert polytope == o_polytope
        return polytope

    @staticmethod
    def simplify_transform_and_round(
        polytope: Polytope, settings: PolyRoundSettings = PolyRoundSettings(),
    ) -> Polytope:
        """
        Conveniently execute simplify_polytope, transform_polytope and round polytope in sequence
        @param polytope:
        @param settings:
        @return:
        """
        polytope = PolyRoundApi.simplify_polytope(polytope, settings=settings,)
        polytope = PolyRoundApi.transform_polytope(polytope, settings=settings,)
        polytope = PolyRoundApi.round_polytope(polytope, settings=settings,)
        return polytope

    @staticmethod
    def cobra_model_to_polytope(model: Model):
        """
        Turn cobrapy model into polytope
        @param model:
        @return:
        """
        return StoichiometryParser.extract_polytope(model)

    @staticmethod
    def polytope_to_hdf5(polytope: Polytope, filename: str):
        HDF5.polytope_to_h5(polytope, filename)

    @staticmethod
    def polytope_to_csvs(polytope: Polytope, dirname: str):
        CSV.polytope_to_csv(polytope, dirname)

    @staticmethod
    def sbml_to_polytope(file_name: str) -> Polytope:
        polytope = StoichiometryParser.parse_sbml_cobrapy(file_name)
        return polytope