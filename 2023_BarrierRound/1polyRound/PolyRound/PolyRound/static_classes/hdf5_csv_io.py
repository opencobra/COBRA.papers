# ©2020-​2021 ETH Zurich, Axel Theorell

import pandas as pd
import os
from PolyRound.mutable_classes.polytope import Polytope
from pathlib import Path
import h5py


class CSV:
    @staticmethod
    def polytope_to_csv(polytope, dirname):
        Path(dirname).mkdir(parents=True, exist_ok=True)
        name = dirname.rstrip("/").split("/")[-1]
        for attribute in dir(polytope):
            tentative_df = getattr(polytope, attribute)
            if isinstance(tentative_df, pd.DataFrame) or isinstance(
                tentative_df, pd.Series
            ):

                if attribute == "transformation":
                    zero_solution_df = pd.Series(0, index=tentative_df.columns)
                    zero_solution_df.to_csv(
                        os.path.join(dirname, "start_" + name + "_rounded.csv"),
                        header=False,
                        index=False,
                    )
                    tentative_df.to_csv(
                        os.path.join(dirname, "N_" + name + "_rounded.csv"),
                        header=False,
                        index=False,
                    )
                elif attribute == "shift":
                    tentative_df.to_csv(
                        os.path.join(dirname, "p_shift_" + name + "_rounded.csv"),
                        header=False,
                        index=False,
                    )
                    name_series = pd.Series(tentative_df.index)
                    name_series.to_csv(
                        os.path.join(
                            dirname, "reaction_names_" + name + "_rounded.csv"
                        ),
                        header=False,
                        index=False,
                    )
                else:
                    tentative_df.to_csv(
                        os.path.join(dirname, attribute + "_" + name + "_rounded.csv"),
                        header=False,
                        index=False,
                    )


class HDF5:
    @staticmethod
    def polytope_to_h5(polytope, filename):
        """
        Writes a Polytope object to an HDF5 file
        :param polytope:
        :param filename:
        :return:
        """

        # hf.create_dataset('start', data=x)

        if os.path.exists(filename):
            os.remove(filename)
        # hf = h5py.File(filename, 'w')
        for attribute in dir(polytope):
            tentative_df = getattr(polytope, attribute)
            if isinstance(tentative_df, pd.DataFrame) or isinstance(
                tentative_df, pd.Series
            ):
                tentative_df.to_hdf(filename, attribute, mode="a")
                # hf.create_dataset(attribute, data=tentative_df)
        #
        # hf.close()

    @staticmethod
    def h5_to_polytope(filename):
        """
        Reads a Polytope object from an HDF5 file. The HDF5 file should have the attributes:
        A (DataFrame), b (Series) Ax < b
        and optionally
        S (DataFrame), h (Series) Sx = h
        :param filename:
        :return: Polytope object
        """
        A = pd.read_hdf(filename, key="A")
        b = pd.read_hdf(filename, key="b")
        try:
            S = pd.read_hdf(filename, key="S")
            h = pd.read_hdf(filename, key="h")
        except KeyError:
            return Polytope(A, b)
        return Polytope(A, b, S=S, h=h)
