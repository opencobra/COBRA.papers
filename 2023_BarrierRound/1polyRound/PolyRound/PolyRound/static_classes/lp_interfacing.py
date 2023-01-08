# ©2020-​2021 ETH Zurich, Axel Theorell

import numpy as np
import pandas as pd
from scipy import sparse as sp
from optlang.symbolics import Zero
from sympy.core import Add, Mul
from cobra.util.solver import solvers
from swiglpk import glp_adv_basis
from PolyRound.mutable_classes.polytope import Polytope
import uuid
from PolyRound.default_settings import default_accepted_tol_violation
import warnings


check_solutions = True


class OptlangInterfacer:
    @staticmethod
    def make_sparse_constraint_system(m, A, b, x, backend, equality=False):
        interface = solvers[backend]
        A = sp.csr_matrix(A)
        constr_dict = dict()
        for row_ind, bound in enumerate(b):
            if equality:
                lb = bound
            else:
                lb = None
            if isinstance(b, pd.Series):
                uid = b.index[row_ind]
                if not isinstance(uid, str):
                    uid = uuid.uuid4().hex
            else:
                uid = uuid.uuid4().hex

            constr_dict[interface.Constraint(Zero, ub=bound, lb=lb, name=uid)] = {
                x[var_ind]: A[row_ind, var_ind]
                for var_ind in A[row_ind, :].nonzero()[1]
            }

            # constr_dict[interface.Constraint(Zero, ub=bound, lb=lb)] = {
            #     x[var_ind]: A[row_ind, var_ind]
            #     for var_ind in A[row_ind, :].nonzero()[1]
            #     }
        m.add(constr_dict.keys())
        m.update()
        for constr, terms in constr_dict.items():
            constr.set_linear_coefficients(terms)

    @staticmethod
    def configure_gurobi_model(m, settings):
        if not settings.sgp:
            m.problem.setParam("OutputFlag", 0)
        else:
            m.problem.setParam("OutputFlag", 1)
        if settings.verbose:
            print("Using the hp flags: " + str(settings.hp_flags))
        for key, val in settings.hp_flags.items():
            m.problem.setParam(key, val)
        # Never let the solver multi thread
        m.problem.setParam("Threads", 1)

    @staticmethod
    def configure_optlang_model(m, hp_flags, verbose=False):
        if "FeasibilityTol" in hp_flags:
            m.configuration.tolerances.feasibility = hp_flags["FeasibilityTol"]
        if "OptimalityTol" in hp_flags:
            try:
                m.configuration.tolerances.optimality = hp_flags["OptimalityTol"]
            except AttributeError:
                warnings.warn(
                    "Flag OptimalityTol ignored, since backend is glpk and optlang 1.5+ does not support optimality tolerances. To use optimality tolerance, switch backend to gurobi or downgrade optlang to 1.4"
                )

    @staticmethod
    def gurobi_solve(obj, A, b, settings, S=None, h=None):
        interface = solvers[settings.backend]
        m = interface.Model()
        if settings.backend == "gurobi":
            OptlangInterfacer.configure_gurobi_model(m, settings)
        r_names = [str(r) for r in range(A.shape[1])]
        x = np.array([interface.Variable(name, lb=None) for name in r_names])
        m.add(x)
        m.objective = interface.Objective(obj @ x, direction="min")
        # m.setObjective(obj @ x, gp.GRB.MINIMIZE)
        # sp_ratio = np.sum(np.abs(A) > 1e-12) / A.size
        # if sp_ratio < 0.1:
        # A = sp.csr_matrix(A)
        OptlangInterfacer.make_sparse_constraint_system(m, A, b, x, settings.backend)
        # m.addConstr(A @ x <= np.squeeze(b), name="ineq")
        if S is not None:
            assert h is not None
            # sp_ratio = np.sum(np.abs(S) < 1e-12) / S.size
            # if sp_ratio < 0.1:
            # print("Sparse mode")
            # S = sp.csr_matrix(S)
            # m.addConstr(S @ x == np.squeeze(h), name="eq")
            OptlangInterfacer.make_sparse_constraint_system(
                m, S, h, x, settings.backend, equality=True
            )
        # m.update()
        m.optimize()
        if m.status == "optimal":
            return pd.Series(m.primal_values).values, m
        else:
            x = np.zeros(A.shape[1])
            x[:] = np.nan
            return x, m

    @staticmethod
    def gurobi_solve_model(obj, m, backend):
        interface = solvers[backend]
        # m.reset()
        # m.update()
        x = m.variables
        # m.setObjective(obj @ x, gp.GRB.MINIMIZE)
        m.objective = interface.Objective(obj @ x, direction="min")
        # m.update()
        m.optimize()
        if m.status == "optimal":
            return pd.Series(m.primal_values).values, m
        else:
            x = np.zeros(len(m.variables))
            x[:] = np.nan
            return x, m

    @staticmethod
    def gurobi_regularize_chebyshev_center(obj_val, m, backend):
        interface = solvers[backend]
        x = np.array(m.variables)
        last_var = x[-1]
        m.add(interface.Constraint(last_var, lb=obj_val / 2))
        # m.addConstr(last_var, gp.GRB.GREATER_EQUAL, obj_val / 2)
        expr = x @ x
        m.objective = interface.Objective(expr, direction="min")
        # obj = np.eye(len(x))
        # obj[-1, -1] = 0
        # m.setMObjective(obj, None, 0, sense=gp.GRB.MINIMIZE)
        m.update()
        m.optimize()
        return pd.Series(m.primal_values).values, m

    @staticmethod
    def polytope_to_optlang(polytope, settings):
        interface = solvers[settings.backend]
        m = interface.Model()
        m.configuration.presolve = settings.presolve
        r_names = [str(r) for r in polytope.A.columns]
        x = np.array([interface.Variable(name, lb=None) for name in r_names])
        m.add(x)

        # make all inequality constraints
        OptlangInterfacer.make_sparse_constraint_system(
            m, polytope.A, polytope.b, x, settings.backend
        )
        # m.addConstr(A @ x <= np.squeeze(polytope.b.values), name="constr")

        if polytope.S is not None:
            # S = sp.csr_matrix(polytope.S)
            # m.addConstr(S @ x == np.squeeze(polytope.h.values), name="eq")
            OptlangInterfacer.make_sparse_constraint_system(
                m, polytope.S, polytope.h, x, settings.backend, equality=True
            )
        m.update()
        return m

    @staticmethod
    def optlang_to_polytope(m):
        A, b = OptlangInterfacer.constraints_as_mat(m, sense="<")
        S, h = OptlangInterfacer.constraints_as_mat(m, sense="=")
        if S.size > 0:
            p = Polytope(A, b, S, h)
        else:
            p = Polytope(A, b)
        return p

    @staticmethod
    def constraints_as_mat(m, sense="<"):
        r_names = [v.name for v in m.variables]
        constrs = m.constraints
        if sense == "<":
            constrs = [c for c in constrs if c.ub != c.lb]
        elif sense == "=":
            constrs = [c for c in constrs if c.ub == c.lb]
        else:
            raise ValueError
        c_array = np.array([(constr.name, constr.ub) for constr in constrs]).T
        if c_array.size == 0:
            return pd.DataFrame(dtype=np.float64), pd.Series(dtype=np.float64)
        b = pd.Series(c_array[1, :], index=c_array[0, :], dtype=np.float64)
        c_df = pd.DataFrame(index=c_array[0, :], columns=r_names, dtype=np.float64)
        for constr in constrs:

            expr_dict = OptlangInterfacer.expr_to_dict(constr.expression)
            for var, coeff in expr_dict.items():
                c_df.loc[constr.name, var.name] = float(coeff)
            # if isinstance(constr.expression, Add):
            #     for mul in constr.expression.args:
            #         assert isinstance(mul, Mul)
            #         mul_dict = mul.as_coefficients_dict()
            #         for var, coeff in mul_dict.items():
            #             c_df.loc[constr.name, var.name] = float(coeff)
            # elif isinstance(constr.expression, Mul):
            #     mul_dict = constr.expression.as_coefficients_dict()
            #     for var, coeff in mul_dict.items():
            #         c_df.loc[constr.name, var.name] = float(coeff)
            # else:
            #     raise ValueError(
            #         "Invalid constraint expression detected, probably caused by a zero row in a constraint matrix."
            #     )
        c_df.fillna(0.0, inplace=True)
        return c_df, b

    @staticmethod
    def expr_to_dict(expr):
        expr_dict = dict()
        if isinstance(expr, Add):
            for mul in expr.args:
                assert isinstance(mul, Mul)
                mul_dict = mul.as_coefficients_dict()
                expr_dict.update(mul_dict)
        elif isinstance(expr, Mul):
            expr_dict = dict(expr.as_coefficients_dict())
        else:
            raise ValueError(
                "Invalid constraint expression detected, probably caused by a zero row in a constraint matrix."
            )
        return expr_dict

    @staticmethod
    def check_tolerances(m):
        # get solution
        # solution = pd.Series(m.primal_values).values
        # # get absolute errors from original system
        # constr_vals = np.matmul(polytope.A.values, solution).T
        # # temp_b = polytope.b.copy().values
        # perturbed = False
        # for ind, c_id in enumerate(polytope.b.index):
        #     if c_id in m.constraints:
        #         if polytope.b[c_id] != m.constraints[c_id].ub:
        #             perturbed = True
        #             dist = np.nan
        #             rel_dist = np.nan
        #             break
        #
        #         # temp_b[ind] = m.constraints[c_id].ub
        # if not perturbed:
        #     dists = (polytope.b.values - constr_vals).T
        #     dist = np.min(dists)
        #     # get relative errors
        #     bigger_than_one = np.abs(constr_vals) > 1
        #     rel_dist = np.min(dists[(bigger_than_one)]/np.abs(constr_vals[bigger_than_one]))
        # get constraint violations
        ub_viol = np.array(
            [
                constr.ub - constr.primal
                for constr in m.constraints
                if constr.ub is not None
            ]
        )
        lb_viol = np.array(
            [
                -constr.lb + constr.primal
                for constr in m.constraints
                if constr.lb is not None
            ]
        )
        viol = np.concatenate([ub_viol, lb_viol])
        worst = np.min(viol)
        thresh = m.configuration.tolerances.feasibility * default_accepted_tol_violation
        if (-worst) > thresh:  # or ((-dist) > thresh and not np.isinf(dist)):
            raise ValueError(
                "Feasibility tolerance violated"
            )  # . The relative violation is: " + str(rel_dist))
            # p_new = OptlangInterfacer.optlang_to_polytope(m)
            # dists_new = (p_new.b.values - np.matmul(p_new.A.values, solution).T).T
            # dist_new = np.min(dists_new)
            # logging.warning(', '.join([str(worst), str(dist), str(rel_dist)]))

        # get optimality tolerance
        r_costs = np.array(list(m.reduced_costs.values()))
        sense = float(m.objective.direction == "max") * 2 - 1
        opt_violation = np.max(r_costs * sense)

        if (
            opt_violation
            > m.configuration.tolerances.optimality * default_accepted_tol_violation
        ):
            raise ValueError("Optimality tolerance violated")
            # logging.warning(opt_violation)

    @staticmethod
    def get_opt(m, settings):
        if m.status == "optimal":
            if settings.check_lps:
                OptlangInterfacer.check_tolerances(m)
            return m.objective.value
        elif (
            m.status == "infeasible"
            or m.status == "check_original_solver_status"
            or m.status == "undefined"
        ):
            print("model in infeasible state, resetting lp")
            if settings.backend == "gurobi":
                m.problem.reset()
            elif settings.backend == "glpk" or settings.backend == "glpk_exact":
                glp_adv_basis(m.problem, 0)
            else:
                raise ValueError(
                    "PolyRound does currently not support system reset for backend "
                    + settings.backend
                )
            m.optimize()
            # m.setParam("BarHomogeneous", -1)
            if m.status == "optimal":
                return m.objective.value
            else:
                print("Solver status: " + str(m.status))
                raise ValueError("Optimization fails despite resetting")
        else:
            print("Solver status: " + str(m.status))
            raise ValueError
