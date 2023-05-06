#
# Copyright (c) 2021, ATA Engineering, Inc.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import pyomo.opt
from pyomo.core.util import sum_product
from pyomo.environ import Objective, Constraint, Expression, Var, Param, ConcreteModel, SolverFactory, RangeSet,\
    Reals, PositiveReals
import math


class TrajectoryResults:

    def __init__(self, results) -> None:
        self.results = results


class TrajectoryGenerator:
    def __init__(self, distance, dt, mass, velocity_bounds, acceleration_bounds, jerk_bounds, power_bounds, tf=-1) -> None:
        self.velocity_bounds = velocity_bounds
        self.acceleration_bounds = acceleration_bounds
        self.jerk_bounds = jerk_bounds
        self.power_bounds = power_bounds
        self.distance = distance
        self.mass = mass
        self.dt = dt
        self.tf = tf
        self.time = [dt * i for i in range(0, int(tf / dt) + 1)]
        self.optimization = self.OptimizationInitialize()

    @staticmethod
    def _VariableBounds(minimum, maximum):

        """ This function sets the bounds on a given variable. It first checks to see if vectors were provided and
        then returns the appropriate bounds. This is in accordance with the Pyomo documentation for setting bounds n
        variables.

        if a lower bound (lb) and an upper bound (ub) are provided, then the function returns a lb,
        ub which corresponds to the constraint

                    lb <= x <= ub
        """

        def f(model):
            lb = True
            ub = True
            try:
                minimum
            except:
                lb = False
            try:
                maximum
            except:
                ub = False
            result = minimum if lb else None, maximum if ub else None
            return result

        return f

    @staticmethod
    def _EqualityConstr(var, const, offset=0):
        """
        This is a general wrapper for initial constraints. the index set i needs to be the index of one set of the
        variables that needs to be initialized.

        This will set a constraint so that
                                [[x_0],          [[x_0(0)],
                                  ...,      =     ...,
                                 [x_n-1]]         [x_n-1(0)]]

        :param var:
        :param initVar:
        :return: constraint
        """

        def f(model, i):
            return var[i] == const[i - offset]

        return f

    def _DynamicConstr(self, var, var_deriv, dt):
        """
        This is a general wrapper for initial constraints. the index set i needs to be the index of one set of the
        variables that needs to be initialized.

        This will set a constraint so that
                                [[x_0],          [[x_0(0)],
                                  ...,      =     ...,
                                 [x_n-1]]         [x_n-1(0)]]

        :param var:
        :param initVar:
        :return: constraint
        """

        def f(model, i):
            return var[i] == var[i - 1] + var_deriv[i - 1] * dt

        return f

    def _PeakConst(self, var, max):
        def f(model, i):
            return var[i] <= max

        return f

    def _UpperConstr(self, var, val):
        def f(model, i):
            return var[i] <= val

        return f

    def _LowerConstr(self, var, val):
        def f(model, i):
            return var[i] >= val

        return f

    def _UpperAbsConstr(self, var, abs_var):
        def f(model, i):
            return var[i] <= abs_var[i]

        return f

    def _LowerAbsConstr(self, var, abs_var):
        def f(model, i):
            return -var[i] <= abs_var[i]

        return f

    def _Objective(self, norms, vars, weights):

        # Error checking on input lengths
        if len(norms) != len(vars) != len(weights):
            raise Exception("Error in GenerateTrajectory: input lengths must all be the same")

        def f(model):
            objective = 0

            for iobj, (norm, var, weight) in enumerate(zip(norms, vars, weights)):

                if norm == 'peak':  # Peak
                    name1 = 'abs_var' + str(iobj)
                    name2 = 'upper_abs_constraint' + str(iobj)
                    name3 = 'lower_abs_constraint' + str(iobj)
                    name4 = 'max' + str(iobj)
                    name5 = 'PeakConstraint' + str(iobj)

                    # Define Xk for each variable
                    model.add_component(name1, Var(model.N, domain=Reals))

                    # Define each constraint of the form xk <= Xk <= xk
                    model.add_component(name2, Constraint(model.N, rule=self._UpperAbsConstr(getattr(model, var),
                                                                                             getattr(model, name1))))
                    model.add_component(name3, Constraint(model.N, rule=self._LowerAbsConstr(getattr(model, var),
                                                                                             getattr(model, name1))))
                    # Instantiate the max variable
                    # WILL NEED TO MAX THIS THE NUMBER OF MAX VARIABLES FOR NUMBER OF MAX IN MULTI-OBJECTIVE use the
                    # getattr and setattr functions
                    model.add_component(name4, Var([1], domain=Reals))

                    # Create the peak variable constraints
                    model.add_component(name5, Constraint(model.N, rule=self._PeakConst(getattr(model, name1),
                                                                                        getattr(model, name4)[1])))

                    # Add the max variable to the total objective
                    objective += weight * getattr(model, name4)[1]

                elif norm == 'rms':  # Root-mean-square
                    objective += weight * sum(getattr(model, var)[i] ** 2 for i in model.N)

                elif norm == 'abs':  # Absolute value
                    name1 = 'abs_var' + str(iobj)
                    name2 = 'upper_abs_constraint' + str(iobj)
                    name3 = 'lower_abs_constraint' + str(iobj)

                    # Define Xk for each variable
                    model.add_component(name1, Var(model.N, domain=Reals))

                    # Define each constraint of the form xk <= Xk <= xk
                    model.add_component(name2, Constraint(model.N, rule=self._UpperAbsConstr(getattr(model, var),
                                                                                             getattr(model, name1))))
                    model.add_component(name3, Constraint(model.N, rule=self._LowerAbsConstr(getattr(model, var),
                                                                                             getattr(model, name1))))

                    # Define the sum of the absolute values
                    objective += weight * sum(getattr(model, name1)[i] for i in model.N)

            return objective

        return f

    @staticmethod
    def _timeObj():
        def f(model):
            return model.dt[1]

        return f

    @staticmethod
    def _energyObj(dt):
        def f(model):
            # Can't do this here because it evaluates the expression
            # val = TrajectoryGenerator._EnergyCalculation(model.v, model.a, dt)

            # This just creates the expression
            val = sum_product(abs(model.a), model.v) * model.dt
            return val

        return f

    @staticmethod
    def _ExtractData(variable, index_set):
        data = []

        # If it is an expression then we need to evaluate it as such
        if type(variable) == pyomo.core.base.expression.IndexedExpression:
            for i in index_set:
                val = pyomo.core.expr.current.evaluate_expression(variable[i])
                if val is None:
                    data.append(0)
                else:
                    data.append(val)

        # Otherwise it is a variable and we can just get the value
        else:
            for i in index_set:
                if variable[i].value is None:
                    data.append(0)
                else:
                    data.append(variable[i].value)

        return data

    @staticmethod
    def _PowerExpression(model, i):
        return model.a[i] * model.v[i] * model.mass

    @staticmethod
    def _TotalEnergyExpression(model, i):
        val = sum_product(model.a, model.v)*model.mass*model.dt
        return val

    def OptimizationInitialize(self):
        # This initializes the model and sets as many of the system parameters and variables as possible.
        v_min = self.velocity_bounds[0]
        v_max = self.velocity_bounds[1]
        a_min = self.acceleration_bounds[0]
        a_max = self.acceleration_bounds[1]
        j_min = self.jerk_bounds[0]
        j_max = self.jerk_bounds[1]
        p_min = self.power_bounds[0]
        p_max = self.power_bounds[1]
        N = int(self.tf / self.dt)
        D = self.distance
        dt = self.dt
        mass = self.mass
        optProg = ConcreteModel()

        # Initialize the RangeSet that defines the number
        optProg.N = RangeSet(0, math.ceil(self.tf / self.dt))
        optProg.Nf = RangeSet(1, math.ceil(self.tf / self.dt))

        # Parameters
        optProg.dt = Param(initialize=dt)
        optProg.mass = Param(initialize=mass)

        # Variables
        optProg.d = Var(optProg.N, domain=Reals)
        optProg.v = Var(optProg.N, domain=Reals, bounds=self._VariableBounds(v_min, v_max))
        optProg.a = Var(optProg.N, domain=Reals, bounds=self._VariableBounds(a_min, a_max))
        optProg.j = Var(optProg.N, domain=Reals, bounds=self._VariableBounds(j_min, j_max))
        optProg.p = Expression(optProg.N, rule=self._PowerExpression)
        optProg.e = Expression([1], rule=self._TotalEnergyExpression)

        # t0 and tf constraints
        optProg.d_0 = Constraint([0], rule=self._EqualityConstr(optProg.d, [0]))
        optProg.v_0 = Constraint([0], rule=self._EqualityConstr(optProg.v, [0]))
        optProg.a_0 = Constraint([0], rule=self._EqualityConstr(optProg.a, [0]))
        optProg.d_N = Constraint([N], rule=self._EqualityConstr(optProg.d, [D], N))
        optProg.v_N = Constraint([N], rule=self._EqualityConstr(optProg.v, [0], N))
        optProg.a_N = Constraint([N], rule=self._EqualityConstr(optProg.a, [0], N))

        # Dynamic constraints to ensure dynamic compatibility between d, v, a, j
        optProg.d_k = Constraint(optProg.Nf, rule=self._DynamicConstr(optProg.d, optProg.v, self.dt))
        optProg.v_k = Constraint(optProg.Nf, rule=self._DynamicConstr(optProg.v, optProg.a, self.dt))
        optProg.a_k = Constraint(optProg.Nf, rule=self._DynamicConstr(optProg.a, optProg.j, self.dt))

        # General constraints
        optProg.p_constrmax = Constraint(optProg.N, rule=self._UpperConstr(optProg.p, p_max))
        optProg.p_constrmin = Constraint(optProg.N, rule=self._LowerConstr(optProg.p, p_min))

        return optProg

    def GenerateTrajectory(self, norm=['peak'], var=['j'], weights=[1]):
        """
        This generates a trajectory for a given norm and variable. Weights has been added for extension to
        multi-objective in later work.

        Valid values for the keyword arguments are:

         norm = peak, rms, abs
         var = v, a, j (velocity, acceleration, and jerk repsectively)

        All values are meant to be lists. This is not required not but will be used later for additional functionality

        :param norm:
        :param var:
        :param weights:
        :return: Vectors for time, distance, velocity, acceleration, and jerk
        """

        solProg = self.optimization.clone()

        if var[0] == 't':
            solProg.dt = Var([1], domain=PositiveReals, initialize=self.dt, bounds=self._VariableBounds(0, 0.1))
            solProg.obj = Objective(rule=self._timeObj())

            # Redefine dynamics constraints to only use the first time step (Why?)
            solProg.d_k = Constraint(solProg.Nf, rule=self._DynamicConstr(solProg.d, solProg.v, solProg.dt[1]))
            solProg.v_k = Constraint(solProg.Nf, rule=self._DynamicConstr(solProg.v, solProg.a, solProg.dt[1]))
            solProg.a_k = Constraint(solProg.Nf, rule=self._DynamicConstr(solProg.a, solProg.j, solProg.dt[1]))

        elif var[0] == 'e':
            # The total energy objective is equivalent to the 1-norm of power multiplied by dt
            solProg.obj = Objective(rule=self._Objective(['abs'], ['p'], [self.dt]))

        else:
            solProg.obj = Objective(rule=self._Objective(norm, var, weights))

        # IPOPT doesn't like the time optimization with the power constraint, but works fine without the power
        # constraint and for the other objective functions. This is probably because of the highly nonlinear nature of
        # the power constraint when delta_t is changing, which I believe is not a strength of IPOPT
        # opt = SolverFactory('ipopt')
        # opt.options['bound_relax_factor'] = 0
        # opt.options['honor_original_bounds'] = 'yes'
        # opt.options['nlp_scaling_max_gradient'] = 1
        # opt.options['max_iter'] = 20000
        # opt.options['halt_on_ampl_error'] = 'yes'
        #
        # # Write optimization history?
        # opthist = True
        # results = opt.solve(solProg, tee=opthist)
        # solProg.pprint()

        # The multistart optimizer works fine for all objective functions (but spews a bunch of warnings)
        # It's a bit slower that IPOPT, but achieves slightly better results
        # opt = SolverFactory('multistart')
        # results = opt.solve(solProg)

        # The gdpopt optimizer works fine for all objective functions and is fast for all
        # It achieves the same results as multistart
        opt = SolverFactory('gdpopt')
        results = opt.solve(solProg)

        d = self._ExtractData(solProg.d, solProg.N)
        v = self._ExtractData(solProg.v, solProg.N)
        a = self._ExtractData(solProg.a, solProg.N)
        j = self._ExtractData(solProg.j, solProg.N)
        P = self._ExtractData(solProg.p, solProg.N)

        if var[0] == 't':
            dt = solProg.dt[1].value
        else:
            dt = self.dt

        ret = [d, v, a, j, P, dt]
        return ret, solProg
