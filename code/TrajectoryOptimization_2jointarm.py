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
from pyomo.environ import cos, sin
import math


class TrajectoryResults:

    def __init__(self, results) -> None:
        self.results = results


class TrajectoryGenerator:
    def __init__(self, dt,
                 theta1i, theta2i, theta1f, theta2f, w1f, w2f, tipVelXf, tipVelYf,
                 tipmass, arm1mass, arm2mass, J1, J2,
                 L1, L2,
                 theta1_bounds, theta2_bounds,
                 w1_bounds, w2_bounds,
                 z1_bounds, z2_bounds,
                 jerk1_bounds, jerk2_bounds,
                 tf=-1) -> None:

        # Initial/final conditions
        self.theta1i = theta1i
        self.theta2i = theta2i
        self.theta1f = theta1f
        self.theta2f = theta2f
        self.w1f = w1f
        self.w2f = w2f
        self.tipVelXf = tipVelXf
        self.tipVelYf = tipVelYf

        # Side constraints
        self.theta1_bounds = theta1_bounds
        self.theta2_bounds = theta2_bounds
        self.w1_bounds = w1_bounds
        self.w2_bounds = w2_bounds
        self.z1_bounds = z1_bounds
        self.z2_bounds = z2_bounds
        self.jerk1_bounds = jerk1_bounds
        self.jerk2_bounds = jerk2_bounds

        # Parameters
        self.tipmass = tipmass
        self.m1 = arm1mass
        self.m2 = arm2mass
        self.J1 = J1
        self.J2 = J2
        self.L1 = L1
        self.L2 = L2
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

                if norm == 'peak':
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

                elif norm == 'rms':
                    objective += weight * sum(getattr(model, var)[i] ** 2 for i in model.N)

                elif norm == 'abs':
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
        # return [a*v for a, v in zip(model.a, model.v)]
        # TrajectoryGenerator._PowerCalculation(model.a, model.v)

    @staticmethod
    def _TotalEnergyExpression(model, i):
        val = sum_product(model.a, model.v)*model.mass*model.dt
        return val

    @staticmethod
    def _Tip1PosXExpression(model, i):
        return model.L1*cos(model.theta1[i])

    @staticmethod
    def _Tip2PosXExpression(model, i):
        return model.L1*cos(model.theta1[i]) + model.L2*cos(model.theta1[i]+model.theta2[i])

    @staticmethod
    def _Tip1PosYExpression(model, i):
        return model.L1*sin(model.theta1[i])

    @staticmethod
    def _Tip2PosYExpression(model, i):
        return model.L1*sin(model.theta1[i]) + model.L2*sin(model.theta1[i]+model.theta2[i])

    @staticmethod
    def _Tip1VelXExpression(model, i):
        return -model.L1*sin(model.theta1[i])*model.w1[i]

    @staticmethod
    def _Tip2VelXExpression(model, i):
        return -model.L1*sin(model.theta1[i])*model.w1[i] - \
               model.L2*sin(model.theta1[i]+model.theta2[i])*(model.w1[i]+model.w2[i])

    @staticmethod
    def _Tip1VelYExpression(model, i):
        return model.L1*cos(model.theta1[i])*model.w1[i]

    @staticmethod
    def _Tip2VelYExpression(model, i):
        return model.L1*cos(model.theta1[i])*model.w1[i] + \
               model.L2*cos(model.theta1[i]+model.theta2[i])*(model.w1[i]+model.w2[i])

    @staticmethod
    def _Torque1Expression(model, i):
        a = model.m1*(model.L1/2)**2 + model.m2*model.L1**2 + model.J1 + \
            0.5*model.m2*model.L1*model.L2*cos(model.theta2[i]-model.theta1[i])
        b = 0.5*model.m2*model.L1*model.L2*cos(model.theta2[i]-model.theta1[i]) + model.m2*(model.L2/2)**2 + model.J2
        c = 0.5*model.m2*model.L1*model.L2*model.w1[i]**2*sin(model.theta2[i]-model.theta1[i])
        d = 0.5*model.m2*model.L1*model.L2*model.w2[i]**2*sin(model.theta2[i]-model.theta1[i])

        return a*model.z1[i] + b*model.z2[i] + c - d

    @staticmethod
    def _Torque2Expression(model, i):
        a = 0.5*model.m2*model.L1*model.L2*cos(model.theta2[i]-model.theta1[i])
        b = model.m2*(model.L2/2)**2 + model.J2
        c = 0.5*model.m2*model.L1*model.L2*model.w1[i]**2*sin(model.theta2[i]-model.theta1[i])

        return a*model.z1[i] + b*model.z2[i] + c

    @staticmethod
    def _Power1Expression(model, i):
        return model.t1[i]*model.w1[i]

    @staticmethod
    def _Power2Expression(model, i):
        return model.t2[i]*model.w2[i]

    def OptimizationInitialize(self):
        # This initializes the model and sets as many of the system parameters and variables as possible.
        N = int(self.tf / self.dt)
        dt = self.dt
        optProg = ConcreteModel()

        # Initialize the RangeSet that defines the number
        optProg.N = RangeSet(0, math.ceil(self.tf / self.dt))
        optProg.Nf = RangeSet(1, math.ceil(self.tf / self.dt))

        # Parameters
        optProg.dt = Param(initialize=dt)
        optProg.tipmass = Param(initialize=self.tipmass)
        optProg.m1 = Param(initialize=self.m1)
        optProg.m2 = Param(initialize=self.m2)
        optProg.J1 = Param(initialize=self.J1)
        optProg.J2 = Param(initialize=self.J2)
        optProg.L1 = Param(initialize=self.L1)
        optProg.L2 = Param(initialize=self.L2)

        # Variables
        optProg.theta1 = Var(optProg.N, domain=Reals, bounds=self._VariableBounds(self.theta1_bounds[0], self.theta1_bounds[1]))
        optProg.w1 = Var(optProg.N, domain=Reals, bounds=self._VariableBounds(self.w1_bounds[0], self.w1_bounds[1]))
        optProg.z1 = Var(optProg.N, domain=Reals, bounds=self._VariableBounds(self.z1_bounds[0], self.z1_bounds[1]))
        optProg.j1 = Var(optProg.N, domain=Reals, bounds=self._VariableBounds(self.jerk1_bounds[0], self.jerk1_bounds[1]))
        optProg.theta2 = Var(optProg.N, domain=Reals, bounds=self._VariableBounds(self.theta2_bounds[0], self.theta2_bounds[1]))
        optProg.w2 = Var(optProg.N, domain=Reals, bounds=self._VariableBounds(self.w2_bounds[0], self.w2_bounds[1]))
        optProg.z2 = Var(optProg.N, domain=Reals, bounds=self._VariableBounds(self.z2_bounds[0], self.z2_bounds[1]))
        optProg.j2 = Var(optProg.N, domain=Reals, bounds=self._VariableBounds(self.jerk2_bounds[0], self.jerk2_bounds[1]))

        # Expressions (dependent variables)
        optProg.t1 = Expression(optProg.N, rule=self._Torque1Expression)
        optProg.t2 = Expression(optProg.N, rule=self._Torque2Expression)
        optProg.p1 = Expression(optProg.N, rule=self._Power1Expression)
        optProg.p2 = Expression(optProg.N, rule=self._Power2Expression)
        optProg.tip1PosX = Expression(optProg.N, rule=self._Tip1PosXExpression)
        optProg.tip1PosY = Expression(optProg.N, rule=self._Tip1PosYExpression)
        optProg.tip1VelX = Expression(optProg.N, rule=self._Tip1VelXExpression)
        optProg.tip1VelY = Expression(optProg.N, rule=self._Tip1VelYExpression)
        optProg.tip2PosX = Expression(optProg.N, rule=self._Tip2PosXExpression)
        optProg.tip2PosY = Expression(optProg.N, rule=self._Tip2PosYExpression)
        optProg.tip2VelX = Expression(optProg.N, rule=self._Tip2VelXExpression)
        optProg.tip2VelY = Expression(optProg.N, rule=self._Tip2VelYExpression)

        # Initial conditions
        optProg.theta1_0 = Constraint([0], rule=self._EqualityConstr(optProg.theta1, [self.theta1i]))
        optProg.w1_0 = Constraint([0], rule=self._EqualityConstr(optProg.w1, [0]))
        optProg.z1_0 = Constraint([0], rule=self._EqualityConstr(optProg.z1, [0]))
        optProg.theta2_0 = Constraint([0], rule=self._EqualityConstr(optProg.theta2, [self.theta2i]))
        optProg.w2_0 = Constraint([0], rule=self._EqualityConstr(optProg.w2, [0]))
        optProg.z2_0 = Constraint([0], rule=self._EqualityConstr(optProg.z2, [0]))

        # Final conditions
        # optProg.theta1_N = Constraint([N], rule=self._EqualityConstr(optProg.theta1, [self.theta1f], N))
        # optProg.w1_N = Constraint([N], rule=self._EqualityConstr(optProg.w1, [self.w1f], N))
        # optProg.theta2_N = Constraint([N], rule=self._EqualityConstr(optProg.theta2, [self.theta2f], N))
        # optProg.w2_N = Constraint([N], rule=self._EqualityConstr(optProg.w2, [self.w2f], N))
        optProg.tipVelX_N = Constraint([N], rule=self._EqualityConstr(optProg.tip2VelX, [self.tipVelXf], N))
        optProg.tipVelY_N = Constraint([N], rule=self._EqualityConstr(optProg.tip2VelY, [self.tipVelYf], N))

        # Dynamic constraints to ensure dynamic compatibility between theta, w, z, j
        optProg.w1_k = Constraint(optProg.Nf, rule=self._DynamicConstr(optProg.theta1, optProg.w1, self.dt))
        optProg.z1_k = Constraint(optProg.Nf, rule=self._DynamicConstr(optProg.w1, optProg.z1, self.dt))
        optProg.j1_k = Constraint(optProg.Nf, rule=self._DynamicConstr(optProg.z1, optProg.j1, self.dt))
        optProg.w2_k = Constraint(optProg.Nf, rule=self._DynamicConstr(optProg.theta2, optProg.w2, self.dt))
        optProg.z2_k = Constraint(optProg.Nf, rule=self._DynamicConstr(optProg.w2, optProg.z2, self.dt))
        optProg.j2_k = Constraint(optProg.Nf, rule=self._DynamicConstr(optProg.z2, optProg.j2, self.dt))

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

        theta1 = self._ExtractData(solProg.theta1, solProg.N)
        theta2 = self._ExtractData(solProg.theta2, solProg.N)
        w1 = self._ExtractData(solProg.w1, solProg.N)
        w2 = self._ExtractData(solProg.w2, solProg.N)
        z1 = self._ExtractData(solProg.z1, solProg.N)
        z2 = self._ExtractData(solProg.z2, solProg.N)
        j1 = self._ExtractData(solProg.j1, solProg.N)
        j2 = self._ExtractData(solProg.j2, solProg.N)
        t1 = self._ExtractData(solProg.t1, solProg.N)
        t2 = self._ExtractData(solProg.t2, solProg.N)
        p1 = self._ExtractData(solProg.p1, solProg.N)
        p2 = self._ExtractData(solProg.p2, solProg.N)
        tip1PosX = self._ExtractData(solProg.tip1PosX, solProg.N)
        tip1PosY = self._ExtractData(solProg.tip1PosY, solProg.N)
        tip1VelX = self._ExtractData(solProg.tip1VelX, solProg.N)
        tip1VelY = self._ExtractData(solProg.tip1VelY, solProg.N)
        tip2PosX = self._ExtractData(solProg.tip2PosX, solProg.N)
        tip2PosY = self._ExtractData(solProg.tip2PosY, solProg.N)
        tip2VelX = self._ExtractData(solProg.tip2VelX, solProg.N)
        tip2VelY = self._ExtractData(solProg.tip2VelY, solProg.N)

        if var[0] == 't':
            dt = solProg.dt[1].value
        else:
            dt = self.dt

        ret = [dt, theta1, theta2, w1, w2, z1, z2, j1, j2, t1, t2, p1, p2,
               tip1PosX, tip1PosY, tip1VelX, tip1VelY, tip2PosX, tip2PosY, tip2VelX, tip2VelY]
        return ret, solProg
