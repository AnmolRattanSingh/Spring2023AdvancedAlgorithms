import picos as pic
from picos import RealVariable
from copy import deepcopy
from heapq import *
import heapq as hq
import numpy as np
import itertools
import math

counter = itertools.count()

class BBTreeNode():
    """
    Creates and handles a BBTreeNode object that can branch and bound
    to determine the optimal result and corresponding best variable
    values.

    Attributes:
        vars (list of picos RealVariable objects): variables in the
            problem.
        constraints (list of constraints): list of problem constraints.
            ex: [z == x+y, -5*x+4*y <= 0, 6*x+2*y <= 17, x>=0, y>=0].
        objective (picos RealVariable object): variable that is being
            maximized.
        prob (picos Problem object): problem created by buildProblem
            using constraints, vars, and objective.
    """

    def __init__(self, vars = [], constraints = [], objective='', prob=None):
        """
        Initializes BBTreeNode.
        """
        self.vars = vars
        self.constraints = constraints
        self.objective = objective
        self.prob = prob

    def __deepcopy__(self, memodict={}):
        """
        Deepcopies the picos problem.
        This overrides the system's deepcopy method bc it doesn't work
        on classes by itself.

        Returns:
            (BBTreeNode object): copy of BBTreeNode.
        """
        newprob = pic.Problem.clone(self.prob)
        return BBTreeNode(self.vars, newprob.constraints, self.objective, newprob)

    def buildProblem(self):
        """
        Builds the initial Picos problem.

        Returns:
            self.prob (picos Problem object): problem created from
                constraints, objective, and vars.
        """
        prob=pic.Problem()
        prob.add_list_of_constraints(self.constraints)
        prob.set_objective('max', self.objective)
        self.prob = prob
        return self.prob

    def is_integral(self):
        """
        Checks if all variables (excluding the one we're maxing) are
        integers.

        Returns:
            (bool): returns True if all variables (excluding the one
                we're maximizing) are integers, otherwise False.
        """
        for v in self.vars[:-1]:
            if v.value == None or abs(round(v.value) - float(v.value)) > 1e-4:
                return False
        return True

    def branch_floor(self, branch_var):
        """
        Makes a child where xi <= floor(xi).

        Args:
            branch_var (float): variable to branch on.

        Returns:
            n1 (BBTreeNode object): child where xi <= floor(xi).
        """
        n1 = deepcopy(self)
        # add in the new binary constraint
        n1.prob.add_constraint(branch_var <= math.floor(branch_var.value))
        return n1

    def branch_ceil(self, branch_var):
        """
        Makes a child where xi >= ceiling(xi).

        Args:
            branch_var (float): variable to branch on.

        Returns:
            n2 (BBTreeNode object): child where xi >= ceiling(xi).
        """
        n2 = deepcopy(self)
        # add in the new binary constraint
        n2.prob.add_constraint(branch_var >= math.ceil(branch_var.value))
        return n2

    def bbsolve(self):
        """
        Uses the branch and bound method to solve an integer program.

        Returns:
            bestres (float): value of the maximized objective function.
            bestnode_vars (list of floats): list of variables that
                create bestres.
        """
        # these lines build up the initial problem and add it to a heap
        root = self
        res = root.buildProblem().solve(solver='cvxopt')
        heap = [(res, next(counter), root)]
        # set bestres to an arbitrary small initial best objective value
        bestres = -1e20
        # initialize bestnode_vars to the root vars
        bestnode_vars = root.vars

        while heap:
            # get the node with the best value so far
            # _, _, cur = hq._heappop_max(heap)
            _, _, cur = hq.heappop(heap)
            # print hq
            print([x.value for x in cur.vars])

            # check if the current node is integral
            if cur.is_integral():
                # if it is, update the best result if necessary
                if cur.objective.value > bestres:
                    bestres = cur.objective.value
                    bestnode_vars = [v.value for v in cur.vars]
            else:
                # if the node is not integral, we need to branch
                # get the next branching variable
                branch_var = None
                for v in cur.vars[:-1]:
                    if abs(round(v.value) - float(v.value)) > 1e-4:
                        branch_var = v
                        break

                # branch on the next branching variable
                n1 = cur.branch_floor(branch_var)
                n2 = cur.branch_ceil(branch_var)
                # check if the first branch is feasible
                try:
                    if n1.prob.solve(solver='cvxopt') != 'infeasible':
                        hq.heappush(heap, (n1.objective.value, next(counter), n1))
                except:
                    pass

                # check if the second branch is feasible
                try:
                    if n2.prob.solve(solver='cvxopt') != 'infeasible':
                        hq.heappush(heap, (n2.objective.value, next(counter), n2))
                except:
                    pass
                # hq._heapify_max(heap)
                print("---", [x[0] for x in list(heap)])

        return bestres, bestnode_vars

if False:
    x = RealVariable("x")
    y = RealVariable("y")
    z = RealVariable("z")
    a= RealVariable("a")
    b= RealVariable("b")
    c = RealVariable("c")
    vars = [x, y, a, b, c,  z]
    objective = z
    constraints = [z == 15*x + 20 *y + 18* a + 13 * b + 12* c, 18*x+10*y+21*a+11*b+11*c <= 50, x>= 0, y >= 0, a >= 0, b >= 0, c >= 0, x <= 1, y <= 1, a <= 1, b <= 1, c <= 1]
    root = BBTreeNode(constraints = constraints, objective = objective, vars = vars)
    res, sol_vars = root.bbsolve()
else:
    x = RealVariable("x")
    y = RealVariable("y")
    z = RealVariable("z")
    vars = [x, y, z]
    objective = z
    constraints = [x+y <= 7, 12*x+ 5*y <= 60, x >= 0, y >=0, z == 80 * x + 45 * y]
    root = BBTreeNode(constraints = constraints, objective = objective, vars = vars)
    res, sol_vars = root.bbsolve()
print(sol_vars)
