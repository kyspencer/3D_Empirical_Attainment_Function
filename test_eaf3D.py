''' test_eaf3D.py

    This python scripts tests methods in eaf3D.py.
        - insert()
        - adjustBalances()
        - case1()
        - case2()
        - case3()

'''

import unittest
import avltree_eaf3d as bst
import eaf3D
import numpy as np
from operator import attrgetter
from stack import Stack


class EAF3DTransformTests(unittest.TestCase):

    def setUp(self):
        folder = 'example/'
        sets = eaf3D.retrieve_input_sequences(folder)
        self.eaf_maker = eaf3D.EAF_3D(sets)

    def test_ensure_setUp_worked(self):
        self.assertEqual(self.eaf_maker.lstar[0].root.point.x, 3)
        self.assertAlmostEqual(self.eaf_maker.lstar[0].root.point.y, 22.59977497)

    def test_transform(self):
        self.eaf_maker.transform()
        for i in range(self.eaf_maker.n):
            self.check_node_balances(self.eaf_maker.lstar[i])

    def test_find_attainment_point(self):
        # p will be point at (4, 15.42, 6131.4) from set 3
        p = self.eaf_maker.qstack.pop()
        j = p.input_set()
        # q will be point at (-inf, inf, -inf)
        q = self.eaf_maker.xstar[j].floor_x(p)
        # Here, t=0, so only one while loop encountering r at (3, 22.6, 6131),
        # resulting in new point s[0] at (4, 22.6, 6131)
        t, tmin = self.eaf_maker.tmax, 0
        s, tmin = self.eaf_maker.find_attainment_point(p, q, t, tmin)
        self.assertEqual(s[0].x, 4)
        self.assertAlmostEqual(s[0].y, 22.59977497)
        self.assertEqual(tmin, 0)

    def test_compare_p_to_surfaces(self):
        # p will be point at (4, 15.42, 6131.4) from set 3
        p = self.eaf_maker.qstack.pop()
        j = p.input_set()
        # q will be point at (-inf, inf, -inf)
        q = self.eaf_maker.xstar[j].floor_x(p)
        # Here, t=0, so only one while loop encountering r at (3, 22.6, 6131),
        # resulting in new point s[0] at (4, 22.6, 6131)
        t, tmin = self.eaf_maker.tmax, 0
        s, tmin = self.eaf_maker.find_attainment_point(p, q, t, tmin)
        s = self.eaf_maker.compare_p_to_surfaces(s, p, q, j, tmin)
        self.assertEqual(len(s), self.eaf_maker.n)

    def test_lstar_higher_x(self):
        # p will be point at (4, 15.42, 6131.4) from set 3
        p = self.eaf_maker.qstack.pop()
        j = p.input_set()
        # q will be point at (-inf, inf, -inf)
        q = self.eaf_maker.xstar[j].floor_x(p)
        newq = self.eaf_maker.xstar[j].higher_x(q)
        self.assertEqual(newq.x, 10E10)

    def test_submit_to_lstar(self):
        st = self.eaf_maker.lstar[0].root.point
        self.eaf_maker.submit_to_lstar(st, 1)
        (pivot, theStack, parent, found) = self.eaf_maker.lstar[1].search(st)
        self.assertTrue(found)

    def test_submit_to_lstar_wdelete(self):
        # This tests is EAF3D correctly deletes nodes from trees
        # Loop 1
        p = self.eaf_maker.qstack.pop()
        j = p.input_set()
        q = self.eaf_maker.xstar[j].floor_x(p)
        if p.y < q.y:
            t, tmin = self.eaf_maker.tmax, 0
            s, tmin = self.eaf_maker.find_attainment_point(p, q, t, tmin)
            s = self.eaf_maker.compare_p_to_surfaces(s, p, q, j, tmin)
            self.eaf_maker.submit_points_lstar(s, p, q, tmin)
            self.eaf_maker.submit_to_xstar(p, j)
        if j not in self.eaf_maker.a_tracker:
            self.eaf_maker.a_tracker.append(j)
            self.eaf_maker.tmax = min(self.eaf_maker.tmax + 1, self.eaf_maker.n - 2)
        # Loop 1.5
        p = self.eaf_maker.qstack.pop()
        j = p.input_set()
        q = self.eaf_maker.xstar[j].floor_x(p)
        if p.y < q.y:
            t, tmin = self.eaf_maker.tmax, 0
            s, tmin = self.eaf_maker.find_attainment_point(p, q, t, tmin)
            s = self.eaf_maker.compare_p_to_surfaces(s, p, q, j, tmin)

    def check_node_balances(self, tree):
        # This module performs a check of all the balances in the tree
        leaves = [tree.root]
        while any(leaves):
            for f in range(len(leaves)):
                if leaves[f]:
                    correct_balance = tree.recalculate_balance(leaves[f])
                    self.assertEqual(leaves[f].balance, correct_balance)
            leaves = tree.next_tree_row(leaves)



if __name__ == '__main__':
    unittest.main()
