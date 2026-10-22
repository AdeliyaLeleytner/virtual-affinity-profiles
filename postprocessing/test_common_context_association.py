import unittest
import numpy as np
from common_context_association import estimate_pair


class CommonAssociation(unittest.TestCase):
    def test_macro_uses_intersection_not_separate_denominators(self):
        a={'qy':np.array([[0.],[1.]]),'yy':np.ones((2,1)),'qq':np.array([[0.],[1.]])}
        b={'qy':np.array([[-1.],[1.]]),'yy':np.ones((2,1)),'qq':np.ones((2,1))}
        arm,base,common=estimate_pair(a,b,np.ones((2,1)),'actual_assay_pair_equal')
        self.assertEqual((arm[0],base[0]),(1,1));self.assertEqual(common.sum(),1)

    def test_intersection_is_separate_for_each_draw(self):
        a={'qy':np.array([[0.,1.],[1.,0.]]),'yy':np.ones((2,2)),'qq':np.array([[0.,1.],[1.,0.]])}
        b={'qy':np.array([[-1.,-1.],[1.,1.]]),'yy':np.ones((2,2)),'qq':np.ones((2,2))}
        arm,base,_=estimate_pair(a,b,np.ones((2,2)),'actual_assay_pair_equal')
        np.testing.assert_array_equal(arm,[1,1]);np.testing.assert_array_equal(base,[1,-1])

    def test_raw_keeps_same_rectangle_pool(self):
        a={'qy':np.array([[0.],[1.]]),'yy':np.ones((2,1)),'qq':np.array([[0.],[1.]])}
        b={'qy':np.array([[-1.],[1.]]),'yy':np.ones((2,1)),'qq':np.ones((2,1))}
        arm,base,_=estimate_pair(a,b,np.ones((2,1)),'raw_rectangle_equal')
        self.assertAlmostEqual(arm[0],1/np.sqrt(2));self.assertEqual(base[0],0)


if __name__=='__main__':unittest.main()
