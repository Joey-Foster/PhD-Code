import numpy as np
from TwoDimStaticAdvDiffFESolver import *
from TwoDimTimeEvolvedAdvDiffFESolver import mass
import pytest

def test_local2globalCoords():
    
    default = {
                "xe": np.array([[0, 1, 0],
                                [0, 0, 1]]),
                "ans": np.array([0.5, 0.5])
              }
    
    translated = {
                "xe": np.array([[1, 2, 1],
                                [0, 0, 1]]),
                "ans": np.array([1.5, 0.5])
                 }
    
    scaled = {
                "xe": np.array([[0, 2, 0],
                                [0, 0, 2]]),
                "ans": np.array([1, 1])
             }
    
    rotated = {
                "xe": np.array([[1, 0, 1],
                                [1, 1, 0]]),
                "ans": np.array([0.5, 0.5])
              }
    
    xi = np.array([0.5, 0.5])
    
    for t in [default, translated, scaled, rotated]:
        assert np.allclose(local2globalCoords(t["xe"], xi),
                           t["ans"]), f"element\n {t['xe']} is broken"
  
def test_jacobian():
      
    default = {
                "xe": np.array([[0, 1, 0],
                                [0, 0, 1]]),
                "ans": np.array([[1, 0],
                                 [0, 1]])
              }
    
    translated = {
                "xe": np.array([[1, 2, 1],
                                [0, 0, 1]]),
                "ans": np.array([[1, 0],
                                 [0, 1]])
                 }
    
    scaled = {
                "xe": np.array([[0, 2, 0],
                                [0, 0, 2]]),
                "ans": np.array([[2, 0],
                                 [0, 2]])
             }
    
    rotated = {
                "xe": np.array([[1, 0, 1],
                                [1, 1, 0]]),
                "ans": np.array([[-1, 0],
                                 [0, -1]])
              }
    
    for t in [default, translated, scaled, rotated]:
        assert np.allclose(jacobian(t["xe"]),
                           t["ans"]), f"element\n {t['xe']} is broken"
        
def test_globalShapeFunctionDerivatives():
    
   default = {
               "xe": np.array([[0, 1, 0],
                               [0, 0, 1]]),
               "ans": np.array([[-1, 1, 0],
                                [-1, 0, 1]])
             }
   
   translated = {
               "xe": np.array([[1, 2, 1],
                               [0, 0, 1]]),
               "ans": np.array([[-1, 1, 0],
                                [-1, 0, 1]])
                }
   
   scaled = {
               "xe": np.array([[0, 2, 0],
                               [0, 0, 2]]),
               "ans": np.array([[-0.5, 0.5, 0],
                                [-0.5, 0, 0.5]])
            }
   
   rotated = {
               "xe": np.array([[1, 0, 1],
                               [1, 1, 0]]),
               "ans": np.array([[1, -1, 0],
                                [1, 0, -1]])
             }
   
   for t in [default, translated, scaled, rotated]:
       assert np.allclose(globalShapeFunctionDerivatives(t["xe"]),
                          t["ans"]), f"element\n {t['xe']} is broken"
    
def test_localQuadtrature():
    
    constant = {
             "psi": lambda x: 1,
             "ans": 0.5 
             } 
    
    linearx = {
             "psi": lambda x: 6*x[0],
             "ans": 1,
             }
    
    lineary = {
             "psi": lambda x: x[1],
             "ans": 1/6,
             }
    
    product = {
             "psi": lambda x: x[0]*x[1],
             "ans": 1/24,
             }
    for t in [constant, linearx, lineary, product]:
        assert np.allclose(localQuadrature(t["psi"]),
                           t["ans"]), f"function with answer\n {t['ans']} is broken"
        
def test_globalQuadrature():   
    
    translated_linear = {
             "xe": np.array([[1, 2, 1],
                             [0, 0, 1]]),
             "phi": lambda x: 3*x[0],
             "ans": 2,
                        }
    
    translated_product = {
             "xe": np.array([[1, 2, 1],
                             [0, 0, 1]]),
             "phi": lambda x: x[0]*x[1],
             "ans": 5/24,
                         }
    
    scaled_linear = {
             "xe": np.array([[0, 2, 0],
                             [0, 0, 2]]),
             "phi": lambda x: 3*x[0],
             "ans": 4,
                     }
    
    scaled_product = {
             "xe": np.array([[0, 2, 0],
                             [0, 0, 2]]),
             "phi": lambda x: x[0]*x[1],
             "ans": 2/3,
                     }
    
    rotated_linear = {
             "xe": np.array([[1, 0, 1],
                             [1, 1, 0]]),
             "phi": lambda x: 3*x[0],
             "ans": 1,
                     }
    
    rotated_product = {
             "xe": np.array([[1, 0, 1],
                             [1, 1, 0]]),
             "phi": lambda x: x[0]*x[1],
             "ans": 5/24,
                     }

    for t in [translated_linear, translated_product, scaled_linear, scaled_product,
              rotated_linear, rotated_product]:
        assert np.allclose(globalQuadrature(t["xe"],t["phi"]), 
                           t["ans"]), f"element\n {t['xe']} with answer\n {t['ans']} is broken"

def test_diffusion_stiffness():

    default = {
                "xe": np.array([[0, 1, 0],
                                [0, 0, 1]]),
                "ans": np.array([[1, -0.5, -0.5],
                                 [-0.5, 0.5, 0],
                                 [-0.5, 0, 0.5]])
              }
    
    translated = {
                "xe": np.array([[1, 2, 1],
                                [0, 0, 1]]),
                "ans": np.array([[1, -0.5, -0.5],
                                 [-0.5, 0.5, 0],
                                 [-0.5, 0, 0.5]])
                 }
    
    scaled = {
                "xe": np.array([[0, 2, 0],
                                [0, 0, 2]]),
                "ans": np.array([[1, -0.5, -0.5],
                                 [-0.5, 0.5, 0],
                                 [-0.5, 0, 0.5]])
            }
    rotated = {
                "xe": np.array([[1, 0, 1],
                                [1, 1, 0]]),
                "ans": np.array([[1, -0.5, -0.5],
                                 [-0.5, 0.5, 0],
                                 [-0.5, 0, 0.5]])
              }
    
    for t in [default, translated, scaled, rotated]:
        assert np.allclose(diffusion_stiffness(t["xe"]),
                           t["ans"]), f"element\n {t['xe']} is broken"
        

def test_advection_stiffness():
    
    default_east = {
                    "xe": np.array([[0, 1, 0],
                                    [0, 0, 1]]),
                    "u": np.array([1, 0]),
                    "ans": 1/6 * np.array([[-1,  1,  0],
                                           [-1,  1,  0],
                                           [-1,  1,  0]])
                   }
    
    default_north = {
                    "xe": np.array([[0, 1, 0],
                                    [0, 0, 1]]),
                    "u": np.array([0, 1]),
                    "ans": 1/6 * np.array([[-1,  0,  1],
                                           [-1,  0,  1],
                                           [-1,  0,  1]])
                   }
    
    default_north_east = {
                    "xe": np.array([[0, 1, 0],
                                    [0, 0, 1]]),
                    "u": np.array([1, 1]),
                    "ans": 1/6 * np.array([[-2,  1,  1],
                                           [-2,  1,  1],
                                           [-2,  1,  1]])
                        }
    
    translated = {
                "xe": np.array([[1, 2, 1],
                                [0, 0, 1]]),
                "u": np.array([0, 1]),
                "ans": 1/6 * np.array([[-1,  0,  1],
                                       [-1,  0,  1],
                                       [-1,  0,  1]])
                 }
    
    scaled = {
                "xe": np.array([[0, 2, 0],
                                [0, 0, 2]]),
                "u": np.array([0, 1]),
                "ans": 1/3 * np.array([[-1,  0,  1],
                                       [-1,  0,  1],
                                       [-1,  0,  1]])
            }
    rotated = {
                "xe": np.array([[1, 0, 1],
                                [1, 1, 0]]),
                "u": np.array([0, 1]),
                "ans": 1/6 * np.array([[1,  0, -1],
                                       [1,  0, -1],
                                       [1,  0, -1]])
              }
    
    for t in [default_east, default_north, default_north_east, translated, 
              scaled, rotated]:
        assert np.allclose(advection_stiffness(t["xe"],t["u"]),
                           t["ans"]), f"element\n {t['xe']} with velocity\n {t['u']} is broken"
        
def test_global2localcoords():
    
    default = {
                "xe": np.array([[0, 1, 0],
                                [0, 0, 1]]),
                "x": np.array([0.5, 0.5])
              }
    
    translated = {
                "xe": np.array([[1, 2, 1],
                                [0, 0, 1]]),
                "x": np.array([1.5, 0.5])
                 }
    
    scaled = {
                "xe": np.array([[0, 2, 0],
                                [0, 0, 2]]),
                "x": np.array([1, 1])
             }
    
    rotated = {
                "xe": np.array([[1, 0, 1],
                                [1, 1, 0]]),
                "x": np.array([0.5, 0.5])
              }
    
    ans = np.array([0.5, 0.5])
    
    for t in [default, translated, scaled, rotated]:
        assert np.allclose(global2localCoords(t["xe"], t["x"]), 
                           ans), f"element\n {t['xe']} is broken"
        
def test_force():
    
    default_const = {
             "xe": np.array([[0, 1, 0],
                             [0, 0, 1]]),
             "S": lambda x: 1,
             "ans": 1/6 * np.array([1, 1, 1]),
                        }
    
    default_linear = {
             "xe": np.array([[0, 1, 0],
                             [0, 0, 1]]),
             "S": lambda x: x[0],
             "ans": 1/24 * np.array([1, 2, 1])
                        }

       
    translated_const = {
             "xe": np.array([[1, 2, 1],
                             [0, 0, 1]]),
             "S": lambda x: 1,
             "ans": 1/6 * np.array([1, 1, 1]),
                        }
    
    translated_linear = {
             "xe": np.array([[1, 2, 1],
                             [0, 0, 1]]),
             "S": lambda x: x[1],
             "ans": 1/24 * np.array([1, 1, 2])
                         }
    
    scaled_const = {
             "xe": np.array([[0, 2, 0],
                             [0, 0, 2]]),
             "S": lambda x: 1,
             "ans": 2/3 * np.array([1, 1, 1]),
                     }
    
    scaled_linear = {
             "xe": np.array([[0, 2, 0],
                             [0, 0, 2]]),
             "S": lambda x: x[0],
             "ans": 1/3 * np.array([1, 2, 1]),
                     }
    
    rotated_const = {
             "xe": np.array([[1, 0, 1],
                             [1, 1, 0]]),
             "S": lambda x: 1,
             "ans": 1/6 * np.array([1, 1, 1]),
                     }
    
    rotated_linear = {
             "xe": np.array([[1, 0, 1],
                             [1, 1, 0]]),
             "S": lambda x: x[1],
             "ans": 1/4 * np.array([1/2, 1/2, 1/3]),
                     }

    for t in [default_const, default_linear, translated_const, translated_linear, 
              scaled_const, scaled_linear, rotated_const, rotated_linear]:
        assert np.allclose(force(t["xe"],t["S"]),
                           t["ans"]), f"element\n {t['xe']} with answer\n {t['ans']} is broken"

def test_mass():

    default = {
                "xe": np.array([[0, 1, 0],
                                [0, 0, 1]]),
                "ans": 1/24 * np.array([[2, 1, 1],
                                        [1, 2, 1],
                                        [1, 1, 2]])
              }
    
    translated = {
                "xe": np.array([[1, 2, 1],
                                [0, 0, 1]]),
                "ans": 1/24 * np.array([[2, 1, 1],
                                        [1, 2, 1],
                                        [1, 1, 2]])
                 }
    
    scaled = {
                "xe": np.array([[0, 2, 0],
                                [0, 0, 2]]),
                "ans": 1/6 * np.array([[2, 1, 1],
                                        [1, 2, 1],
                                        [1, 1, 2]])
            }
    rotated = {
                "xe": np.array([[1, 0, 1],
                                [1, 1, 0]]),
                "ans": 1/24 * np.array([[2, 1, 1],
                                        [1, 2, 1],
                                        [1, 1, 2]])
              }
    
    for t in [default, translated, scaled, rotated]:
        assert np.allclose(mass(t["xe"]),
                           t["ans"]), f"element\n {t['xe']} is broken"