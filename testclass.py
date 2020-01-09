# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 16:16:20 2020

@author: slauniai
"""
import numpy as np

class A:
    def __init__(self):
        self.a = 1.0
        self.b = 'b'
    def outer(self):
        def inner(self):
            print(self.b)

        print(self.a)
        inner(self)
        self.outer2(55.0)
        
    def outer2(self, arg):
        print('outer2', arg)
        

def modifier(w):
    f = -0.04 * w + 1.38
    ix = (w < 7.0)
    f[ix] = -0.45 + 0.4 * w[ix] - 0.0273 * w[ix] ** 2

    return f