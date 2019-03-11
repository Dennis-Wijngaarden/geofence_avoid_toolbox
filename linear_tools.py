#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 15:34:22 2019

@author: dennis
"""

import numpy as np

def line_eq_from_points(point0, point1):
    # a x_0 + b y_0 + c = 0
    # a x_1 + b y_1 + c = 0
    # a x_0 + b y_0 = a x_1 + b y_1 
    try:
        a_over_b = (point1[1] - point0[1]) / (point0[0] - point1[0])
        if abs(a_over_b) <= 1.:
            b = 1.
            a = a_over_b
            
        else:
            b_over_a = (point1[0] - point0[0]) / (point0[1] - point1[1])
            a = 1.
            b = b_over_a
    except RuntimeWarning:
        b_over_a = (point1[0] - point0[0]) / (point0[1] - point1[1])
        a = 1.
        b = b_over_a
    
    c = -a * point0[0] - b * point0[1]
    return (a, b, c)