#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 14:06:06 2021

@author: arun
"""
import DigitizeSpline as DS

Fig=DS.DigitizeSpline("cut_spline.txt")
Fig.distribute_nodes(100)
injectionlocations=[[0,0]]
Fig.set_injection_points(injectionlocations)
# Fig.PlotBezier()
# Fig.PlotNodes()

