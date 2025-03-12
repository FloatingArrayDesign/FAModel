# -*- coding: utf-8 -*-
"""
Created on Wed May 15 12:55:16 2024

@author: lsirkis
"""
from famodel.famodel_base import Node




class Substation(Node):
    
    def __init__(self, dd, id):
        
        Node.__init__(self, id)  # initialize Edge base class
        
        self.dd = dd
        
        # cost dictionary
        self.cost = {}
        # dictionary of failure probability
        self.failure_probability = {}
        
        # further functionality to be added later
        
        