# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 00:39:45 2021

@author: Mihir Gadgi
"""
from manim import *
import Dipole_Sim

class Field(ThreeDScene):
    def construct(self):
        #v = ArrowVectorField(lambda x: Dipole_Sim.B(x,0), delta_x=2, delta_y=2)
        v = VGroup(
            *[ Vector(Dipole_Sim.B(x,0)).move_to(x) for x in
             np.transpose(np.array(np.meshgrid(np.arange(-1,1,0.1),np.arange(-1,1,0.1),np.arange(-1,1,0.1))),
                          [1,2,3,0]).reshape(-1,3)
            ]
            )
    
        self.add(v)
        #self.move_camera(45DEGREES, 45DEGREES, 4,)
        #self.play(self.camera.animate.scale(0.5))

%manim -ql Field